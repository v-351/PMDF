#include "../include/mpi-wrapper.hpp"

#include <algorithm>

/** Конструктор. 
 * Инициализирует MPI. Определяет количество процессов и номер данного.
 */
MPIw::MPIw(int* argc, char*** argv){
    MPI_Init(argc, argv);
    MPI_Comm_size( MPI_COMM_WORLD , &size);
    MPI_Comm_rank( MPI_COMM_WORLD , &rank);
}

/** Деструктор. 
 * Закрывает MPI-сессию.
 */
MPIw::~MPIw(){ MPI_Finalize(); }

int MPIw::getRank() const{
    return rank;
}

int MPIw::getSize() const{
    return size;
}

/// Отправка частиц на другой процессор. Передаются итератор начала, конца, и адрес.
void MPIw::send(std::vector<Particle>::iterator begin,std::vector<Particle>::iterator end, int address){
    ex.clear(); // Очистка буфера
    // Цикл по всем элементам от begin до end
    for(begin; begin < end; begin++){
        begin->serialize(ex); // Запись свойств частицы одним блоком в числовой массив.
    }
    // Отправка числовго массива на другой процессор.
    MPI_Send( &ex[0] , ex.size() , MPI_DOUBLE , address , 0 , MPI_COMM_WORLD);
}

void MPIw::recv(std::vector<Particle>& particles, int address, ParticleSubstance* substance){
    // Получение информации о посылке
    MPI_Probe( address , MPI_ANY_TAG , MPI_COMM_WORLD , &status);
    // Определение количества элементов
    MPI_Get_count( &status , MPI_DOUBLE , &count);
    ex.clear(); // Очистка буфера
    ex.resize(count); // Резервирование места
    // Получение и запись числового массива в буфер
    MPI_Recv( &ex[0] , count , MPI_DOUBLE , address , 0 , MPI_COMM_WORLD , &status);
    d_it = ex.begin();
    // Проход по буферу и восстановление объектов Particle из их параметров.
    while(d_it != ex.end()){
        particles.emplace_back(d_it, substance);
    }
}

void MPIw::send(const double& d, int address){
    MPI_Send( &d , 1 , MPI_DOUBLE , address , 1 , MPI_COMM_WORLD);
}

void MPIw::recv(double* d, int address){
    MPI_Recv(&d , 1 , MPI_DOUBLE , address , 1 , MPI_COMM_WORLD, &status);
}

void MPIw::send_vec(const std::vector<double>& vec, int address, int tag_){
    MPI_Send( &vec[0], vec.size() , MPI_DOUBLE , address , tag_ , MPI_COMM_WORLD);
}

void MPIw::recv_vec(std::vector<double>& vec, int address, int tag_){
    MPI_Probe( address , tag_ , MPI_COMM_WORLD , &status);
    MPI_Get_count( &status , MPI_DOUBLE , &count);
    ex.clear();
    ex.resize(count);
    MPI_Recv( &ex[0] , count , MPI_DOUBLE , address , tag_ , MPI_COMM_WORLD , &status);
    vec.insert(vec.end(), ex.begin(), ex.end()); //rewrite to move
}

void MPIw::exchangeBlocking(std::vector<Particle>& particles, ParticleSubstance* particle_substance){
    // Выделение частиц не требующих перемещения на другой процесс. 
    iter_0 = std::partition(particles.begin(), particles.end(), [](Particle& p){ return p.MPI_exchange == -1; });
    // Все элементу удовлетворявшие условию расположены до iter_0
    iter_1 = iter_0;
    // Частицы массива до iter_0 в обмене не участвуют
    mpi_income.clear();
    // Цикл по всем возможным комбинациям отправителя-приемника
    // s_rank - отправитель
    // r_rank - приемник
    for(int s_rank = 0; s_rank < size; s_rank++){
        for(int r_rank = 0; r_rank < size; r_rank++){
            if(s_rank == r_rank){
                continue;
            }
            if(s_rank == rank){
                // Рассмытривается область от iter_1 до конца массива
                // Выделение частиц для которых адрес совпадает с текущим значением адреса приемника.
                auto iter_2 = std::partition(iter_1, particles.end(), [&r_rank](Particle& p){return p.MPI_exchange == r_rank;});
                // Все элементу удовлетворявшие условию расположены от iter_1 до iter_2.
                // Отправка частиц
                send(iter_1, iter_2, r_rank);
                // Смещение области для следующей итерации
                iter_1 = iter_2;
            }
            if(r_rank == rank){
                // Прием частиц в буфер
                recv(mpi_income, s_rank, particle_substance);
            }
        }
    }
    // Удаление отправленых частиц
    particles.erase(iter_0,particles.end());
    if(mpi_income.size()){
        // Перенос частиц из буфера в основной вектор
        particles.insert(particles.end(), mpi_income.begin(), mpi_income.end());
    }
    // Очистка буфера
    mpi_income.clear();
}

// Расчет глубины проникновения спрея. Информация собирается со всех процессов на процессе с рангом 0.
double MPIw::pencalc(const std::vector<Particle>& particles, const Inlet& inlet){
    penvec.clear(); // Расстояния от форсунки до частиц
    massvec.clear(); // Массы частиц
    for(const Particle& p: particles){
        // Заполнение векторов
        penvec.push_back((p.cGet()-inlet.cGet()).length());
        massvec.push_back(p.mGet()*p.nGet());
    }
    if(rank != 0){
        // Отправка результатов на ранг 0.
        send_vec(penvec, 0, 3);
        send_vec(massvec, 0, 4);
        penvec.clear();
        massvec.clear();
    }else{
        // Сбор результатов
        for(int i = 1; i < size; i++){
            recv_vec(penvec, i, 3);
            recv_vec(massvec, i, 4);
        }
        // Вектор для сортировки пар по первому элементу пары
        // Т.е. для сортировки по удалению от форсунки
        sortvec.clear();
        sortvec.reserve(penvec.size());
        // Конструирование
        for(int i = 0; i < penvec.size(); i++){
            sortvec.emplace_back(penvec[i], massvec[i]);
        }
        std::sort(sortvec.begin(), sortvec.end());
        all_mass = 0; // Сумарная масса всех частиц в расчете
        mass = 0;
        pen = 0;
        for(const std::pair<double&,double&> p: sortvec){
            all_mass += p.second;
        }
        all_mass *= 0.95;
        for(const std::pair<double&,double&> p: sortvec){
            if(mass >= all_mass){
                break;
            }
            pen = p.first;
            mass += p.second;
        }
        sortvec.clear();
        penvec.clear();
        massvec.clear();
    }
    return pen;
}

// Проверка по всем процессам наличия частиц с неизрасходованым временем шага. 
bool MPIw::particlesToSolve(std::vector<Particle>& particles){
    // Флаг true если хоть одна частица имеет время больше чем 1е-10
    here = std::any_of(particles.begin(), particles.end(), [](Particle& p){return p.timeGet()>1e-10;});
    // Для однопроцессорного запуска
    if(size == 1){
        return here;
    }
    // Сбор значения флага со всех процессов и применение "логического ИЛИ". Запись результата в result.
    MPI_Allreduce(&here, &result, 1, MPI_CXX_BOOL, MPI_LOR, MPI_COMM_WORLD);
    return result;
}