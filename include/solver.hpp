#pragma once

#include "init.hpp"
#include "mpi-wrapper.hpp"

#include "particle.hpp"
#include "cell.hpp"
#include "models.hpp"
#include "metrics.hpp"
#include "timing.hpp"
#include "kdtree.hpp"

#include <ctime>

class Solver{
private:
	Init& init;
	MPIw& mpi;

	ParticleSubstance particle_substance;

	Air air; // Воздух

	std::vector<Surface> surface; // STL поверхности
	Surface* ss_crossed; // Указатель на пересеченную поверхность
	double dt_to_cross; // Время до персечения
	Vec cross_point; // Точка пересечения
	
	std::vector<Cell> cells; // Ячейки расчетной сетки
	HyperCell hpc; // Интерполяционная ячейка
	double lg;
	std::unordered_map<Cell*, SourceStorage> source_storage; // Временное хранилище источниковых членов
	CellsMap cells_map;
	std::vector<Cell*> cells_near;

	TurbulentDispersion turbulent_dispersion;
	Momentum momentum;
	Dukowicz dukowicz;
	Wave wave;

	// Вектор частиц
	std::vector<Particle> particles;

	Inlet inlet;

	int total_iterations; // Общее количество итераций
	double current_time = 0;
	double Cu; ///< Число Куранта
	double Cu_treshold = 1; ///< Множитель числа Куранта

	Vec drag; // Сила сопротивления

	bool solving_flag; // Флаг наличия в расчете по всем процессам хотя бы одной недосчитаной частицы.

	LogSingle logger;
	LogSingle penetration;

	Timing calc_timer;
	Timing cells_timer;

	std::vector<Cell*> search_results;
	
	std::clock_t start;
	std::clock_t end;
	double duration = 0;

public:
	Solver( Init& init_, MPIw& mpi_): init(init_), mpi(mpi_), 
		turbulent_dispersion(TurbulentDispersion(air, init_.fq)), momentum(Momentum(air, init_.fq, "../resources/sina_k_table")),
		dukowicz(Dukowicz(air, init_.fq, init_.duk_C1)), wave(Wave(air, init_.fq, init_.wave_C1, init_.wave_C2)),
		particle_substance(ParticleSubstance("../resources/substances/" + init.substance_name + "/"))
		{
		init.readSurface(surface, mpi.getRank(), mpi.getSize());
		init.readCells(cells, mpi.getRank());
		source_storage.reserve(cells.size());
		cells_map.makeMap(cells);
		particles.reserve(1000);
		total_iterations = std::ceil(init.total_time  / init.dt_global); // Общее количество итераций

		penetration.set_values_names("time,pen");	

		search_results.reserve(50);	
	};

	void run(){
		// ОСНОВНОЙ ЦИКЛ

		KDTree kdtree(cells);
		/*
		std::vector<Cell*> test = kdtree.search(Vec(0,0,0), init.connection_range);
		for(auto c: test){
			if((Vec(0,0,0) - c->cGet()).length() <= init.connection_range){
				std::cout << c->cGet().print() << "\n";

			}
		}
		return;
		*/
		

		Inlet inlet(init.inlet_coord, init.inlet_alpha, init.inlet_d, init.inlet_G,
			init.inlet_T, &particle_substance, init.fq, init.inletFlag(cells), init.inlet_n);
		calc_timer.start();
		for(int i = 0; i < total_iterations; i++){
			if(!init.single_particle_flag){
				inlet.generateParcel(particles, init.dt_global); // Генерация частиц в расчет
			}
			// Добавление каждой частице глобального шага интегрирования.
			// Он задается как параметр расчета и расходуется каждой частицой индивидульно в зависимости от ее локаьного шага.
			std::for_each(particles.begin(), particles.end(), [this](Particle& p){p.timeSet() += this->init.dt_global;});

			solving_flag = true; // Установка значения для входа в первую итерацию
			while(solving_flag){
				for(Particle& p: particles){
					
					while(p.timeGet()>1e-10){ // Проверка нерасходованого глобального шага
						
						cells_timer.start();
						start = std::clock();

						
						// Поиск по KD-Tree
						search_results.clear();
						search_results = kdtree.search(p.cGet(), init.connection_range);
						//std::cout << p.cGet().print() << " | size: " << search_results.size() << "\n";
						for(Cell* cell: search_results){
							lg = (p.cGet() - cell->cGet()).length(); // Расстояние от частицы до центра ячейки
							if(lg <= init.connection_range){
								// Добавление адреса ячейки в список, вместе с весовым коэффициентом = 1/расстояние.
								hpc.add(cell, 1/lg);
							}
						}
						/*
						// Перебор всех ячеек для поиска входящих в окрестность connection_range
						for(Cell& cell: cells){
							lg = (p.cGet() - cell.cGet()).length(); // Расстояние от частицы до центра ячейки
							if(lg <= init.connection_range){
								// Добавление адреса ячейки в список, вместе с весовым коэффициентом = 1/расстояние.
								hpc.add(&cell, 1/lg); 
							}
						}
						*/
						/*
						cells_map.search(p.cGet(), init.connection_range);
						for(Cell* cell: cells_map.result){
							lg = (p.cGet() - cell->cGet()).length(); // Расстояние от частицы до центра ячейки
							if(lg <= init.connection_range){
								// Добавление адреса ячейки в список, вместе с весовым коэффициентом = 1/расстояние.
								hpc.add(cell, 1/lg); 
							}
						}
						*/
						cells_timer.pause();
						end = std::clock();
						duration += (double)(end - start) / CLOCKS_PER_SEC;
						try{
							// Конструирование ячейки путем весовой интерполяции ячеек из окрестности.
							// Выбрасывает исключение если список пуст.
							hpc.construct(); 
						}
						// Обработка исключения пустого списка.
						catch(int){
							std::cerr << "n=" << p.count << " on proc " << mpi.getRank() <<  "\n";
							p.setOutOfBounds(); // Флаг частицы за границами расчета
							break;
						}
						
						if(init.turbulent_dispersion_flag){
							turbulent_dispersion.solve(p, hpc); // Не зависит от dt, поэтому может быть выполнено до коррекции шага.
						}
						// Установка локального шага.
						if( p.timeGet() >= init.dt_global ){
							p.dtSet() = init.dt_global;
						} else{
							p.dtSet() = p.timeGet();
						}
						drag = momentum.drag(p, hpc); // Ускорение от сопротивления

						// Провека пересечения траетории частицы с поверхностями. Возвращает nullprt если пересечений нет.
						std::tie(ss_crossed, cross_point) = momentum.checkCross(p, surface, p.iskip); 

						if(ss_crossed != nullptr){ 
							dt_to_cross = ( cross_point - p.cGet() ).length() / p.vGet().length(); // Времемя до столкнования
							switch(ss_crossed->get_kind()){ // Выбор алгоритма в зависимости от типа поверзности
								case SURFACEKIND_WALL: // Стенка
									p.vSet() += drag * dt_to_cross; // Расчет скорости для перемещения до стенки
									p.cSet() += p.vGet() * dt_to_cross; // Расчет координаты для перемещения до стенки
									hpc.recvMomentum(drag * p.nGet()); // Запись источника в уравнении сохранения импульса (нормировано на секунду)
									momentum.wallInteraction(p, *ss_crossed, p.iskip); // Отражение от стенки (изменяет скорсть частицы)
									p.cSet() += p.vGet() * (p.dtGet() - dt_to_cross); // Пермещение от стенки
									break;
								case SURFACEKIND_MPI: // Граница MPI		
									// Движение до границы
									p.dtSet() = dt_to_cross;
									p.vSet() += drag * p.dtGet();
									p.cSet() = cross_point;
									// Запись источника в уравнении сохранения импульса (нормировано на секунду)
									hpc.recvMomentum(drag * p.nGet());
									// Установка адреса процесса на который пермещается частица в качестве флага перемещения
									p.MPI_exchange = ss_crossed->getAddress();
									break;								
								case SURFACEKIND_INLET: // Вход
									// Движение до границы
									p.dtSet() = dt_to_cross;
									p.vSet() += drag * dt_to_cross;
									p.cSet() += p.vGet() * dt_to_cross;
									hpc.recvMomentum(drag * p.nGet()); // norm by 1 second
									p.setOutOfBounds(); // Установка флага выхода за границу расчета
									break;
								case SURFACEKIND_OUTLET: // Выход
									// Движение до границы
									p.dtSet() = dt_to_cross;
									p.vSet() += drag * dt_to_cross;
									p.cSet() += p.vGet() * dt_to_cross;
									hpc.recvMomentum(drag * p.nGet()); // norm by 1 second
									p.setOutOfBounds();// Установка флага выхода за границу расчета
									break;
							}
						}else{
							p.vSet() += drag * p.dtGet(); // Расчет скорости
							p.cSet() += p.vGet() * p.dtGet(); // РАсчет коордианты
							hpc.recvMomentum(drag * p.nGet()); // norm by 1 second
							p.iskip = -1;
						}
						
						if(init.evap_flag){
							dukowicz.solve(p, hpc);
						}
						
						if(init.brup_flag){
							wave.solve(p, hpc);
						}

						hpc.releaseToStorage(source_storage); // Расчет по весовым коэф. и запись в хранилище источниковых членов

						p.makeTimeStep(); // Совершение локального шага по глобальному. (Глобальный шаг частицы -= локальный шаг частицы)
						
						if(p.OOB()){ break; } // Выход если установлен флаг outOfBounds
						if(p.MPI_exchange != -1){ break; } // Выход если записан адрес пересылки частицы
					}
				}
				// Удаление всех outOfBounds
				particles.erase(std::remove_if(particles.begin(), particles.end(), OutsidePredicate), particles.end()); 
				// Межпроцессорный обмен частицами
				mpi.exchangeBlocking(particles, &particle_substance); 
				// Проверка по всем процессам наличия частиц с неизрасходованным глобальным шагом.
				solving_flag = mpi.particlesToSolve(particles); 
				break;
			}
			// Запись источниковых членов из хранилища в ячейки
			for(auto& s: source_storage){
				s.second.release(s.first);
			}
			source_storage.clear();
		}
		calc_timer.pause();
		std::cout << "time spend on cells search: " << duration << "\n";
		std::cout << "DONE on "<< mpi.getRank() <<" with particle count " << particles.size() << "\n";
	};

	void results(){
		
		std::string string_rank = std::to_string(mpi.getRank());
		// Запись координат частиц в конце расчета
		logger.write_particle(particles, "../results-test/" + init.taskname + "-coords" + string_rank + ".csv");
		
		if(mpi.getRank() == 0){
			// Запись глубины проникновения. Результаты собираются на ранге = 0.
			penetration.write_value("../results-test/" + init.taskname + "-pen.csv");
		}
	};
};