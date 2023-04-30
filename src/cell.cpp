#include "../include/cell.hpp"

#include <fstream>
#include <sstream>

/** Чтение ячеек. 
 * Считывает данные из файла и конструирует массив ячеек.
 */
void readCellsAll(std::string filepath, std::vector<Cell>& cells){
    std::ifstream file(filepath);
    std::vector<double> s;
    double val;
    std::stringstream ss;
    std::string buffer;
    std::getline(file, buffer);
    while(std::getline(file, buffer)){
        ss.clear();
        ss << buffer;
        s.clear();
        while(ss >> val){
            s.push_back(val);
        }
        cells.emplace_back(
            Vec(s[0], s[1], s[2]),
            Vec(s[8], s[9], s[10]),
            s[3], s[11], s[12], s[13]);
    }
    file.close();
};


HyperCell::HyperCell(){
    list.reserve(10);
};

void HyperCell::clear(){
    Cell::clear();
    list.clear();
    wSum = 0;
};

std::size_t HyperCell::size(){
    return list.size();
};

void HyperCell::add(Cell* c, double w){
    // с - адрес ячейки, w - вес
    list.emplace_back(c,w);
};


void HyperCell::construct(){
    // Если список пуст, выбрасывается исключение. Обрабатывается в main.cpp.
    // ИСключение прерывает дальнейшее выполнение метода.
    if(list.size() == 0){
        throw 0;
    }
    // Так как класс наследуется от Cell, то он содержит в себе все поля и методы Cell.
    // Другими словами, внутри HyperCell есть объект Cell, к которому здесь обращаемся.
    for(std::pair<Cell*, double> p: list){
        // Умножение ячейки на ее вес и добавление к суммарной ячейке.
        Cell::operator+=((*p.first) * p.second);
        wSum += p.second; // Сумма весов
    }
    // Деление суммы на сумму весов.
    Cell::operator/=(wSum);
};

// Распределение долей источниковых сленов в соответствии с весами.
// Запись в промежуточное хранилище.
void HyperCell::releaseToStorage(std::unordered_map<Cell*, SourceStorage>& storage){
    momentum_source = momentum_source *  wSum;
    heat_source *= wSum;
    mass_source *= wSum;
    for(std::pair<Cell*, double> p: list){
        storage[p.first].momentum += momentum_source * p.second;
        storage[p.first].heat += heat_source * p.second;
        storage[p.first].mass += mass_source * p.second;
    }
    clear();
};

void HyperCell::releaseToCells(){
    momentum_source = momentum_source *  wSum;
    heat_source *= wSum;
    mass_source *= wSum;

    for(std::pair<Cell*, double> p: list){
        p.first->recvHeat(heat_source * p.second);
        p.first->recvMass(mass_source * p.second);
        p.first->recvMomentum(momentum_source * p.second);
    }

    clear();
};