#pragma once

#include <unordered_map>
#include <algorithm>
#include <iostream>

#include "vec.hpp"

/** Рассчетная ячейка.
 */
class Cell{
friend class CellsMap;
protected:
    Vec c; ///< Координата центра ячейки
    Vec v; ///< Скорость в ячейке
    double V; ///< Объем ячейки
    double p; ///< Давление в ячейке
    double ro; ///< Плотность
    double T; ///< Температура
    double Yv = 0; ///< Концентрация паров топлива
    double k = 0;
    double eps = 0;
    
    Vec momentum_source;
    double heat_source = 0;
    double mass_source = 0;

    bool exchange_marker = false; // Маркер состоявшегося обмена источниками

public:
    const Vec& cGet() const{
        return c;
    }
    const Vec& vGet() const{
        return v;
    }
    const double& VGet() const{
        return V;
    }
    const double& pGet() const{
        return p;
    }
    const double& roGet() const{
        return ro;
    }
    const double& TGet() const{
        return T;
    }
    const double& YvGet() const{
        return Yv;
    }
    void recvMass(const double& d){
        mass_source += d;
        exchange_marker = true;
    }
    void recvHeat(const double& d){
        heat_source += d;
        exchange_marker = true;
    }
    void recvMomentum(const Vec& d){
        momentum_source += d;
        exchange_marker = true;
    }
    const double& massGet() const {
        return mass_source;
    }
    const double& heatGet() const {
        return heat_source;
    }
    const Vec& momentumGet() const {
        return momentum_source;
    }
    void cleanSources(){
        momentum_source = Vec();
        heat_source = 0;
        mass_source = 0;
        exchange_marker = false;
    }

    const bool& markerGet() const {
        return exchange_marker;
    }
    
    /** Конструктор по умолчанию. 
     * Создает ячейку в нулевой координате, со всеми нулевыми параметрам.
     */
    Cell(): c(Vec(0,0,0)), v(Vec(0,0,0)), V(0), p(0), ro(0), T(0), Yv(0), k(0), eps(0) {}

    /// Конструктор.
    Cell(Vec c_, Vec v_, double vol_, double pres_, double ro_, double temp_, 
        double k_ = 0, double eps_ = 0, Vec momentum_source_ = {0,0,0}, double heat_source_ = 0, double mass_source_ = 0): 
        c(c_), v(v_), V(vol_), p(pres_), ro(ro_), T(temp_), k(k_), eps(eps_), 
        momentum_source(momentum_source_), heat_source(heat_source_), mass_source(mass_source_) {}

    void clear(){
        c = Vec(0,0,0);
        v = Vec(0,0,0); 
        V = 0;
        p = 0;
        ro = 0;
        T = 0;
        k = 0;
        eps = 0;
        momentum_source = Vec(0,0,0);
        heat_source = 0;
        mass_source = 0;
    }

    /** Сумма. 
     * Возвращает ячейку с параметрами равными суммам параметров слагаемых.
     */
    Cell operator+(const Cell& rv) const{
        return Cell(c + rv.c, 
            v + rv.v,
            V + rv.V,
            p + rv.p,
            ro + rv.ro,
            T + rv.T,
            k + rv.k,
            eps + rv.eps,
            momentum_source + rv.momentum_source,
            heat_source + rv.heat_source,
            mass_source + rv.mass_source
        );
    }

    /** Сумма. 
     * Присваивает данной ячейке параметры равные сумме параметров.
     */
    void operator+=(const Cell& rv) {
        c += rv.c;
        v += rv.v;
        V += rv.V;
        p += rv.p;
        ro += rv.ro;
        T += rv.T;
        k += rv.k;
        eps += rv.eps;
        momentum_source += rv.momentum_source;
        heat_source += rv.heat_source;
        mass_source += rv.mass_source;
    }

    /** Умножение. 
     * Возвращает ячейку с параметрами равными произведению параметров и числа.
     */
    template <typename Type>
    Cell operator*(const Type& d) const{
        return Cell(c * d, 
            v * d,
            V * d,
            p * d,
            ro * d,
            T * d,
            k * d,
            eps * d,
            momentum_source * d,
            heat_source * d,
            mass_source * d
        );
    }

    /** Деление. 
     * Возвращает ячейку с параметрами равными частному параметра и числа.
     */
    template <typename Type>
    Cell operator/(const Type& d) const{
        return Cell(c / d, 
            v / d,
            V / d,
            p / d,
            ro / d,
            T / d,
            k / d,
            eps / d,
            momentum_source / d,
            heat_source / d,
            mass_source / d
        );
    }

    /** Деление. 
     * Присваивает данной ячейке параметры равные частному параметра и числа.
     */
    template <typename Type>
    void operator/=(const Type& d) {
        c /= d;
        v /= d;
        V /= d;
        p /= d;
        ro /= d;
        T /= d;
        k /= d;
        eps /= d;
        momentum_source /= d;
        heat_source /= d;
        mass_source /= d;
    }

    static bool coordCompare(const Cell& a, const Cell& b){
        return a.c < b.c;
    }
    static bool coordComparePointers(const Cell* a, const Cell* b){
        return a->c < b->c;
    }

    static bool compareByX(const Cell* a, const Cell* b){
        return a->c[0] < b->c[0];
    }
    static bool compareByY(const Cell* a, const Cell* b){
        return a->c[1] < b->c[1];
    }
    static bool compareByZ(const Cell* a, const Cell* b){
        return a->c[2] < b->c[2];
    }

    static bool coordCompareVec(const Cell& a, const Vec& b){
        return a.c < b;
    }
};

// Чтение ячеек из файла.
void readCellsAll(std::string filepath, std::vector<Cell>& cells);

/// Хранилище источников. Позволяет накапливать и отдавать источники в ячейку.
struct SourceStorage{
    Vec momentum;
    double heat;
    double mass;
    
    /// Обнуление.
    void clear(){
        momentum.clear();
        heat = 0;
        mass = 0;
    }

    /// Возврат в ячейку
    void release(Cell* c){
        c->recvMomentum(momentum);
        c->recvHeat(heat);
        c->recvMass(mass);
        clear();
    }
};

class CellsMap{
private:
    std::vector<Cell*> mapX;
    std::vector<Cell*> mapY;
    std::vector<Cell*> mapZ;
    int left, right, mid;

    int axis_binary_search(const double coord, const int axis){
        left = 0;
        right = mapX.size() - 1;
        while (left < right) {
            mid = (left + right) / 2;
            if((mapX[mid-1]->c[axis] < mapX[mid]->c[axis]) && (mapX[mid]->c[axis] >= coord)){
                return mid;
            }
            if(mapX[mid]->c[axis] > coord){
                right = mid;
            }
            else{
                left = mid;
            }
        }
        return -1;
    }
public:
    std::vector<Cell*> result;

    void makeMap(std::vector<Cell>& cells){
        mapX.reserve(cells.size());
        mapY.reserve(cells.size());
        mapZ.reserve(cells.size());
        for(Cell& a: cells){
            mapX.push_back(&a);
        }
        mapY = mapX;
        mapZ = mapX;
        std::sort(mapX.begin(), mapX.end(), Cell::compareByX);
        std::sort(mapY.begin(), mapY.end(), Cell::compareByY);
        std::sort(mapZ.begin(), mapZ.end(), Cell::compareByZ);
    }

    void search(const Vec& v, const double& range){
        int x_left = axis_binary_search(v[0] - range, 0);
        result.clear();
        if(x_left != -1){
            int x_right = x_left;
            while(mapX[x_right]->c[0] <= (v[0]+range)){
                x_right++;
            }
            for(x_left; x_left <= x_right; x_left++){
                if((mapX[x_left]->c[1] <= v[1] + range) && (mapX[x_left]->c[1] >= v[1] - range)
                    && (mapX[x_left]->c[2] <= v[2] + range) && (mapX[x_left]->c[2] >= v[2] - range)){
                        result.push_back(mapX[x_left]);
                }
            }
        }
    }
};


/**
 * @brief Сверхячейка.
 * @details Абстрактная ячейка, полученная интерполяцией ячеек с области в окрестности частицы.
 */
class HyperCell: public Cell{
private:
    std::vector<std::pair<Cell*, double>> list; // Список: адрес ячейки, вес.
    double wSum = 0;

public:
    HyperCell();

    void clear(); // Очистка

    std::size_t size(); // Размер

    void add(Cell* c, double w); // Добавление ячейки и ее веса в список интерполяции

    void construct(); // Интерполяция 

    void releaseToStorage(std::unordered_map<Cell*, SourceStorage>& storage); // Перенос значений в хранилище

    void releaseToCells(); // Перенос значений в ячейку
};