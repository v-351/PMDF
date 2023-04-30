#pragma once

#include <vector>
#include <map>
#include <string>
#include <fstream>
#include <cmath>


/**
 * @brief Физические свойства. 
 * Читает из csv-файлов. Предоставляет доступ по индексации через температуру.
 */
class Property{
private:
    std::vector<double> pr; // Массив значений
    int Tmin, Tmax; // Минимальное и максимальное значения индекса
    int Tstep = 1; // Шаг индекса
public:
    Property(){};

    // Чтение данных из csv файла
    void readCSV(std::string filepath, char sep = ',', int Tstep_ = 1);
    
    Property(std::string filepath, char sep = ',', int Tstep_ = 1) { readCSV(filepath, sep, Tstep_); };

    /// Доступ через температуру.
    const double& operator[](const double& T) const;

    void print() const;
};

/**
 * @brief Таблица.
 * Читает csv и предоставляет доступ.
 * 
 */
class Table{
private:
    double index, value;
    std::map<double, double> valmap;
public:
    Table(std::string filepath){
        parse(filepath);
    }
    /// Чтение файла
    void parse(std::string filepath);
    /// Возврат значения
    double k(const double&) const;
};

/**
 * @brief Вещество частицы.
 * Представляет вещество через набор его свойств.
 * 
 */
class ParticleSubstance{
public:
    double M; ///< Молярная масса
    Property p_sat; ///< Давление насыщеных паров
    Property cp; ///< Теплоемкость жидкости
    Property cp_vap; ///< Теплоемкость пара
    Property ro; ///< Плотность вещества частицы.
    Property L; ///< Теплота парообразования
    Property sigma; ///< Коэф. поверхностного натяжения
    Property lambda; ///< Коэф. теплопроводности жидкости
    Property lambda_vap; ///< Коэф. теплопроводности пара
    Property dyn_vis; ///< Динамическая пара жидкости
    Property dyn_vis_vap; ///< Динамическая вязкость пара
    ParticleSubstance(std::string folder);
};

/**
 * @brief Вещество ячейки.
 * Представляет вещество через набор его свойств.
 * 
 */
class CellSubstance{
public:
    const double Yv; // Паросодержание
    const double dyn_vis; // Динамическая вязкость
    const double lambda; // Теплопроводность
    const double M; ///< Молярная масса
    Property cp; // Теплоемкость

    CellSubstance(const double& Yv_, const double& dyn_vis_, const double& lambda_, const std::string cp_filepath, const double& M_):
        Yv(Yv_), dyn_vis(dyn_vis_), lambda(lambda_), cp(cp_filepath), M(M_) {}
};