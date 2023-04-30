#pragma once

#include "particle.hpp"
#include "cell.hpp"
#include "surface.h"

#include <random>

/// Интерфейс взаимосдействия свойств и моделей
class Environment{
public:
    virtual double dyn_vis(const double& T) const = 0;
    virtual double lambda(const double& T) const = 0;
    virtual double cp_o2(const double& T) const = 0;
    virtual double cp_n2(const double& T) const = 0;
};

/// Класс свойств воздуха
class Air: public Environment{
private:
    /// Константы полинома теплопроводности
    const double lambda_c[3] = {3.5100e-04, 7.6520e-01, 2.5767e+01};

    /// Константы полинома вязкости
    const double dyn_vis_c[4] = {6.6055e-06, 4.5297e-08, -1.2064e-11, 1.6092e-15};

    /// Константы полинома теплоемкости молекулярного кислорода
    const double o2_gas_spec_heat[5] = {1.2668945e+03, 0.0024062095, -3.6538126e+03, 6.2504293e+03, -2.9740421e+03};

    /// Константы полинома теплоемкости молекулярного азота
    const double n2_gas_spec_heat[5] = {1.3864722e+03, 0.0020343171, -4.093335e+03, 7.4680536e+03, -3.7721195e+03};
    
public:
    /// Молярная масса
    const double M = 0.02897;

    /// Динамическая вязкость
    virtual double dyn_vis(const double& T) const;

    /// Теплопроводность
    virtual double lambda(const double& T) const;

    /// Теплоемоксть молекулярного кислорода
    virtual double cp_o2(const double& T) const;

    /// Теплоемоксть молекулярного азота
    virtual double cp_n2(const double& T) const;
};


/**
 * @brief Абстрактная модель.
 * @details Базовый класс, от которого наследуются все модели.
 * 
 */
class Model{
protected:
    const Air* env; // Адрес объекта свойств воздуха
    const double fq; // Величина обратная шагу интегрирования модели
    Model(const Air& env_, const double& fq_): env(&env_), fq(fq_) {};

public:
    virtual void solve(Particle& p, Cell& c) = 0; ///< Решение модели на данном шаге.
};

/**
 * @brief Модель движения.
 * @details Отвечает за движение частицы, пересечение с границами, взаимодействие со стенками.
 * 
 */

class Momentum: public Model{
private:
    Vec a;     ///< Ускорение частицы
    Vec dv; ///< Изменение скорости 
    Vec dx; ///< Изменение координаты
    double Re; ///< Число Рейнольдса
    double Cd; ///< Коэффициент сопротивления
    double count;

    Table sin_alpha_k_table; // Таблица для расчета взаимодействия со стенкой.

    /// STL CROSS
    std::vector<double> x0, x1, xn;
    Vec x1v, xnv;
	std::mt19937 gen; //Standard mersenne_twister_engine seeded with rd()
	std::uniform_real_distribution<double> dis0;  // distribution from 0 to 1

public:
    Momentum(const Air& env_, const double& fq_, std::string filepath): Model(env_, fq_), sin_alpha_k_table(filepath) {
        x0.reserve(3);
        x1.reserve(3);
        xn.reserve(3);
        gen.seed(0);
        //gen.seed(time(0));
        dis0 = std::uniform_real_distribution<double>(0.0, 1.0);
    };
    virtual void solve(Particle& p, Cell& c) {};
    /// Ускорение сопротивления
    Vec drag(Particle& p, Cell& c); 
    /// Пересечение с поверхностью 
    std::tuple<Surface*, Vec> checkCross(Particle& p, std::vector<Surface>& ss_vec, int& iskip);
    /// Взаимодействие со стенкой
    void wallInteraction(Particle& p, Surface& ss, int& ielem);
};

/**
 * @brief Модель испарения (Дукович).
 * @details Управляет нагревом и испранием частиц.
 * 
 */
class Dukowicz: public Model{
friend class LogSingle;
private:
    double q;  ///< Тепловой поток
    double Re; ///< Число Рейнольдса
    double Pr; ///< Число Прандтля
    double Nu; ///< Число Нуссельта
    double f_q; ///< Отношение потока массы к потоку тепла
    double By; 
    double Yvs; ///< Массовая доля пара на поверхности капли
    double f_vs; ///< Отношение давления насыщенного пара капли при данной температуре к давлению в ячейке
    double dT, dm;
    double T;
    double o, a;
    double fmol, omol, aimol, pmol;
    double cp_film;
    double lambda_film;
    double dyn_vis_film;

    const double C1;

public:
    Dukowicz(const Air& env_, const double& fq_, double C1_ = 1): Model(env_, fq_), C1(C1_) {};
    double get_cp_film(const double& T,const double& f, const double& o, const double& a, const double& ppp, const Particle& p);
    double get_lambda_film(const double& T,const double& f, const double& o, const double& a, const double& ppp, const Particle& p);
    double get_dyn_vis_film(const double& T,const double& f, const double& o, const double& a, const double& ppp, const Particle& p);
    virtual void solve(Particle& p, Cell& c);
};

/**
 * @brief Модель дробления (Wave)
 * @details Управляет дроблениме частиц. Модель расчитывает новый диаметр частицы и изменяет количество частиц в парцеле исходя из сохранения массы.
 * Количество логических частиц в расчете отстается неизменным.
 */
class Wave: public Model{
private:
    double Re; // Число Рейнольдса
    double We, We1, We2; ///< Число Вебера
    double Lambda;
    double Omega;
    double Oh; // Число Онезорге
    double T_;
    double tau_a;
    double dd;
    // Эмперические константы
    const double C1;
    const double C2;

public:
    Wave(const Air& env_, const double& fq_, double C1_ = 0.61, double C2_ = 15): Model(env_, fq_), C1(C1_), C2(C2_) {};
    virtual void solve(Particle& p, Cell& c);
};

/**
 * @brief Модель турбулентной дисперсии.
 * @details Добавляет к скорости частицы случайный вектор. Модель расчитывает вектор тубулентной скорости и время жизни вихря.
 * На длительность времени жизни вихря скорость остается постоянной, по истечении - пересчитывается вместе со временем.
 */
class TurbulentDispersion: public Model{
private:
    // Генратор случайных чисел
    std::mt19937_64 gen;
    // Равномерное распределение
    std::normal_distribution<double> normal_dist{0.,1.};
    double Ct = 0.164;
    double u_;
    double u, k, eps;
public:
    TurbulentDispersion(const Air& env_, const double& fq_): Model(env_, fq_) {
        gen.seed(0);
        //gen.seed(time(0));
    }
    virtual void solve(Particle& p, Cell& c);
};