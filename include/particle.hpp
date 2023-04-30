#pragma once

#include "vec.hpp"
#include "property.hpp"

/** Частица. 
 * Класс моделирующий физические частицы через парцелы. 
 * Ведут себя как одна частица, но представляют некоторое количество реальных одинковых физических частиц.
 */
class Particle {
protected:
    Vec c; ///< Координата.
    Vec v; ///< Скорость.
    double d = 0; ///< Диаметр частицы.
    double m = 0; ///< Масса
    double T = 0; ///< Температура капли
    double n = 0; ///< Число физических частиц в логической

    bool outOfBounds = false; ///< Флаг выхода частицы за границу рассчетной зоны.

    // Для модели турбулентной дисперсии.
    Vec v_pulse = {0,0,0}; ///< Вектор скорости турбелентной пульсации
    double time_pulse = 0; ///< Время нахождения в вихре

    double time; ///< Нерасчитанное время для данной частицы (глобальный шаг)
    double dt; ///< Шаг интегрирования по времени для частицы (локальный шаг интегрирования по времени)


public:
    /**
     * @brief Флаг межпроцессорного обмена.
     * По умолчанию равен -1, то есть парцел в переносе 
     * на другой процесс не нуждается.
     * При необходимости перехода на другой процесс
     * в данную переменную записывается ранг процесса [0, size).
     * 
     */
    int MPI_exchange = -1;

    int count;
    int iskip = -1;

    static double all_evap;
    static double all_inj;

    ParticleSubstance* props; 

    const Vec& vGet() const{
        return v;
    }
    Vec& vSet() {
        return v;
    }
    const double& dGet() const{
        return d;
    }
    double& dSet(){
        return d;
    }
    const double& mGet() const{
        return m;
    }
    const double& TGet() const{
        return T;
    }

    const double& nGet() const{
        return n;
    }
    Vec& cSet(){
        return c;
    }
    const Vec& cGet() const {
        return c;
    }
    Vec& v_pulseSet(){
        return v_pulse;
    }
    const Vec& v_pulseGet() const {
        return v_pulse;
    }
    const double& time_pulseGet() const{
        return time_pulse;
    }
    double& time_pulseSet(){
        return time_pulse;
    }
    const double& timeGet() const{
        return time;
    }
    double& timeSet(){
        return time;
    }
    const double& dtGet() const{
        return dt;
    }
    double& dtSet(){
        return dt;
    }

    /** Временной шаг. 
     * Совершение локального шага по глобальному. (Глобальный шаг частицы -= локальный шаг частицы)
     */
    void makeTimeStep(){
        time -= dt;
        if(time_pulse != 0){
            time_pulse -= dt;
        }
    }

    /// Устанавливает флаг выхода за границу рассчетной зоны
    void setOutOfBounds(){
        outOfBounds = true;
    }

    /// Возвращает флаг outOfBounds
    const bool& OOB(){
        return outOfBounds;
    }

    /// Предикат выхода за границу рассчетной зоны
    friend bool OutsidePredicate(const Particle& p) {
        return p.outOfBounds;
    }

    /// Предикат наличия частицы в рассчетной зоне
    friend bool InsidePredicate(const Particle& p) {
        return !p.outOfBounds;
    }

    /// Конструктор
    Particle(ParticleSubstance* props_, const double& d_, const Vec& c_, const Vec& v_, const double& temp_, const double& n_, double time_ = 0, int count_ = 0): 
        props(props_), d(d_), c(c_), v(v_), T(temp_), n(n_), time(time_), count(count_) { 
        m = props->ro[T] * pi * d * d * d / 6;
        all_inj += m * n;
        iskip = -1;
    }

    /** Десериализация. 
     * Используется при восстановлении объектов Particle из массива чисел при межпроцессорном обмене.
     * Должен быть обратен методу Particle::serialize().
     * @param it Итератор указывающий на первый элемент массива чисел для конструирования частиц
     * @param props_ указатель на свойства вещества частицы
     */
    Particle(std::vector<double>::iterator& it, ParticleSubstance* props_):
        c(Vec(*it, *(it+1), *(it+2))), 
        v(Vec(*(it+3), *(it+4), *(it+5))), 
        d(*(it+6)), T(*(it+7)), n(*(it+8)), props(props_),
        v_pulse(Vec(*(it+9), *(it+10), *(it+11))), time_pulse(*(it+12)), time(*(it+13)), iskip((int)*(it+14)), count((int)*(it+15)) { 
            m = props->ro[T] * pi * d * d * d / 6;
            it = it + 16;
            outOfBounds = false;
            MPI_exchange = -1;
    }
    

    /** Сериализация.
     * Используется при межпроцессорном обмене для представления объектов 
     * в виде строгой ограниченой последовательности чисел.
     * Метод должен быть обратен конструктору из вектора чисел.
     * @param[out] ex Вектор в который дописываются сериализованные частицы
     */
    void serialize(std::vector<double>& ex) const {
        c.serialize(ex);
        v.serialize(ex);
        ex.insert(ex.end(), {d, T, n});
        v_pulse.serialize(ex);
        ex.insert(ex.end(), {time_pulse, time, (double)iskip, (double)count});
    }

    void evaporation(const double& dT, const double& dm);

    void breakUp(const double& dd);

    std::string properties() const;

    std::string propertiesNames() const;

    bool operator<(const Particle& p) const;
};

// Внешнее объявление.

bool InsidePredicate(const Particle&); 
bool OutsidePredicate(const Particle&);