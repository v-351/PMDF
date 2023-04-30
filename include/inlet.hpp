#pragma once

#include <random>

#include "particle.hpp"
#include "vec.hpp"

/** Форсунка.
 * Генерирует и вводит частицы в поток в соответствии с конфигурацией расчета. 
 * Может эмулировать работу реальной форсунки используя генератор случайных чисел.
 * Генерирует частицы эквивалентные физическим, либо логические, представляющие некоторое, 
 * возможно нецелое, число частиц расчитанное на основе заданного расхода.
 * 
 */
class Inlet {
private:
    double d; ///< Диаметр форсунки
    double d_droplet; ///< Диаметр капель
    double G; ///< Расход
    Vec c; ///< Точка входа частиц в поток
    Vec c_return; ///< Рандомизированая координата генерируемого парцела
    Vec normal, nB, nC;
    double T;
    double v; ///< Модуль скорости
    ParticleSubstance* props; // Ссылка на свойства вещества частиц
    double n_parcels = 1; ///< Количество парцелов генерируемых за шаг
    int i = 0; ///< Счетчик

    std::mt19937_64 gen; ///< Генератор случайных чисел
    std::uniform_real_distribution<double> alpha_rand;
    std::uniform_real_distribution<double> eta_rand;
    std::uniform_real_distribution<double> r_rand;

    bool flag; ///< Флаг вхождения координаты форсунки в данную рассчетную область
    double n; ///< Число физических частиц приходящихся на одну логческую за один шаг
    double n_sec;
    double alpha_gen, eta_gen, v_side;
    Vec v_return;

    int count = 0;

public:

    Inlet(){};

    /** Конструктор. 
     * Создает объект, проверяет вхождение координаты форсунки в расчетную область, инициализирует генератор случайных чисел и распределения, 
     * вычисляет количество физических частиц в одной логической.
     */
    Inlet(const Vec& c_, const double& alpha, const double& d_, const double& G_, const double& T_,
        ParticleSubstance* props_, const double& fq, bool flag_, Vec normal_ = {1,0,0}, 
        double v_ = 0, double d_droplet_ = 0, int n_parcels_ = 1): 
        c(c_), normal(normal_), d(d_), T(T_), props(props_), n_parcels(n_parcels_) {
            flag = flag_; // Флаг включения форсунки. Определяется наличием форсунки на данном процессе.
            if(flag){
                d_droplet = d;
                if(d_droplet_ != 0){
                    d_droplet = d_droplet_;
                }
                n_sec = 6. * G_ / (props->ro[T_] * pi * std::pow(d_droplet, 3)); // Число частиц в секунду
                n = n_sec / fq; // Число частиц за шаг
                // Если скорость не задана вручную, то производится расчет.
                if(v_ == 0){
                    v = G_  / props->ro[T_] / (pi * d_ * d_ / 4); // Начальная скорость частицы.
                }else{
                    v = v_;
                }

                // Инициализация параметров распределений
                alpha_rand = std::uniform_real_distribution<double>(-alpha*pi/180,alpha*pi/180);
                eta_rand = std::uniform_real_distribution<double>(0,2*pi);
                r_rand = std::uniform_real_distribution<double>(0,d/2);
                //gen.seed(time(0));
                gen.seed(0);
                // Подготовка векторов
                normal = normal.normalize(); // Нормализация
                nB = normal.norm(); // Получение перпендикулярного вектора
                nC = dot(normal, nB); // Получение перпендикулярного вектора
            }
    }

    /// Генерация логических частиц в поток.
    void generateParcel(std::vector<Particle>& ex, const double& dt){
        if(flag){
            for(i = 0; i<n_parcels; i++){
                // Рандомизация углов
                alpha_gen = alpha_rand(gen); 
                eta_gen = eta_rand(gen);
                // Получение рандомного единичного вектора в границах заданого телесного угла
                v_return = normal * std::cos(alpha_gen) + (nB * std::cos(eta_gen) + nC * std::sin(eta_gen)) * std::sin(alpha_gen);
                // Умножение вектора на модуль скорости
                v_return = v_return * v;
                // Добавление парцела в массив
                ex.emplace_back(props, d_droplet, c, v_return, T, n/n_parcels, -i*dt/n_parcels, ++count);
            }
        }
    }

    /// Предикат для сортировки частиц по удалению от форсунки.
    bool penetration_predicate(const Particle& p1, const Particle& p2){
        return (p1.cGet() - c).length() > (p2.cGet() - c).length();
    }

    const Vec& cGet() const{
        return c;
    }
};