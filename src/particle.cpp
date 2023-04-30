#include "../include/particle.hpp"

/** Испарение. 
 * Применение изменений в температуре и массе парцела.
 * Проверка корректности значений. Изменение диаметра в соответствии с
 * изменением массы. 
 */
void Particle::evaporation(const double& dT, const double& dm){
    T += dT;
    m += dm;
    all_evap -= dm * n;
    if(m > 0){
        d = std::pow(6 * m / props->ro[T] / pi, 1./3);
    }else{
        outOfBounds = true;
    }

    if(m <= 0 || d <= 0){
        outOfBounds = true;
    }
}

/** Дробление. 
 * Применение изменения диаметра парцела.
 * Проверка корректности значений. Изменение массы и количества
 * частиц в парцеле в соответствии с изменением диаметра частиц.
 */
void Particle::breakUp(const double& dd){
    if(dd > 0){
        double d_old = d;
        d -= dd;
        m = props->ro[T] * pi * d * d * d / 6;
        n = std::pow(d_old / d, 3) * n;
    }
    if(m <= 0 || d <= 0){
        outOfBounds = true;
    }
}

std::string Particle::properties() const{
    return std::to_string(T) + "," + std::to_string(m) + "," + std::to_string(d)+ "," + std::to_string(n);
}

std::string Particle::propertiesNames() const{
    return "T,m,d,n";
}

bool Particle::operator<(const Particle& p) const{
    return c.length() < p.c.length();
}

double Particle::all_evap = 0;
double Particle::all_inj = 0;