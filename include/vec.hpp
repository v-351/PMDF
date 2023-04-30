#pragma once

#include <cmath>
#include <vector>
#include <array>
#include <string>
#include <tuple>

const double pi = 3.14159; ///< Число Пи.

/** Математический вектор. 
 * Класс моделирующий математический трехмерный вектор. Используется для хранения скоростей и координат.
 * Математические операторы (+, -, *, /) перегружены для обеспечения корректного математического поведения. 
 * 
 */
class Vec {
private:
    double x, y, z; ///< Координаты
public:
    /** Конструктор по умолчанию. 
     * Создает нулевой вектор (0, 0, 0).
     */
    Vec(): x(0), y(0), z(0) {};

    /** Конструктор.
     * @param x_ Координата по оси x.
     * @param y_ Координата по оси y.
     * @param z_ Координата по оси z.
     */
    Vec(double x_, double y_, double z_): x(x_), y(y_), z(z_) {}

    void set(double x_, double y_, double z_){
        x = x_; 
        y = y_; 
        z = z_; 
    }

    void clear(){
        x = 0;
        y = 0;
        z = 0;
    }
    
    /// Оператор доступа. [x,y,z] = [0,1,2]
    const double& operator[](const int& axis) const{

        if(axis == 0){
            return x;
        }
        if(axis == 1){
            return y;
        }
        if(axis == 2){
            return z;
        }
        return x;
    } 

    /// Векторное сложение.
    Vec operator+(const Vec& v) const{
        return Vec(x + v.x, y + v.y, z + v.z);
    }

    friend Vec operator+(const double& n, const Vec& v){
        return Vec(n + v.x, n + v.y, n + v.z);
    }

    void operator+=(const Vec& v){
        x += v.x;
        y += v.y;
        z += v.z;
    }

    /// Векторное вычитание.
    Vec operator-(const Vec& v) const{
        return Vec(x - v.x, y - v.y, z - v.z);
    }

    Vec operator-(const double& i) const{
        return Vec(x - i, y - i, z - i);
    }

    friend Vec operator-(const double& n, const Vec& v){
        return Vec(n - v.x, n - v.y, n - v.z);
    }

    void operator-=(const Vec& v){
        x -= v.x;
        y -= v.y;
        z -= v.z;
    }
    

    /// Поэлементное произведение.
    Vec operator*(const Vec& p) const{
        return Vec(x * p.x, y * p.y, z * p.z);
    }

    /// Скалярное произведение.
    friend double scalar(const Vec& v1, const Vec& v2) {
        return (v1.x * v2.x + v1.y * v2.y + v1.z * v2.z);
    }

    template <typename T>
    Vec operator*(const T& d) const{
        return Vec(x * d, y * d, z * d);
    }
    
    Vec normalize(){
        double len = length();
        return Vec(x/len, y/len, z/len);
    }
    
    /// Векторное произведение
    friend Vec dot(const Vec& v1, const Vec& v2){
        return Vec(v1[1]*v2[2] - v1[2]*v2[1], v1[2]*v2[0] - v1[0]*v2[2], v1[0]*v2[1] - v1[1]*v2[0]);
    }

    /**
     * @brief Нормаль.
     * 
     * @return Вектор перпендикулярный данному.
     */
    Vec norm() const{
        double a[3] = {x,y,z};
        int imin = 0;
        double min = a[0];
        for(int i = 0; i < 3; i++){
            if(a[i] < min){
                min = a[i];
                imin = i;
            }
        }
        a[imin] = 0;
        if(imin == 0){
            std::swap(a[1],a[2]);
            a[1] *= -1;
        }
        if(imin == 1){
            std::swap(a[0],a[2]);
            a[0] *= -1;
        }
        if(imin == 2){
            std::swap(a[0],a[1]);
            a[0] *= -1;
        }
        return Vec(a[0], a[1], a[2]);
    }
    
    /// Поэлементное деление.
    Vec operator/(const Vec& p) const{
        return Vec(x / p.x, y / p.y, z / p.z);
    }
    template <typename T>
    Vec operator/(const T& d) const{
        return Vec(x / d, y / d, z / d);
    }
    
    template <typename T>
    void operator/=(const T& d) {
        x /= d;
        y /= d; 
        z /= d;
    }

    friend Vec operator/(int left, const Vec& p){
        return Vec(left / p.x, left / p.y, left / p.z);
    }

    /// Оператор сравнения
    bool operator<(const Vec& v) const{
        if(x > v.x) { return false; }
        if(y > v.y) { return false; }
        if(z >= v.z) { return false; }
        return true;
    }

    /// Длина.
    double length() const{
        if(x == 0 && y==0 && z==0){
            return 1e-10;
        }
        return std::sqrt((x*x) + (y*y) + (z*z));
    }

    double length(Vec& c){
        if((c.x - x) == 0 && (c.y - y)==0 && (c.z - z)==0){
            return 1e-10;
        }
        return std::sqrt(((c.x - x)*(c.x - x)) + ((c.y - y)*(c.y - y)) + ((c.z - z)*(c.z - z)));
    }

    /// Поэлементный логарифм.
    friend Vec log(const Vec& c){
        return Vec(std::log(c.x), std::log(c.y), std::log(c.z));
    }

    bool within(const Vec& begin, const Vec& end) const{
        return ( (x >= begin.x && x < end.x) &&
        (y >= begin.y && y < end.y) &&
        (z >= begin.z && z < end.z) );
    }

    /// Сериализация
    void serialize(std::vector<double>& v) const{
        v.push_back(x);
        v.push_back(y);
        v.push_back(z);
    }

    /// Возвращает координату точки строкой из координат разделенных запятой.
    std::string print() const{
        return (std::to_string(x) + "," + std::to_string(y) + "," + std::to_string(z));
    }

    /**
     * @brief Приведение к стандартной библиотеке.
     * 
     * @return std::vector<double> (x,y,z)
     */
    std::vector<double> to_std() const {
        return {x,y,z};
    }
};

/// Векторное умножение
Vec dot(const Vec& v1, const Vec& v2);

/// Скалярное умножение
double scalar(const Vec& v1, const Vec& v2);