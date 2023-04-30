#pragma once

#include "particle.hpp"
#include "cell.hpp"
#include "models.hpp"
#include "inlet.hpp"

class LogSingle{
private:
    std::vector<std::vector<double>> v;
    std::vector<std::pair<double,Vec>> time_coord;
    std::vector<std::pair<double,double>> time_pen;
    std::vector<std::pair<double,double>> time_value;
    std::vector<std::vector<double>> values;
    std::string values_names;

public:
    LogSingle& log(const Particle& p, const Dukowicz& duk, double time);
    LogSingle& log_coord(const Particle& p, double time);
    LogSingle& write_log(std::string filepath);
    LogSingle& write_log_coord(std::string filepath);
    LogSingle& log_pen(std::vector<Particle>& vec, const Inlet& inlet, double time);
    LogSingle& write_pen(std::string filepath);
    LogSingle& write_particle(const std::vector<Particle>& particles, std::string filepath);
    LogSingle& log_value(const double& p, double time);
    LogSingle& write_value(std::string filepath);
    void write_cells(std::vector<Cell>& cells, std::string filepath);
    void write_to_values(std::vector<double> vec);
    void write_values(std::string filepath);
    void set_values_names(std::string s);
};

/**
 * @brief Абстрактный класс для логгера.
 * @details Устанавливает интерфейс для логгеров в виде метода сохранения и метода записи в файл.
 */
class LoggerAbstract{
public:
    virtual void push_back() = 0;
    virtual void writeToCsv(std::string) = 0;
};


class Logger: public LoggerAbstract{
private:
    std::vector<std::pair<double, double>> values;
public:
    Logger(){
        values.reserve(1e3);
    }
    
    virtual void push_back(double v1, double v2){
        values.emplace_back(v1, v2);
    }

    void writeToCsv(std::string filepath){
        std::ofstream file(filepath);
        file << ("time,value\n");
        for(const std::pair<double,double>& v: values){
            file << v.first << "," << v.second << "\n";
        }
        file.close();
    }
};

class LoggerVector: public LoggerAbstract{

};
