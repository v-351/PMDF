#include "../include/metrics.hpp"
#include "../include/inlet.hpp"

#include <fstream>
#include <algorithm>
#include <functional>

LogSingle& LogSingle::log(const Particle& p, const Dukowicz& duk, double time){
    v.push_back({time, p.mGet(), p.nGet(), p.dGet(), p.cGet()[1], p.all_evap, p.TGet(), duk.q, duk.Re}); //
    return *this;
}

LogSingle& LogSingle::write_log(std::string filepath){
    std::ofstream file(filepath);
    file << ("time,m,n,d,pen,evap,T,q,Re,x,y,z\n");
    for(const std::vector<double>& i: v){
        for(const double& d: i){
            file << d;
            if(&d != &*i.rbegin()){
                file << ",";
            }
        }
        file << "\n";
    }
    file.close();
    return *this;
}

LogSingle& LogSingle::log_coord(const Particle& p, double time){
    time_coord.push_back(std::make_pair(time, p.cGet()));
    return *this;
}

LogSingle& LogSingle::log_value(const double& p, double time){
    time_value.push_back(std::make_pair(time, p));
    return *this;
}

LogSingle& LogSingle::write_log_coord(std::string filepath){
    std::ofstream file(filepath);
    file << ("time,x,y,z\n");
    for(const std::pair<double,Vec>& p: time_coord){
        file << p.first << "," << p.second.print() << "\n";
    }
    file.close();
    return *this;
}

LogSingle& LogSingle::log_pen(std::vector<Particle>& particles, const Inlet& inlet, double time){
    double mass = Particle::all_inj - Particle::all_evap;
    auto cmp = std::bind<bool>(&Inlet::penetration_predicate, inlet, std::placeholders::_1, std::placeholders::_2);
    std::sort(particles.begin(), particles.end(), cmp);
    double sum_mass = 0;
    double pen = 0;
    mass *= 0.95;
    std::vector<Particle>::reverse_iterator riter = particles.rbegin();
    for(riter; riter < particles.rend(); riter++){
        if(sum_mass >= mass){
            break;
        }
        sum_mass += riter->mGet() * riter->nGet();
        pen = (riter->cGet() - inlet.cGet()).length();
    }

    time_pen.emplace_back(std::make_pair(time, pen));

    return *this;
}

LogSingle& LogSingle::write_pen(std::string filepath){
    std::ofstream file(filepath);
    file << ("time,pen\n");
    for(const std::pair<double,double>& p: time_pen){
        file << p.first << "," << p.second << "\n";
    }
    file.close();
    return *this;
}

LogSingle& LogSingle::write_value(std::string filepath){
    std::ofstream file(filepath);
    file << ("time,value\n");
    for(const std::pair<double,double>& p: time_value){
        file << p.first << "," << p.second << "\n";
    }
    file.close();
    return *this;
}

LogSingle& LogSingle::write_particle(const std::vector<Particle>& particles, std::string filepath){
    std::ofstream file(filepath);
    file << ("x,y,z,vx,vy,vz,tvmod,d,m,n,T\n");
    for(const Particle& p: particles){
        file << p.cGet().print() << "," << p.vGet().print()<< "," << p.v_pulseGet().length() << "," << 
            p.dGet() << "," << p.mGet()<< "," << p.nGet()<< "," << p.TGet() << "\n";
    }
    file.close();
    return *this;
}

void LogSingle::write_cells(std::vector<Cell>& cells, std::string filepath){
    std::ofstream file(filepath);
    file << ("x,y,z,moment_x,moment_y,moment_z,m,q,mark\n");
    for(const Cell& c: cells){
        file << c.cGet().print() << "," << c.momentumGet().print() << "," << c.massGet() << "," << c.heatGet() << "," << c.markerGet() << "\n";
    }
    file.close();
}

void LogSingle::set_values_names(std::string s){
    values_names = s;
}

void LogSingle::write_to_values(std::vector<double> vec){
    values.push_back(vec);
}

void LogSingle::write_values(std::string filepath){
    std::ofstream file(filepath);
    file << (values_names + "\n");
    for(const std::vector<double>& i: values){
        for(const double& d: i){
            file << d;
            if(&d != &*i.rbegin()){
                file << ",";
            }
        }
        file << "\n";
    }
    file.close();
}

