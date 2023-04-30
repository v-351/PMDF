#include "../include/property.hpp"

const double& Property::operator[](const double& T) const{
    // Обработка выхода за границы диапазона через возврат граничных значений.
    if(T <= Tmin){
        return *pr.begin();
    }
    if(T >= Tmax){
        return *pr.rbegin();
    }
    // Возврат значения
    return pr[int(T-Tmin)];
};

void Property::readCSV(std::string filepath, char sep /* = ',' */, int Tstep_/* = 1 */ ){
    Tstep = Tstep_;
    pr.reserve(1000);
    double V, dV, Vlast, a;
    int T, dT, Tlast;
    int pos;
    std::ifstream file(filepath);
    std::string line;
    std::getline(file, line); //skip header
    std::getline(file, line);
    pos = line.find(sep);
    Tlast = (int)std::stod(line.substr(0, pos));
    Vlast = std::stod(line.substr(pos+1));
    Tmin = Tlast;
    pr.push_back(Vlast);
    while(std::getline(file, line)){
        pos = line.find(sep);
        T = (int)std::stod(line.substr(0, pos));
        V = std::stod(line.substr(pos+1));
        dT = T - Tlast;
        dV = V - Vlast;
        if(dT == 1){
            pr.push_back(V);
        }else{
            a = dV / dT;
            for(int Tincr = 1; Tincr < dT; Tincr += Tstep){
                Vlast += a * Tstep;
                pr.push_back(Vlast);
            }
            pr.push_back(V);
        }
        Tlast = T;
        Vlast = V;
    }
    file.close();
    Tmax = Tmin + pr.size()*Tstep;
};

/// Контруктор. Требует путь до папки с перечнем файлов.
ParticleSubstance::ParticleSubstance(std::string folder){
    p_sat.readCSV(folder + "p_sat");
    ro.readCSV(folder + "ro");
    L.readCSV(folder + "L");
    sigma.readCSV(folder + "sigma");
    lambda.readCSV(folder + "lambda");
    lambda_vap.readCSV(folder + "lambda_vap");
    cp.readCSV(folder + "cp");
    cp_vap.readCSV(folder + "cp_vap");
    dyn_vis.readCSV(folder + "dyn_vis");
    dyn_vis_vap.readCSV(folder + "dyn_vis_vap");

    std::ifstream file(folder + "M");
    file >> M;
    file.close();
};

void Table::parse(std::string filepath){
    std::ifstream file(filepath);
    std::string line;
    std::getline(file, line); //skip header
    char sep = ' ';
    double pos;
    while(std::getline(file, line)){
        pos = line.find(sep);
        index = std::stod(line.substr(0, pos));
        if(index == 0){
            index = 1e-10;
        }
        value = std::stod(line.substr(pos+1));    
        if(value == 0){
            value = 1e-10;
        }
        valmap[index] = value;
    }
    file.close();
};

double Table::k(const double& index) const{
    return valmap.lower_bound(index)->second;
};