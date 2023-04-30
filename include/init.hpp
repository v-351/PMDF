#pragma once

#include <string>
#include <filesystem> // C++17
#include <unordered_map>

#include "../nlohmann/json.hpp"
#include "surface.h"
#include "vec.hpp"
#include "cell.hpp"

class Init{
public:
	std::string taskname;

	double dt_global;
	double fq;
	double total_time;
	double connection_range;
	bool single_particle_flag;
	bool single_cell_flag;
	bool evap_flag;
	bool brup_flag;
	bool turbulent_dispersion_flag;
	bool wall_flag;
	std::string substance_name;
	std::string stl_dir;
	std::string cells_dir;

    int inlets_number;
    Vec inlet_coord;
    Vec inlet_n;
    double inlet_alpha;
    double inlet_d;
    double inlet_G;
    double inlet_T;

    double duk_C1;

    double wave_C1, wave_C2;

	Init(std::string taskpath){
		std::ifstream file(taskpath); // Чтение имени задачи из заданного файла
		file >> taskname;
		file.close();
		
        std::ifstream conf_file("../configs/" + taskname + "-config.json"); // Чтение конфигурационного файла
        nlohmann::json json = nlohmann::json::parse(conf_file);
        conf_file.close();

        dt_global = json["config"]["dt"]; // Глобальный шаг по времени
        fq = 1/dt_global; 
        total_time = json["config"]["t"]; // Общее время расчета
        connection_range = json["config"]["connection_range"]; // Радиус области взаимодействия частицы с ячейками

        inlet_coord = Vec(json["inlet"]["c"][0], json["inlet"]["c"][1] , json["inlet"]["c"][2]);
        inlet_n = Vec(json["inlet"]["n"][0], json["inlet"]["n"][1], json["inlet"]["n"][2]);
        inlet_alpha = json["inlet"]["alpha"];
        inlet_d = json["inlet"]["d"];
        inlet_G = json["inlet"]["G"];
        inlet_T = json["inlet"]["T"];
        inlets_number = json["inlets_number"];
        substance_name = json["inlet"]["substance"];

        single_particle_flag = json["use_single_particle"]; // Флаг режима единичной частицы
        single_cell_flag = json["use_single_cell"]; // Флаг режима единичной ячейки

        evap_flag = json["evap"]["engage"]; // Флаг модели испарения
        brup_flag = json["brup"]["engage"]; // Флаг модели дробления
        turbulent_dispersion_flag = json["turbulent_dispersion"]["engage"]; // Флаг модели турбулентной дисперсии
        wall_flag = json["walls"]; // Флаг модели взаимедействия со стенками

        duk_C1 = json["evap"]["C1"];
        wave_C1 = json["brup"]["C1"];
        wave_C2 = json["brup"]["C2"];

        stl_dir = json["stl_dir"];
        stl_dir = "../configs/" + stl_dir;
        cells_dir = json["cells_dir"];
        cells_dir = "../configs/" + cells_dir;
		
	}

    void readSurface(std::vector<Surface>& surface, int mpi_rank, int mpi_size) const{
        // Заполнение массива повехностей на основании названия файлов
        for(auto& f: std::filesystem::directory_iterator(stl_dir)){ // C++17
            std::string fname = f.path().filename().string(); // Выделение имени файла из пути до файла
            // Проверка наличия в имени файла номера данного процесса
            if(fname.find(std::to_string(mpi_rank)) != std::string::npos ){ 
                if(fname.find("wall") != std::string::npos){
                    surface.emplace_back(f.path(), SURFACEKIND_WALL);
                    std::cout << f.path() << " | wall | on rank " << mpi_rank << "\n";
                }
                if(fname.find("inlet") != std::string::npos){
                    surface.emplace_back(f.path(), SURFACEKIND_INLET);
                    std::cout << f.path() << " | inlet | on rank " << mpi_rank << "\n";
                }
                if(fname.find("outlet") != std::string::npos){
                    surface.emplace_back(f.path(), SURFACEKIND_OUTLET);
                    std::cout << f.path() << " | outlet | on rank " << mpi_rank << "\n";
                }
                if(fname.find("mpi") != std::string::npos){
                    int from = fname[8] - '0';
                    int to = fname[10] - '0';
                    if(to == mpi_rank){
                        std::swap(from, to);
                    }
                    if(to < mpi_size){
                        surface.emplace_back(f.path(), SURFACEKIND_MPI, to);
                        std::cout << f.path() << " | MPI | from " << mpi_rank << " to " << to << "\n";
                    }else{
                        surface.emplace_back(f.path(), SURFACEKIND_OUTLET);
                        std::cout << f.path() << " | outlet | on rank " << mpi_rank << "\n";
                    }
                }
            }	
        }
        surface.shrink_to_fit(); // Сжатие размеров массива до занятого объема ( capacity = size; )
    }

    void readCells(std::vector<Cell>& cells, int mpi_rank) const {
        for(auto& f: std::filesystem::directory_iterator(cells_dir)){ // C++17
            if(f.path().filename().string()[0] == ('0'+ mpi_rank)){
                std::ifstream file(f.path());
                std::vector<double> s;
                double val;
                std::stringstream ss;
                std::string buffer;
                std::getline(file, buffer);
                while(std::getline(file, buffer)){
                    ss.clear();
                    ss << buffer;
                    s.clear();
                    while(ss >> val){
                        s.push_back(val);
                    }
                    cells.emplace_back(
                        Vec(s[0], s[1], s[2]),
                        Vec(s[8], s[9], s[10]),
                        s[3], s[11], s[12], s[13]);
                }
                file.close();
                std::cout << "cells | on rank " << mpi_rank << " file: " << f.path() << "\n";
            }
	    }
        cells.shrink_to_fit();
    }

    /*
    void constructInlets(std::vector<Inlet>& inlets){
        for(int i = 0; i < inlets_number; i++){
            std::string name = "inlet" + std::to_string(i);
            Vec inlet_coord = Vec(json[name]["c"][0], json["inlet"]["c"][1] , json["inlet"]["c"][2]);
            inlets.emplace_back()


            // Флаг включения форсунки. Ставится исходя из наличия хотя бы одной ячейки расчетной сетки в окрестности координаты форсунки размером (connection_range/2).
            bool inlet_flag = std::any_of(cells.begin(), cells.end(), [inlet_coord, connection_range](Cell& c){return ((inlet_coord - c.cGet()).length() < (connection_range/2));});

            Inlet inlet(Vec(json["inlet"]["c"][0], json["inlet"]["c"][1] , json["inlet"]["c"][2] ), 
                json["inlet"]["alpha"], json["inlet"]["d"], json["inlet"]["G"],
                json["inlet"]["T"], particle_substance, fq, inlet_flag,
                {json["inlet"]["n"][0], json["inlet"]["n"][1], json["inlet"]["n"][2]}
            );
        }
    }
    */

    bool inletFlag(std::vector<Cell>& cells){
        // Флаг включения форсунки. Ставится исходя из наличия хотя бы одной ячейки расчетной сетки в окрестности координаты форсунки размером (connection_range/2).
	    return std::any_of(cells.begin(), cells.end(), [this](Cell& c){return ((this->inlet_coord - c.cGet()).length() < (this->connection_range/2));});
    }
};