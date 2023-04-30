#pragma once
#include <fstream>
#include <iostream>
#include <vector>

enum SurfaceKind{
	SURFACEKIND_WALL,
	SURFACEKIND_MPI,
	SURFACEKIND_INLET,
	SURFACEKIND_OUTLET
};

class Surface
{
private:
// variables
	std::string stlfile; // stl file related to this surface class
	int nelem;  // number of elements
	int npoints;  // number of points
	double volume;  // volume of closed stl surface [m3]
	double sarea;  // total surface area of stl surface [m2]
	std::vector<double> xmin;  // min x,y,z coordinates [3]
	std::vector<double> xmax;  // max x,y,z coordinates [3]
	std::vector<std::vector<int>> xlink;  // indicies of triangle points [nelem][3]
	std::vector<std::vector<int>> icell_min;  // min indicies of triangle points [nelem][3]
	std::vector<std::vector<int>> icell_max;  // min indicies of triangle points [nelem][3]
	std::vector<std::vector<double>> x;  // coordinates of points [npoints][3]
	std::vector<std::vector<double>> n;  // normal vector for element [nelem][3]
	std::vector<double> seq;  // 4th value in surface equation, first 3 in n [nelem]
	std::vector<double> surf;  // surface area of element [nelem]
	SurfaceKind kind;
	int address = -1;
// functions
	void read_stl_file(std::string name);  // read stl file
	void read_stl_xyz(std::string line, double * x);  // read first 3 numbers in string line
	void calc_stl_info();  // calculate information about stl, volume, sarea, surface equations

//  
    std::vector<double> find_double_vec(std::string line); // find double vector in string
	double calc_dist(int & i1, int & i2); // method returns distance between vertecies i1 and i2
    double calc_surf3(double * side); // method returns surface area of triangle with sides side[3]
    std::vector<std::string> split_string(std::string line, std::string delim); // splits string line on string vector using delim
	double check_inside_element_spec(int & i, std::vector<double> & x1);  // checks that point x1 are inside triangle i
//=====================================
// geomerty operations
    std::vector<double> flat_line_cross(std::vector<double> & seq, std::vector<double> & x1, std::vector<double> & x2); // returns cross point between LINE x1->x2 and surface seq(4), returns 0 vector in case of error
    std::vector<double> flat_ray_cross(std::vector<double> & seq, std::vector<double> & x1, std::vector<double> & x2); // returns cross point between RAY x1->x2 and surface seq(4), returns 0 vector in case of error
	int isign_dbl(double & val); // returns sign (int) of double value
    std::vector<double> reflect_ray(std::vector<double> & vec, std::vector<double> & n); // reflect ray, n must be normalized!! // return reflected ray (0:2)
    double get_ss3(double & a, double & b, double & c); // return triangle surface based on 3 sides

public:
	Surface(std::string name, SurfaceKind kind_, int address_ = 0); // constructor
	Surface() {};

	int& getAddress(){
		return address;
	}

    std::tuple<std::vector<double>, int> check_cross(std::vector<double> & x1, std::vector<double> & x2, int iskip=-1); // checks vector x1->x2 cross surface and return cross coo
    std::vector<double> reflect_ray_full(std::vector<double> & vec, std::vector<double> & n); // reflect ray, n must be normalized!! // return reflected ray (0:2), normal component of reflected ray (3:5), tangential component of reflected ray (6:8)
	// vector operations
	std::vector<double> getvec(std::vector<double> & x1, std::vector<double> & x2); // returns vector from to points x1->x2
    void vecmult(std::vector<double> & v1, double mult); // multiply all elements of vector v1 on multiplier mult
	double vecdot(std::vector<double> & v1, std::vector<double> & v2, const int n = -1); // returns dot product for vectors v1 and v2 (optionally one can specify size of dot product)
    double normdist(std::vector<double> & seq, std::vector<double> & x); // returns normal distance between surface seq and point x
    double veccosa(std::vector<double> & v1, std::vector<double> & v2, const int n = -1); // returns cosa between two vctors v1 and v2
	std::vector<double> vecnorm(std::vector<double> v); // returns normalized vector v
	double vecmod(std::vector<double> v, int n=-1); // returns module of the vector v
	std::vector<double> vecrot_surf(std::vector<double> & vec, std::vector<double> & seq, double & phi); // rotate vector vec on angle phi (-pi...pi, [rad]) on surface (seq)
//==================================
//  interface
	int get_nelem();
	double get_volume();
	double get_sarea();
	std::string get_stlname();
	std::vector<int> get_xlink(int i);
	std::vector<double> get_xmin();
	std::vector<double> get_xmax();
	std::vector<double> get_x(int i);
	std::vector<double> get_x3(int i);
	std::vector<double> get_n(int i);
	std::vector<double> get_seq(int i);
	double get_surf(int i);
	SurfaceKind get_kind(){ return kind; }
};