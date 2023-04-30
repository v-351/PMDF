#include "../include/surface.h"
#include "../include/constants.h"
#include <vector>
#include <fstream>
#include <sstream>
#include <iostream>
#include <string>
#include <string.h>
#include <tuple>
#include <random>
using namespace std;

// interface methods
int Surface::get_nelem() { 
	return nelem;
}
std::string Surface::get_stlname() {
	return stlfile;
}
std::vector<int> Surface::get_xlink(int i) {
	std::vector<int> ii(3,0);
	for (int j = 0; j < 3; j++) { ii[j] = xlink[i][j]; }
	return ii;
}
std::vector<double> Surface::get_xmin() {
	return xmin;
};
std::vector<double> Surface::get_xmax() {
	return xmax;
};
std::vector<double> Surface::get_x(int i1) {
	std::vector<double> x1(3,0);
	for (int j = 0; j < 3; j++) { x1[j] = x[i1][j]; }
	return x1;
}
std::vector<double> Surface::get_x3(int i) {
	std::vector<double> x1(9, 0);
	for (int j = 0; j < 3; j++) { x1[j] = x[xlink[i][0]][j]; }
	for (int j = 0; j < 3; j++) { x1[3+j] = x[xlink[i][1]][j]; }
	for (int j = 0; j < 3; j++) { x1[6+j] = x[xlink[i][2]][j]; }
	return x1;
}
std::vector<double> Surface::get_n(int i) {
	std::vector<double> x1(3, 0);
	for (int j = 0; j < 3; j++) { x1[j] = n[i][j]; }
	return x1;
}
std::vector<double> Surface::get_seq(int i) {
	std::vector<double> s(4, 0);
	for (int j = 0; j < 3; j++) { s[j] = n[i][j]; }
	s[3] = seq[i];
	return s;
}
double Surface::get_surf(int i) {
	return surf[i];
}
double Surface::get_volume() {
	return volume;
}
double Surface::get_sarea() {
	return sarea;
}

// main method preparing information about surface
Surface::Surface(std::string name, SurfaceKind kind_, int address_)
{
	nelem = 0;
	npoints = 0;
	read_stl_file(name);
	//cout << "Calculate surface data" << endl;
	calc_stl_info();
	kind = kind_;
	address = address_;
}

// method read *.stl type file
void Surface::read_stl_file(std::string name)
{
	ifstream fstl;
	stlfile = name;
	//cout << "Read stl file " << name << endl;
	fstl.open(name);
	std::string line;
	const char * value = NULL;
	std::vector<double> x1;

//	check file exist
	if (!fstl.is_open()) {
		cout << "ERROR: File " << name << " does not exist" << endl;
		return;
	}

	int icou = 0;
	while (fstl) 	{ 
		std::getline(fstl, line);

		if (icou == 4) { icou = 0; }
		if (icou == 0) { value = "normal"; }
		if (icou > 0 && icou<4) { value = "vertex"; }

		if (line.find(value) != std::string::npos)	{
			icou += 1;
			x1 = find_double_vec(line);
			if (value == "normal") { 
				nelem += 1;
				n.push_back(x1);
				xlink.push_back(vector<int>(3));
			}
			if (value == "vertex") { 
				npoints += 1;
				for (int j = 0; j < 3; j++) { x1[j] = x1[j]; }
				x.push_back(x1);
				xlink[nelem - 1][icou - 2] = npoints - 1;
			}
		}
	}
	x1.clear();
	fstl.close();
	//cout << "Surface " << name << " read complete" << " Nelem=" << nelem << " Npoints=" << npoints << endl;
}

// method reads xyz coordinates from line
void Surface::read_stl_xyz(std::string line, double * x)
{
	char * cline = const_cast<char*>(line.c_str());
	const char * delim = " \t\n";
	char * context = NULL;
	char * tok = NULL;
	char * check;
	double val;

	tok = strtok_r(cline, delim, &context);
	int icou = 0;
	while (tok != NULL) {
		val = strtod(tok, &check);
		if (check != tok) 	{
			if (icou < sizeof(x)) 	{
				x[icou] = val;
				icou += 1;
			}
			else {cout << "number of reading data more than expected" << endl;}
		}
		tok = strtok_r(NULL, delim, &context);
	}
}

// method calculates information about surface elements
void Surface::calc_stl_info()
{
	surf.reserve(nelem);
	seq.reserve(nelem);
	icell_min.reserve(3 * nelem);
	icell_max.reserve(3 * nelem);

	double * side = new double [3];
	std::vector<double> nloc(3,0);
	std::vector<int> imin(3,0), imax(3,0), xlt(3,0);
	volume = 0;
	sarea = 0;
	double sarea1;

	for (int i = 0; i < nelem; i++) {
//      calculate surface are of each element
		side[0] = calc_dist(xlink[i][0], xlink[i][1]);
		side[1] = calc_dist(xlink[i][1], xlink[i][2]);
		side[2] = calc_dist(xlink[i][2], xlink[i][0]);
		if (side[0] < side[1] && side[0] < side[2]) { // first element with lowes angle
			for (int j = 0; j < 3; j++) { xlt[j] = xlink[i][j]; }
			xlink[i][0] = xlt[2];
			xlink[i][1] = xlt[0];
			xlink[i][2] = xlt[1];
		} else if (side[2] < side[0] && side[2] < side[1]) {
			for (int j = 0; j < 3; j++) { xlt[j] = xlink[i][j]; }
			xlink[i][0] = xlt[1];
			xlink[i][1] = xlt[2];
			xlink[i][2] = xlt[0];
		}

		sarea1 = calc_surf3(side);
		surf.push_back(sarea1);
//      calculate total surface area
		sarea += sarea1;
//      calculate volme
		volume += (1.0 / 6.0)* \
			   (x[xlink[i][1]][0] * x[xlink[i][2]][1] * x[xlink[i][0]][2] + \
				x[xlink[i][2]][0] * x[xlink[i][0]][1] * x[xlink[i][1]][2] + \
				x[xlink[i][0]][0] * x[xlink[i][1]][1] * x[xlink[i][2]][2] - \
				x[xlink[i][2]][0] * x[xlink[i][1]][1] * x[xlink[i][0]][2] - \
				x[xlink[i][0]][0] * x[xlink[i][2]][1] * x[xlink[i][1]][2] - \
				x[xlink[i][1]][0] * x[xlink[i][0]][1] * x[xlink[i][2]][2]);
//		calculate surface equation parameters and normalize normal vector
		for (int j = 0; j < 3; j++) { nloc[j] = n[i][j]; }
		double norm = vecmod(nloc, 3);
		for (int j = 0; j < 3; j++) { n[i][j] /= norm; }
		seq.push_back( -n[i][0] * x[xlink[i][0]][0] - \
					    n[i][1] * x[xlink[i][0]][1] - \
						n[i][2] * x[xlink[i][0]][2]);
//      update mesh indicies
//		std::tie(imin,imax) = find_elem_indx(i);
//		icell_min.push_back(imin);
//		icell_max.push_back(imax);
	}
	//cout << "Surface volume: " << volume << endl;
	xmin.assign(3, 1e10);
	xmax.assign(3, -1e10);
	for (int i = 0; i < npoints; i++) {
		for (int j = 0; j < 3; j++) {
			if (x[i][j] < xmin[j]) { xmin[j] = x[i][j]; }
			if (x[i][j] > xmax[j]) { xmax[j] = x[i][j]; }
		}
	}
	delete[] side;
}

// find double vector in string
std::vector<double> Surface::find_double_vec(std::string line)
{
	std::vector<std::string> str;
	std::string delim = " \t\n";
	str = split_string(line, delim);

	double val;
	std::vector<double> vec;
	char * cline;
	char * check;
	for (unsigned int i = 0; i < str.size(); i++) {
		cline = const_cast<char*>(str[i].c_str());
		val = strtod(cline, &check);
		if (check != str[i]) { vec.push_back(val); }
	}
	return vec;
}

// returns module of the vector v
double Surface::vecmod(std::vector<double> v, int n)
{
	double mod = 0.0;
	if (n == -1) {
		for (unsigned int i = 0; i < v.size(); i++) { mod += v[i] * v[i]; }
	}
	else {
		for (int i = 0; i < n; i++) { mod += v[i] * v[i]; }
	}
	return sqrt(mod);
}

// multiply all elements of vector v1 on multiplier mult
void Surface::vecmult(std::vector<double> & v1, double mult)
{
	for (unsigned int i = 0; i < v1.size(); i++) { v1[i] += v1[i] * mult; }
}

// method returns distance between vertecies i1 and i2
double Surface::calc_dist(int & i1, int & i2)
{
	return sqrt((x[i1][0] - x[i2][0])*(x[i1][0] - x[i2][0]) + \
				(x[i1][1] - x[i2][1])*(x[i1][1] - x[i2][1]) + \
				(x[i1][2] - x[i2][2])*(x[i1][2] - x[i2][2]));
}

// method returns surface area of triangle with sides side[3]
double Surface::calc_surf3(double * side)
{
	double per = 0;
	for (int i = 0; i < 3; i++) { per += side[i]; }
	per *= 0.5;
	return sqrt(per*(per-side[0])*(per - side[1])*(per - side[2]));
}

// splits string line on string vector using delim
std::vector<std::string> Surface::split_string(std::string line, std::string delim)
{
	char * cline = const_cast<char*>(line.c_str());
	const char * d = delim.c_str();
	char * context = NULL;
	char * tok = NULL;
	std::vector<std::string> str;

	tok = strtok_r(cline, d, &context);
	int icou = 0;
	while (tok != NULL) {
		str.push_back(tok);
		tok = strtok_r(NULL, d, &context);
	}
	return str;
}

// checks vector x1->x2 cross surface and return cross coo
std::tuple<std::vector<double>, int> Surface::check_cross(std::vector<double> & x1, std::vector<double> & x2, int iskip)
{
	std::vector<double> xc(3, 0);
	std::vector<double> xcl;
	std::vector<double> vref;
	std::vector<double> seqi;
	std::vector<int> indxs(3, 0);
	std::vector<int> indxe(3, -1);
	int icou;

	int ielem = -1;
	double ds,dd;
	double dsmin = dbig10;
	double ddmin = dbig10;
	int i = 0;
	for (i = 0; i < nelem; i++) {
		if (i == iskip) { continue; }
		seqi = get_seq(i);  // surface equation

		xcl = flat_line_cross(seqi, x1, x2);  // find cross point
		if (xcl.size() == 0) { continue; }

		ds = check_inside_element_spec(i, xcl);
		dd = vecmod(getvec(x1, xcl));
		if (dd < ddmin && ds < err6) {
			ddmin = dd;
			ielem = i;
			for (int j = 0; j < 3; j++) { xc[j] = xcl[j]; }
		}
	}
	if (ielem >= 0) { return std::make_tuple(xc, ielem); }
	else { 
		xc.clear();
		return std::make_tuple(xc, -1);
	}
}

//=====================================
// geomerty operations

// returns cross point between LINE x1->x2 and surface seq(4), returns 0 vector in case of error
std::vector<double> Surface::flat_line_cross(std::vector<double> & seq, std::vector<double> & x1, std::vector<double> & x2)
{
	std::vector<double> xc = flat_ray_cross(seq, x1, x2);
	if (xc.size() == 0) { return xc; }
	double ndist1 = normdist(seq, x1);
	double ndist2 = normdist(seq, x2);
	if (isign_dbl(ndist1) == isign_dbl(ndist2)) { xc.assign(0, 0); }
	return xc;
}

// returns cross point between RAY x1(start)->x2 and surface seq(4), returns 0 vector in case of error
std::vector<double> Surface::flat_ray_cross(std::vector<double> & seq, std::vector<double> & x1, std::vector<double> & x2)
{
	std::vector<double> v(3, 0);
	v = getvec(x1, x2);
	double cosa = veccosa(seq, v, 3);
	v = vecnorm(v);
	double ndist1 = normdist(seq, x1);
	double ndist2 = normdist(seq, x2);
	if (abs(ndist2) > abs(ndist1) && isign_dbl(ndist1) == isign_dbl(ndist2)) { return std::vector<double>(0); }
	double vlen = abs(ndist1) / abs(cosa);
	if (abs(cosa) < err6) { return std::vector<double>(0); }
	for (int i = 0; i < 3; i++) { v[i] = x1[i] + vlen*v[i]; }
	return v;
}

// returns normal distance between surface seq and point x
double Surface::normdist(std::vector<double> & seq, std::vector<double> & x)
{
	return seq[0] * x[0] + seq[1] * x[1] + seq[2] * x[2] + seq[3];
}

// returns cosa between two vctors v1 and v2
double Surface::veccosa(std::vector<double> & v1, std::vector<double> & v2, const int n)
{
	double dotp = 0;
	if (n == -1) {
		if (v1.size() != v2.size()) { return -2; }
		for (unsigned int i = 0; i < v1.size(); i++) { dotp += v1[i] * v2[i]; }
		dotp = dotp / vecmod(v1) / vecmod(v2);
	}
	else {
		if ((int)v1.size() <  n || (int)v2.size() < n) { return -2; }
		for (int i = 0; i < n; i++) { dotp += v1[i] * v2[i]; }
		dotp = dotp / vecmod(v1, n) / vecmod(v2, n);
	}
	return dotp;
}

// returns vector from to points x1->x2
std::vector<double> Surface::getvec(std::vector<double> & x1, std::vector<double> & x2)
{
	if (x1.size() != x2.size()) { return std::vector<double>(0); }
	std::vector<double> v(x1.size());
	for (unsigned int i = 0; i < x1.size(); i++) { v[i] = x2[i] - x1[i]; }
	return v;
}

// returns normalized vector v
std::vector<double> Surface::vecnorm(std::vector<double> v)
{
	double mod = vecmod(v);
	for (unsigned int i = 0; i < v.size(); i++) { v[i] /= mod; }
	return v;
}

// returns dot product for vectors v1 and v2 (optionally one can specify size of dot product)
double Surface::vecdot(std::vector<double> & v1, std::vector<double> & v2, const int n)
{
	double dotp = 0;
	if (n == -1) {
		if (v1.size() != v2.size()) { return -2; }
		for (unsigned int i = 0; i < v1.size(); i++) { dotp += v1[i] * v2[i]; }
	}
	else {
		if ((int)v1.size() <  n || (int)v2.size() < n) { return -2; }
		for (int i = 0; i < n; i++) { dotp += v1[i] * v2[i]; }
	}
	return dotp;
}

// returns sign (int) of double value
/// IMPORTANT: >= - ошибочный пролет при координате равной координате поверхности
int Surface::isign_dbl(double & val)
{
	if (val > 0.0) { return 1; } else { return - 1; }
}

// reflect ray, n must be normalized!!
// return reflected ray (0:2)
std::vector<double> Surface::reflect_ray(std::vector<double> & vec, std::vector<double> & n)
{
	std::vector<double> res(3,0);
	double dot = vecdot(vec, n);
	for (int j = 0; j < 3; j++) { res[j] = vec[j] - 2.0*n[j] * dot; }
	return res;
}

// reflect ray, n must be normalized!!
// return reflected ray (0:2), normal component of reflected ray (3:5), tangential component of reflected ray (6:8)
std::vector<double> Surface::reflect_ray_full(std::vector<double> & vec, std::vector<double> & n)
{
	std::vector<double> res(9,0);
	double dot = vecdot(vec, n);
	for (int j = 0; j < 3; j++) {
		res[j] = vec[j] - 2.0*n[j] * dot;  // reflected ray
		res[j + 3] = -dot*n[j];            // normal component
		res[j + 6] = res[j] - res[j + 3];  // tangential component
	}
	return res;
}

// checks that point x1 are inside triangle i
double Surface::check_inside_element_spec(int & i, std::vector<double> & x1)
{
	//cout << "check_inside_element_spec " << i << endl;
	double v0, v1, v2, v01, v12, v20;
	v0 = vecmod(getvec(x1, x[xlink[i][0]]));
	v1 = vecmod(getvec(x1, x[xlink[i][1]]));
	v2 = vecmod(getvec(x1, x[xlink[i][2]]));
	v01 = vecmod(getvec(x[xlink[i][0]], x[xlink[i][1]]));
	v12 = vecmod(getvec(x[xlink[i][1]], x[xlink[i][2]]));
	v20 = vecmod(getvec(x[xlink[i][2]], x[xlink[i][0]]));

	if (v0 < err10) { return true; }
	if (v1 < err10) { return true; }
	if (v2 < err10) { return true; }

	double ss = get_ss3(v01, v12, v20);
	double ss1 = get_ss3(v01, v0, v1);
	double ss2 = get_ss3(v12, v1, v2);
	double ss3 = get_ss3(v20, v2, v0);

	return ss1 + ss2 + ss3 - ss;
}

// return triangle surface based on 3 sides
double Surface::get_ss3(double & a, double & b, double & c)
{
	//cout << a << "  " << b << "  " << c << endl;
	double p = 0.5*(a + b + c);
	return sqrt(p*(p - a)*(p - b)*(p - c));
}

// rotate vector vec on angle phi (-pi...pi, [rad]) on surface (seq)
std::vector<double> Surface::vecrot_surf(std::vector<double> & vec, std::vector<double> & nn, double & phi)
{
	double sinp = std::sin(phi);
	double cosp = std::sin(phi);

	double m11 = cosp + (1.0-cosp)*nn[0]*nn[0];
	double m12 = (1.0-cosp)*nn[0]*nn[1] - sinp*nn[2];
	double m13 = (1.0-cosp)*nn[0]*nn[2] + sinp*nn[1];

	double m21 = (1.0-cosp)*nn[1]*nn[0] + sinp*nn[2];
	double m22 = cosp + (1.0-cosp)*nn[1]*nn[1];
	double m23 = (1.0-cosp)*nn[1]*nn[2] - sinp*nn[0];

	double m31 = (1.0-cosp)*nn[2]*nn[0] - sinp*nn[1];
    double m32 = (1.0-cosp)*nn[2]*nn[1] + sinp*nn[0];
	double m33 = cosp + (1.0-cosp)*nn[2]*nn[2];

	std::vector<double> vres(3,0);
	vres[0] = m11*vec[0] + m12*vec[1] + m13*vec[2];
	vres[1] = m21*vec[0] + m22*vec[1] + m23*vec[2];
	vres[2] = m31*vec[0] + m32*vec[1] + m33*vec[2];

	return vres;
}