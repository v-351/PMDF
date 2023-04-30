#include "../include/init.hpp"
#include "../include/solver.hpp"

int main(int argc, char* argv[]){
	MPIw mpi(&argc, &argv);
	Init init("../taskname.txt");
	Solver solver(init, mpi);
	solver.run();
	solver.results();
	return 0;
};
