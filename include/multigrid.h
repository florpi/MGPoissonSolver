#ifndef multigrid_H 
#define multigrid_H 

#include <math.h>
#include "vector3d.cpp"
#include "particle.h"
#include "grid.h"
#include "omp.h"
class Multigrid{
	private:
		int gridlevel;
		int n;
		int currentstep;
		const int maxlevel;

	public:
	vector3d<double> left;
	vector3d<double> right;
	vector3d<double> residual;
	
	vector<Grid> grids;
	Multigrid(const int maxlevel);
	void Initial_conditions(Grid mother);

	int get_currentstep();
	void restrict(); // restricts a given quantity to a finer grid via interpolation
		void prolong();// prolongs a given quantity to a coarser grid via interpolation
		void gauss(int n_iters); // Performs n_iters gauss_seidel steps
		void gauss_omp(int n_iters); // Performs n_iters gauss_seidel steps
		void compute_residual();
		void result(int n_iters0);
		void vcycle(int n_iters);
		//vector3d<double> get_left();

};
vector<int> applyBC(int i, int j, int MAXGRID); // Applies periodic boundary conditions
vector<int> applyBC_doublestep(int i, int j, int MAXGRID); // Applies periodic boundary conditions


void print(vector3d<double> &matrix, int, int, int);
#endif
