#ifndef multigrid_H 
#define multigrid_H 

#include <math.h>
#include "vector3d.cpp"
#include "particle.h"
#include "grid.h"
#include "omp.h"
class Multigrid{
	private:
		int gridlevel; // Current grid, Each grid has different cell size
		int n; 
		const int maxlevel; // Last grid

	public:
		Multigrid(const int maxlevel);
		vector<Grid> grids; // array of grids with different sizes
		void Initial_conditions(Grid mother); // generates initial array of grids 
		void restrict(); // restricts a given quantity to a finer grid via interpolation
		void prolong();// prolongs a given quantity to a coarser grid via interpolation
		void gauss(int n_iters); // Performs n_iters gauss_seidel steps
		void gauss_omp(int n_iters); // Performs n_iters gauss_seidel steps with parallel red black ordering
		void compute_residual();
		void vcycle(int n_iters, bool paral);

};
vector<int> applyBC(int i, int j, int MAXGRID); // Applies periodic boundary conditions

#endif
