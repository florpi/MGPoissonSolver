#ifndef grid_H
#define grid_H

#include "multigrid.h"

class Grid{
	private:
		const double L; // Lenght of the box to simulate
		const double G; // Gravitational constant, default set to 1.
		int step;
		vector<Particle> particles; // Grid contains array of particles belonging to the grid

	public:
		Grid(int MAXNPART, const double L );
		vector2d<double> rhs; // right hand side in equation to solve, for us  density 
		vector2d<double> lhs;// left hand side in equation to solve, for us potential
		vector2d<double> residual;
		vector2d<double> error;

		int  MAXGRID; // Grid size
		const double h; // Cell separation 
		vector2d<double> fx; // force, x component
		vector2d<double> fy; // force, y component
		void add_particle(Particle part);
		void compute_density(); // Assigns density to the grid using CIC 
		double get_density(int, int);
		void compute_force(); // Computes fx and fy
		vector2d<double> get_fx();
		vector2d<double> get_fy();
};

#endif
