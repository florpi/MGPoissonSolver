#ifndef force_H
#define force_H

#include "potential.h"
extern int MAXGRID;
extern double h;
extern int currentstep;
double invert_force(Particle *part, 
		vector2d<double> &fx,
		vector2d<double> &fy);

// Retrieves grid force using the same CIC scheme applied to compute the grid density

void force(vector2d<double> phi,
		vector2d<double> &fx,
		vector2d<double> &fy);
// Computes force from potential
#endif
