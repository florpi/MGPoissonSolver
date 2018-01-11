//#ifndef computeforce_H
//#define computeforce_H

#include "multigrid.h"
#include "particle.h"
#include "vector2d.cpp"

/*      VARIABLES      */
extern int currentstep ; // Step size, initial = 1
extern int MAXGRID ; // Grid size
extern double L; // Box length
extern double h; // Grid spacing

/*      FUNCTIONS      */
void print(vector2d<double> &matrix, int rows,int cols); // Outputs matrix.

void density(Particle *part, vector3d<double> &rho); // Assigns density to grid using the CIC scheme.
