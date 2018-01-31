#ifndef main_H 
#define main_H 

#include "grid.h"
#include <random>
#include <chrono>
#include <fstream>
#include <tclap/CmdLine.h>

// Default grid parameters
const int NGRID= 256; // Number of cells
const double L = 1.; // Simulation box length (periodic boundary conditions will be applied)
const int n_grids = 6; // Number of grids with fewer cells for the multigrid method
double h = L/NGRID; // Cell length
int n_iters_per_grid = 100; // Number of gauss seidel steps in each grid

// Default particles
int n_particles = 1; // Number of particles in grid
int n_test = 500;  // Number of test particles to compute force
double posx, posy; // position of source particle


bool paral_flag; // False for sequential gauss-seidel,
				// True for parallel version (blackred ordering)

#endif
