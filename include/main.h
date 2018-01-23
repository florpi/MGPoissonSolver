#ifndef main_H 
#define main_H 

#include "grid.h"
#include <random>
#include <chrono>
#include <fstream>

// Default grid parameters
const int NGRID= 256*20; // Number of cells
const double L = 1.; // Simulation box length (periodic boundary conditions will be applied)
double h = L/NGRID; // Cell length

int n_particles; // Number of particles in grid
double posx, posy; // position of particle


#endif
