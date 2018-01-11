#ifndef particle_H
#define particle_H

#include <list>
#include <vector>
#include "vector2d.cpp"

class Particle{
  private:
    double mass;
    double pos[2];
	double acc; // acceleration created by the other particles
  public:
    Particle(double mass, double pos_x, double pos_y);
    double position(const size_t & i);
	void compute_acceleration(int MAXGRID, double h,vector2d<double> fx, vector2d<double> fy);
	double get_acceleration();
};
#endif
