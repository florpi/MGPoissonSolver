#ifndef particle_H
#define particle_H

#include <list>
#include <vector>

class Particle{
  private:
    double mass;
    double pos[2];
  public:
    Particle(double mass, double pos_x, double pos_y);
    double position(const size_t & i);
};
#endif
