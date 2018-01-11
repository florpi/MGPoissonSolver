#include "particle.h"
#include <stdio.h>
#include "iostream"

Particle::Particle(double mass, double pos_x, double pos_y):
  mass(mass)
  {
    pos[0] = pos_x; pos[1] = pos_y;
  }

double Particle::position(const size_t & i){
  if (i!=0 && i!=1) exit(1);
  return pos[i];
}
