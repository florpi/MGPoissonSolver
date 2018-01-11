#include "particle.h"
#include <stdio.h>
#include <math.h>
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

void Particle::compute_acceleration(int MAXGRID, double h,vector2d<double> fx, vector2d<double> fy){
	int i,j;
	double xx = position(0)/h;
	i = (int) xx;
	double yy = position(1)/h;
	j = (int) yy;

	double u, v, ax, ay, a;
	u = xx - i;
	v = yy-j;

	int ii = i +1;
	if (ii>=MAXGRID)
		ii = 0;
	int jj = j + 1;
	if (jj>=MAXGRID)
		jj = 0;
	ax = fx(i,j) * (1-u)*(1-v) + fx(ii,j)*u*(1-v) +	fx(i,jj)*v*(1-u) + fx(ii,jj)*u*v;
	ay = fy(i,j)*(1-v)*(1-u)+ fy(ii,j)*u*(1-v) +
		fy (i,jj)*v*(1-u) + fy(ii,jj)*u*v;
	this->acc=sqrt(ax*ax + ay*ay);
}
double Particle::get_acceleration(){
	return this->acc;
}
