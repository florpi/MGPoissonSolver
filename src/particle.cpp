#include "particle.h"
#include <stdio.h>
#include <math.h>
#include "iostream"
#include  "grid.h"

Particle::Particle(double mass, double pos_x, double pos_y):
  mass(mass)
  {
    pos[0] = pos_x; pos[1] = pos_y;
  }

double Particle::position(const size_t & i){
  if (i!=0 && i!=1) exit(1);
  return pos[i];
}

void Particle::compute_acceleration(Grid &gr){
	int i,j;
	double xx = position(0)/gr.h;
	i = (int) xx;
	double yy = position(1)/gr.h;
	j = (int) yy;

	double u, v, ax, ay, a;
	u = xx - i;
	v = yy-j;

	int ii = i +1;
	if (ii>=gr.MAXGRID)
		ii = 0;
	int jj = j + 1;
	if (jj>=gr.MAXGRID)
		jj = 0;
	ax = gr.fx(i,j) * (1-u)*(1-v) + gr.fx(ii,j)*u*(1-v) +	gr.fx(i,jj)*v*(1-u) + gr.fx(ii,jj)*u*v;
	ay = gr.fy(i,j)*(1-v)*(1-u)+ gr.fy(ii,j)*u*(1-v) +
		gr.fy (i,jj)*v*(1-u) + gr.fy(ii,jj)*u*v;
	this->acc=sqrt(ax*ax + ay*ay);
}
double Particle::get_acceleration(){
	return this->acc;
}
