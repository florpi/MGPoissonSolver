#include "force.h"
using namespace std;
double invert_force(Particle *part, vector2d<double> &fx, vector2d<double> &fy){
	int i,j;
	double u, v, ax, ay, a;
	double xx = part->position(0)/h;
	i = (int) xx;
	double yy = part->position(1)/h;
	j = (int) yy;

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
	a =sqrt(ax*ax + ay*ay);
	return a;
}


void force(vector2d<double> phi,vector2d<double> &fx,  vector2d<double> &fy){
	vector<int> boundary;
	for( int i=0; i<MAXGRID;i+=currentstep)
		for(int j=0; j<MAXGRID;j+=currentstep)
		{
		 	boundary = applyBC(i,j,currentstep,MAXGRID);	
			fx(i,j) = -(phi(boundary[1],j) - phi(boundary[0],j)) /(2.*h);
			fy(i,j) = - (phi(i,boundary[3]) - phi(i,boundary[2]))/(2.*h);
		}
}



