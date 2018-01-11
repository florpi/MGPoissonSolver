#include "potential.h"

#define MAXITERS 50000
#define INITGRID 1
void print(vector2d<double> &matrix, int rows, int cols)
{
	for(int r=0; r<rows; r++)
	{
		for(int c=0; c<cols; c++)
		{
			cout<< matrix(r,c)<<" ";
		}
		cout<<endl;
	}
}
void density(Particle * part, vector3d<double> &rho){
  int i, ii, j, jj; 
  double xx = part->position(0)/h;

  i = (int) xx;
  double yy = part->position(1)/h;

  j = (int) yy; 
  double u = xx - i;
  double v = yy - j;


  ii = i + 1;
  if(ii >= MAXGRID)
	  ii = 0;
  jj = j + 1;
 if(jj >= MAXGRID)
	jj=0;

    rho(0,i,j) += 4*M_PI*(1-u)*(1-v)/(h * h);
    rho(0,ii,j) += 4*M_PI*(u)*(1-v)/(h*h);;
    rho(0,i,jj) += 4*M_PI*v * (1-u)/(h*h); 
    rho(0,ii,jj) = 4*M_PI*u*v/(h*h);
}

