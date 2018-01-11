#include "main.h"
using namespace std;

default_random_engine dre(chrono::steady_clock::now().time_since_epoch().count());
double random(double lim)
{
	uniform_real_distribution<> uid{0.,lim};
	return uid(dre);
}

int main(){
	ofstream myfile;

	clock_t t_initial, t_final;
	float seconds;
	t_initial=clock();
	currentstep = 1;
	MAXGRID = 256;
	myfile.open("../results/acceleration"+ to_string(MAXGRID) +".txt");
	L = 1.;
	h = L/MAXGRID;
	maxlevel = 4;

	gridlevel = 0;
	minlevel =0;
	
	posx = 0.27;
	posy = 0.25;
	vector2d<double> phi(MAXGRID,MAXGRID);
	vector2d<double> fx(MAXGRID,MAXGRID);
	vector2d<double> fy(MAXGRID,MAXGRID);

	vector3d<double> left(maxlevel, MAXGRID, MAXGRID);
	vector3d<double> right(maxlevel, MAXGRID, MAXGRID);
	vector3d<double> residual(maxlevel, MAXGRID, MAXGRID);

	Particle * part = new Particle(1.,posx,posy);
	density(part,right); 
	
	Multigrid * mg = new Multigrid(gridlevel, currentstep, maxlevel, h, MAXGRID);
	mg->Initial_conditions(left,right,residual);
	mg->result(100);
	for(int i=0; i<MAXGRID; i++)
		for(int j=0; j<MAXGRID;j++){
			phi(i,j) = mg->left(0,i,j);
		}
	force(phi, fx, fy);
	int n_test = 1000;
	double acc, rmin, rmax, p, angle, r,alpha,xs,ys,dx,dy;
	for(int i=0;i<n_test;i++){
		rmin = 0.3*h;
		rmax = 0.5;
		p = random(1.);
		r = rmin*pow(rmax/rmin,p);
		angle= random(1.);
		dx = r*cos(2*M_PI*angle);
		dy = r*sin(2*M_PI*angle);
		xs = posx + dx;
		ys = posy + dy;
		while(xs < 0){
			xs += 1.;
		}
		while(xs > 1.){
			xs -= 1.;
		}
		while(ys<0.){
			ys += 1.;
		}
		while(ys>1.){
			ys -= 1.;
		}
		Particle * test_particle = new Particle(1.,xs,ys);
		acc = invert_force(test_particle,fx,fy);
		delete test_particle;
		myfile << r << "    " << acc << "\n" ;
	}

	myfile.close();
	t_final = clock();
	float diff = ((float)t_final - (float) t_initial);
	seconds = diff/CLOCKS_PER_SEC;
	cout << "The program took " << seconds << " seconds to finish" << endl;
	delete part;
	
	delete mg;

	return 0;
}
