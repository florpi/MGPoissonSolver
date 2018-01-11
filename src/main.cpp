#include "main.h"
using namespace std;
// TODO: MIRAR BIEN COMO HACER DEFINICIONES INICIALES, MAX DEPTH DEL MULTIGIRD

// Random seed for particle's positions
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
	int n_particles = 1;
	const int NGRID = 256;
	double L = 1.;
	double h = L/NGRID;
   	Grid gr(NGRID,10,1);
	double posx = 0.27;
	double posy = 0.25;
	for(int p=0; p<n_particles ; ++p){
		Particle part(1., posx,posy);
		gr.add_particle(part);
	}
	gr.compute_density();
	gr.compute_force();

	int n_test = 1000;
	double acc, rmin, rmax, p, angle, r,alpha,xs,ys,dx,dy;
	myfile.open("../results/acceleration" + to_string(NGRID) + ".txt");
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
		Particle test_part(1.,xs,ys);
		test_part.compute_acceleration(NGRID,h,gr.get_fx(), gr.get_fy());
		myfile << r << "    " << test_part.get_acceleration()<< "\n" ;
	}

	myfile.close();
	t_final = clock();
	float diff = ((float)t_final - (float) t_initial);
	seconds = diff/CLOCKS_PER_SEC;
	cout << "The program took " << seconds << " seconds to finish" << endl;
}
