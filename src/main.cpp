#include "main.h"
using namespace std;

// Random seed for particle's positions
default_random_engine dre(chrono::steady_clock::now().time_since_epoch().count());
double random(double lim)
{
	uniform_real_distribution<> uid{0.,lim};
	return uid(dre);
}


int main(int argc, char* argv[]){

	ofstream myfile;
	if (argc ==1){
		paral_flag = 0; // Sequential code activated by default
	}
	else{
		paral_flag = atoi(argv[1]); // 0 means sequential, 1 parallel version
	}
	// Define time measures:
	clock_t t_initial, t_final;
	float seconds;
	t_initial=clock();

	// Define grid :
   	Grid mother(NGRID,L);

	// Define initial particle(s) parameters
	// and add them to the grid (in this case only one
	// since we want to check the force of an individual particle, let's add only one source particle to the grid
	posx = random(1.);
	posy = random(1.);;
	Particle part(1.,posx,posy);
	mother.add_particle(part);
	mother.compute_density();
	// Generate multiple grids 
	Multigrid mg(n_grids);
	mg.Initial_conditions(mother);
	mg.vcycle(n_iters_per_grid, paral_flag);
	for( int i=0; i<NGRID;++i)
		for(int j=0; j<NGRID;++j){
			mother.lhs(i,j) = mg.grids[0].lhs(i,j);
	}
	// Compute force in initial grid
	mother.compute_force();
	
	// Generate test particles to compute the acceleration they suffer
	// due to grid particles
	double acc, rmin, rmax, p, angle, r,alpha,xs,ys,dx,dy;
	// Save them in file
	myfile.open("../results/acceleration" + to_string(NGRID) + ".txt");
	for(int i=0;i<n_test;++i){
		// Randomly generate test particles at distances between 0.3*h and 0.5
		// from grid particle
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
		test_part.compute_acceleration(mother);
		myfile << r << "    " << test_part.get_acceleration()<< "\n" ;
	}

	myfile.close();
	t_final = clock();
	float diff = ((float)t_final - (float) t_initial);
	seconds = diff/CLOCKS_PER_SEC;
	cout << "The program took " << seconds << " seconds to finish" << "\n";
}
