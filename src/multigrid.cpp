#include "multigrid.h"
using namespace std;

#define TOL 0.0003 

Multigrid::Multigrid( const int ml, const int mg,const int n_particles):

	MAXGRID(mg),
	h(1./mg),
	gridlevel(0),
	currentstep(1),
	n(0), // number of iterations carried
	n_particles(n_particles),
	maxlevel(ml)
	{}	



void Multigrid::Initial_conditions(double posx, double posy)
{
		Grid gr(MAXGRID, 1.);
		grids.push_back(gr);
		for(int p=0; p<n_particles; ++p){
			Particle part(1.,posx,posy);
			grids[0].add_particle(part);
		}
		grids[0].compute_density();
		cout << "esto va" << endl;
		for(int k=1; k<maxlevel; ++k){
			int new_grid = MAXGRID / (2*k);
			Grid gr(new_grid, 1.);
			grids.push_back(gr);
		}
		
}

vector<int> applyBC(int i, int j, int MAXGRID)
{
	int im, ip, jm, jp;
	im = i - 1;
	ip = i + 1;
	jm = j - 1;
	jp = j + 1;
	if(im<0){
		im += MAXGRID;
	}
	if(ip>MAXGRID-1){
		ip -= MAXGRID;
	   } 
	if(jm<0){
		jm += MAXGRID;
	}
	if(jp>MAXGRID-1){
		jp -= MAXGRID;
	}
	vector<int> boundary{im, ip, jm, jp};
	return boundary;

}
int Multigrid::get_currentstep()
{
	return this->currentstep;
}
void Multigrid::restrict()
{
	vector<int> boundary;
	for(int j=0; j<grids[gridlevel+1].MAXGRID; ++j)
		for(int i=0; i<grids[gridlevel+1].MAXGRID; ++i)
		{
		boundary = applyBC(2*i,2*j,grids[gridlevel].MAXGRID);
			grids[gridlevel+1].rhs(i,j) = 0.25*grids[gridlevel].residual(2*i,2*j) 
				+ 1./8. * (grids[gridlevel].residual(2*i,boundary[3]) + grids[gridlevel].residual(2*i,boundary[2]) + grids[gridlevel].residual(boundary[1],2*j) + grids[gridlevel].residual(boundary[0],2*j))
				   + (grids[gridlevel].residual(boundary[1],boundary[3])+grids[gridlevel].residual(boundary[1],boundary[2]) + grids[gridlevel].residual(boundary[0],boundary[3])
				+ grids[gridlevel].residual(boundary[0],boundary[2]))/16.; 
		}
	gridlevel += 1;
}

void Multigrid::prolong()
{
	vector<int> boundary;
	for(int j=0; j < grids[gridlevel-1].MAXGRID; ++j)
		for(int i=0; i < grids[gridlevel-1].MAXGRID; ++i){
		boundary = applyBC(i/2,j/2,grids[gridlevel].MAXGRID);
			if(i%2 == 0 && j%2 == 0){
				grids[gridlevel-1].error(i,j) = grids[gridlevel].lhs(i/2,j/2); // Direct injection
			}
			else if( i%2 == 0 && j%2 == 1){
				grids[gridlevel-1].error(i,j) = 0.5*(grids[gridlevel].lhs(i/2, j/2) + grids[gridlevel].lhs(i/2,boundary[3]));

			}
			else if ( i%2 == 1 && j%2==0){
				grids[gridlevel-1].error(i,j) = 0.5*(grids[gridlevel].lhs(i/2, j/2) + grids[gridlevel].lhs(boundary[1], j/2));
			}
			else if( i%2 == 1 && j%2 == 1){
				int ip = (i+1)/2;
			  	if(ip >= grids[gridlevel].MAXGRID){
					ip = 0;
				}
				int jp = (j+1)/2;
				if(jp >= grids[gridlevel].MAXGRID){
					jp = 0;
				}
				grids[gridlevel-1].error(i,j) = 0.25*(grids[gridlevel].lhs((i-1)/2 , (j-1)/2) + grids[gridlevel].lhs(ip, (j-1)/2)
					   							+ grids[gridlevel].lhs((i-1)/2, jp) + grids[gridlevel].lhs(ip,jp));
		
			}
	}	


gridlevel -=1;
}

void Multigrid::gauss(int n_iters)
{
	vector<int> boundary;
	double maxChange;
	double change;
	double leftold;
	omp_set_dynamic(0);
	for(int k=0;k<n_iters;k++){
	maxChange = 0.0;
	bool cancel = false;
	#pragma omp parallel for schedule(static) num_threads(4)
			for(int i=0; i<grids[gridlevel].MAXGRID; ++i){
				for(int j=0; j<grids[gridlevel].MAXGRID; ++j){
					
					boundary = applyBC(i,j, grids[gridlevel].MAXGRID);		
					leftold = grids[gridlevel].lhs(i,j);
					grids[gridlevel].lhs(i,j) = 0.25*(grids[gridlevel].lhs(boundary[0],j) + grids[gridlevel].lhs(boundary[1],j)
							+ grids[gridlevel].lhs(i,boundary[2])+ grids[gridlevel].lhs(i,boundary[3])
							- h*h*grids[gridlevel].rhs(i,j));
					change = fabs(grids[gridlevel].lhs(i,j)/leftold -1.);
					if (change > maxChange) maxChange = change;
					if (maxChange < TOL && k>10 ) {
						cout << "Converged after " << k << " iterations." << "\n" ;
						return;
					}
					} }
	}
	cout << "Gauss-Seidel did not coverged in " << n_iters << "  iterations. The maximum difference is = " << maxChange << "\n";
}

void Multigrid::compute_residual(){
	vector<int> boundary;
	vector2d<double> ddleft(grids[gridlevel].MAXGRID,grids[gridlevel].MAXGRID);
	for(int i=0; i<grids[gridlevel].MAXGRID; ++i)
		for(int j=0; j<grids[gridlevel].MAXGRID;++j){
			boundary = applyBC(i,j, grids[gridlevel].MAXGRID);
			ddleft(i,j) = 1./h/h * (grids[gridlevel].lhs(boundary[0],j) + grids[gridlevel].lhs(boundary[1],j) 
					+ grids[gridlevel].lhs(i,boundary[2])  
					+ grids[gridlevel].lhs(i,boundary[3]) 
					- 4.*grids[gridlevel].lhs(i,j));
			grids[gridlevel].residual(i,j) = grids[gridlevel].rhs(i,j)
				- ddleft(i,j);
	}
}
			

void Multigrid::result(int n_iters0){

	gauss(10); // First initial guess for level 0 (finest grid(
	
	vcycle(n_iters0);
}

void Multigrid::vcycle(int n_iters){
	if(maxlevel -1  > n){

		compute_residual();
		restrict(); // +1 gridlevel
		cout << ">>>>> GRIDLEVEL  = " << gridlevel << "\n"; 
		gauss(n_iters); // Improves error on current level
		n+=1;
		vcycle(n_iters);
	}
	else if(this->gridlevel>0)
	{
		
//		print(left, MAXGRID, MAXGRID, gridlevel);
		prolong(); // -1 gridlevel
		cout << ">>>>> GRIDLEVEL  = " << gridlevel << "\n"; 
		for(int i=0; i<grids[gridlevel].MAXGRID; ++i)
		   for(int j=0; j<grids[gridlevel].MAXGRID; ++j )
		   {

				grids[gridlevel].lhs(i,j) += grids[gridlevel].error(i,j); 	   

			}
		gauss(n_iters);
	n += 1;
		vcycle(n_iters);
	}

}
