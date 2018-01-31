#include "multigrid.h"
using namespace std;

#define TOL 0.0001

Multigrid::Multigrid(const int ml):

	gridlevel(0),
	n(0), // number of iterations carried
	maxlevel(ml)
	{}	



void Multigrid::Initial_conditions(Grid mother){
	/* Sets initial grid and creates bigger cell ones*/
		grids.push_back(mother);
		for(int k=1; k<maxlevel; ++k){
			int new_grid = mother.MAXGRID / (2*k);
			Grid gr(new_grid, 1.);
			grids.push_back(gr);
		}
		
}

vector<int> applyBC(int i, int j, int MAXGRID){
	/* Applies periodic boundary conditions */
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
void Multigrid::restrict(){
	/* Interpolates right hand side to coarser grid */
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

void Multigrid::prolong(){
	/* Interpolates error to thinner grid */ 
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
	int i,j;
	double change;
	double leftold;
	int iter = 0;
	double maxChange= 2*TOL;
	while (maxChange > TOL && iter < n_iters){
		maxChange = 0.;
		for(i=0; i<grids[gridlevel].MAXGRID; ++i){
			for(j = 0; j<grids[gridlevel].MAXGRID ; ++j){
					leftold = grids[gridlevel].lhs(i,j);
					int im = i - 1;
					int ip = i + 1;
					int jm = j - 1;
					int jp = j + 1;	
					if(im<0){
						im = grids[gridlevel].MAXGRID -1;
					}
					if(ip>grids[gridlevel].MAXGRID-1){
						ip = 0;
					}
					if(jm<0){
						jm = grids[gridlevel].MAXGRID -1;
					}
					if(jp>grids[gridlevel].MAXGRID - 1){
						jp = 0;
					}
					grids[gridlevel].lhs(i,j) = 0.25*(grids[gridlevel].lhs(im,j) + grids[gridlevel].lhs(ip,j)
							+ grids[gridlevel].lhs(i,jm)+ grids[gridlevel].lhs(i,jp)
							- grids[gridlevel].h*grids[gridlevel].h*grids[gridlevel].rhs(i,j));
					change = fabs(grids[gridlevel].lhs(i,j)/leftold -1.);

					if (change > maxChange) maxChange = change;

					} 
	}
	
	iter += 1;
}
/*
if(maxChange < TOL)
{
	cout << "Converged after  " << iter << " iterations." << endl;
}
else{
 	cout << "Gauss-Seidel did not coverged in " << n_iters << "  iterations. The maximum difference is = " << maxChange << "\n";
}
*/

}


void Multigrid::gauss_omp(int n_iters){
	/* Parallel gauss seidel with black red ordering*/
	int i,j;
	double maxChangeB, maxChangeR, maxChange;
	double change;
	double leftold;
	int iter = 0;
	double diff = 1.1*TOL;
	while (diff > TOL && iter < n_iters){
		diff = 0.;
		#pragma omp parallel private(j, maxChange, maxChangeR, maxChangeB)
		{
		maxChangeB = 0.;
		#pragma omp for schedule(static) 
			for(i=0; i<grids[gridlevel].MAXGRID; ++i){
				for(j = i%2; j<grids[gridlevel].MAXGRID - (i+1)%2; j+=2){
					leftold = grids[gridlevel].lhs(i,j);
					int im = i -2;
					int ip = i + 2;
					int jm = j - 2;
					int jp = j + 2;	
					if(im<0){
						im = grids[gridlevel].MAXGRID -1;
					}
					if(ip>grids[gridlevel].MAXGRID-1){
						ip = 0;
					}
					if(jm<0){
						jm = grids[gridlevel].MAXGRID -1;
					}
					if(jp>grids[gridlevel].MAXGRID-1){
						jp = 0;
					}
					grids[gridlevel].lhs(i,j) = 0.25*(grids[gridlevel].lhs(im,j) + grids[gridlevel].lhs(ip,j)
							+ grids[gridlevel].lhs(i,jm)+ grids[gridlevel].lhs(i,jp)
							- 4*grids[gridlevel].h*grids[gridlevel].h*grids[gridlevel].rhs(i,j));
					change = fabs(grids[gridlevel].lhs(i,j)/leftold -1.);

					if (change > maxChangeB) maxChangeB = change;
					} 
	}
		maxChangeR = 0.;
	#pragma omp for schedule(static) 

			for(int i=0; i<grids[gridlevel].MAXGRID; ++i){
				for(int j= (i+1)%2; j<grids[gridlevel].MAXGRID-i%2; j+=2){
					
					leftold = grids[gridlevel].lhs(i,j);
					int im = i -2;
					int ip = i + 2;
					int jm = j - 2;
					int jp = j + 2;	
					if(im<0){
						im = grids[gridlevel].MAXGRID -1;
					}
					if(ip>grids[gridlevel].MAXGRID-1){
						ip = 0;
					}
					if(jm<0){
						jm = grids[gridlevel].MAXGRID -1;
					}
					if(jp>grids[gridlevel].MAXGRID-1){
						jp = 0;
					}
					grids[gridlevel].lhs(i,j) = 0.25*(grids[gridlevel].lhs(im,j) + grids[gridlevel].lhs(ip,j)
							+ grids[gridlevel].lhs(i,jm)+ grids[gridlevel].lhs(i,jp)
							- grids[gridlevel].h*grids[gridlevel].h*grids[gridlevel].rhs(i,j));
	

					change = fabs(grids[gridlevel].lhs(i,j)/leftold -1.);
					if (change > maxChangeR) maxChangeR = change;
					}
	}
	maxChange = maxChangeB < maxChangeR ? maxChangeR : maxChangeB;
	#pragma omp barrier
	if(diff < maxChange)
	#pragma omp atomic write
		diff = maxChange;
	} // Ends parallelisation


	iter += 1;
}

//cout << "Gauss-Seidel did not coverged in " << n_iters << "  iterations. The maximum difference is = " << diff << "\n";
	
}

void Multigrid::compute_residual(){
	/* Computes difference between approximated result and real */
	vector<int> boundary;
	vector2d<double> ddleft(grids[gridlevel].MAXGRID,grids[gridlevel].MAXGRID);
	for(int i=0; i<grids[gridlevel].MAXGRID; ++i)
		for(int j=0; j<grids[gridlevel].MAXGRID;++j){
			boundary = applyBC(i,j, grids[gridlevel].MAXGRID);
			ddleft(i,j) = 1./grids[gridlevel].h/grids[gridlevel].h * (grids[gridlevel].lhs(boundary[0],j) + grids[gridlevel].lhs(boundary[1],j) 
					+ grids[gridlevel].lhs(i,boundary[2])  
					+ grids[gridlevel].lhs(i,boundary[3]) 
					- 4.*grids[gridlevel].lhs(i,j));
			grids[gridlevel].residual(i,j) = grids[gridlevel].rhs(i,j)
				- ddleft(i,j);
	}
}
			


void Multigrid::vcycle(int n_iters, bool paral){
	/* Solves equation in coarser grids and uses solutions to accelerate convergence in thinner grids*/
	gauss(10); // First initial guess for level 0 (finest grid)

	if(maxlevel -1  > n){

		compute_residual();

		restrict(); // +1 gridlevel
		//cout << ">>>>> GRIDLEVEL  " << gridlevel << "  >>>>>>"  << "\n"; 
		if(paral){
			gauss_omp(n_iters);
		}
		else{
			gauss(n_iters);
		}
		n+=1;
		vcycle(n_iters,paral);
	}
	else if(this->gridlevel>0)
	{
		
		prolong(); // -1 gridlevel
		//cout << ">>>>> GRIDLEVEL  " << gridlevel << "  >>>>>>>" << "\n"; 
		for(int i=0; i<grids[gridlevel].MAXGRID; ++i)
		   for(int j=0; j<grids[gridlevel].MAXGRID; ++j )
		   {

				grids[gridlevel].lhs(i,j) += grids[gridlevel].error(i,j); 	   

			}
		if(paral){
			gauss_omp(n_iters);
		}
		else{
			gauss(n_iters);
		}
	
	n += 1;
		vcycle(n_iters, paral);
	}

}
