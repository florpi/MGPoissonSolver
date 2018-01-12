#include "multigrid.h"
using namespace std;

#define TOL 0.0003 

Multigrid::Multigrid( const int ml, double h, const int mg):

	MAXGRID(mg),
	h(h),
	gridlevel(0),
	currentstep(1),
	n(0), // number of iterations carried
	maxlevel(ml),
	left(ml,mg,mg),
	right(ml,mg,mg),
	residual(ml,mg,mg)
	{}	


void Multigrid::Initial_conditions(vector2d<double> r )
{
		for(int i=0; i<MAXGRID; ++i)
			for(int j=0; j<MAXGRID;++j){
				this->right(0,i,j) = r(i,j);
			}
			
}

vector<int> applyBC(int i, int j, int step, int MAXGRID)
{
	int im, ip, jm, jp;
	im = i - step;
	ip = i + step;
	jm = j - step;
	jp = j + step;
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
void Multigrid::restrict(vector3d<double> &aux)
{
	int coarsestep = 2*currentstep;
	vector<int> boundary;
	for(int j=0; j<MAXGRID; j+=coarsestep)
		for(int i=0; i<MAXGRID; i+=coarsestep)
		{
		boundary = applyBC(i, j, currentstep,MAXGRID);
			//    boundary[0] = im boundary[1]= ip boundary[2]=jm boundary[3]=jp
			aux(gridlevel,i,j) = 0.25*aux(gridlevel,i,j) 
				+ 1./8. * (aux(gridlevel,i,boundary[3]) + aux(gridlevel,i,boundary[2]) + aux(gridlevel,boundary[1],j) + aux(gridlevel,boundary[0],j))
				   + (aux(gridlevel,boundary[1],boundary[3])+aux(gridlevel,boundary[1],boundary[2]) + aux(gridlevel,boundary[0],boundary[3])
				+ aux(gridlevel,boundary[0],boundary[2]))/16.; 
		}
	currentstep = coarsestep;
	gridlevel += 1;
	h = 2*h;

}
void Multigrid::prolong(vector3d<double> &aux)
{
	int finestep = currentstep/2;
	vector<int> boundary;
	vector<int> fineboundary;
	
	for(int j=finestep; j<MAXGRID; j+=currentstep)
		for(int i=0; i<MAXGRID; i+=currentstep)
		{
			fineboundary = applyBC(i,j,finestep,MAXGRID);
			aux(gridlevel,i,j) = 0.5*(aux(gridlevel,i,j+finestep) + aux(gridlevel,i,j-finestep));
		}
	for(int i = finestep; i<MAXGRID; i+= currentstep)
		for(int j = 0; j<MAXGRID; j+= currentstep)
		{
			fineboundary = applyBC(i,j,finestep,MAXGRID);
			aux(gridlevel,i,j) = 0.5*(aux(gridlevel,fineboundary[1],j) + aux(gridlevel,fineboundary[0],j));
		}	
	for(int i = finestep; i<MAXGRID; i+= currentstep)
		for(int j = finestep; j<MAXGRID; j+= currentstep)
		{
			fineboundary = applyBC(i,j,finestep,MAXGRID); 
			aux(gridlevel,i,j) = 0.25*(aux(gridlevel,fineboundary[0],fineboundary[2]) + aux(gridlevel,fineboundary[1],fineboundary[2]) + aux(gridlevel,fineboundary[0],fineboundary[3]) + aux(gridlevel,fineboundary[1],fineboundary[3]));
		}

currentstep = finestep;
gridlevel -=1;
h = h/2.;

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
		for(int i=0; i<MAXGRID; i+=currentstep){
				for(int j=0; j<MAXGRID; j+=currentstep){
					if((i+j)%2==0){	

						boundary = applyBC(i,j, currentstep,MAXGRID);		
						leftold = left(gridlevel,i,j);
						left(gridlevel,i,j) = 0.25*(left(gridlevel,boundary[0],j) + left(gridlevel,boundary[1],j)
								+ left(gridlevel,i,boundary[2]) + left(gridlevel,i,boundary[3]) 
								- h*h*right(gridlevel,i,j));
						change = fabs(left(gridlevel,i,j)/leftold -1.);
						if (change > maxChange) maxChange = change;
						if (maxChange < TOL && k>10 ) {
							cout << "Black converged after " << k << " iterations." << "\n" ;
							cancel=true;
						}
					}
				}
			}
	
//		cout << "Black Gauss-Seidel did not coverged in " << n_iters << "  iterations. The maximum difference is = " << maxChange << "\n";
		#pragma omp parallel for schedule(static) num_threads(4)
		for(int i=0; i<MAXGRID; i+=currentstep){
				for(int j=0; j<MAXGRID; j+=currentstep){
					if((i+j)%2==1){	
						boundary = applyBC(i,j, currentstep,MAXGRID);		
						leftold = left(gridlevel,i,j);
						left(gridlevel,i,j) = 0.25*(left(gridlevel,boundary[0],j) + left(gridlevel,boundary[1],j)
								+ left(gridlevel,i,boundary[2]) + left(gridlevel,i,boundary[3]) 
								- h*h*right(gridlevel,i,j));
						change = fabs(left(gridlevel,i,j)/leftold -1.);
						if (change > maxChange) maxChange = change;
						if (maxChange < TOL && k>10 ) {
							cout << "Red converged after " << k << " iterations." << "\n" ;
							cancel=true;
						}
					}
				}
			}
	}
	cout << "Gauss-Seidel did not coverged in " << n_iters << "  iterations. The maximum difference is = " << maxChange << "\n";

}

void Multigrid::compute_residual(){
	vector<int> boundary;
	vector3d<double> ddleft(maxlevel,MAXGRID,MAXGRID);
	for(int i=0; i<MAXGRID; i+= currentstep)
		for(int j=0; j<MAXGRID; j+= currentstep){
			boundary = applyBC(i,j,currentstep, MAXGRID);
			ddleft(gridlevel,i,j) = 1./h/h * (left(gridlevel,boundary[0],j) + left(gridlevel,boundary[1],j) 
					+ left(gridlevel,i,boundary[2])  
					+ left (gridlevel,i,boundary[3]) 
					- 4.*left(gridlevel,i,j));
			residual(gridlevel,i,j) = right(gridlevel,i,j)
				- ddleft(gridlevel,i,j);
	}
}
			

void Multigrid::result(int n_iters0){

	gauss(10); // First initial guess for level 0 (finest grid(
	
	vcycle(n_iters0);
}
//vector<vector<vector<double>>> Multigrid::get_left()
//{
	//return this->left;
//}

void Multigrid::vcycle(int n_iters){
	if(maxlevel -1  > n){

		compute_residual();
		restrict(residual); // +1 gridlevel

		cout << ">>>>> GRIDLEVEL  = " << gridlevel << "\n"; 
		for(int i=0; i<MAXGRID;i+= currentstep)
		   for(int j=0; j<MAXGRID; j+= currentstep){
				right(gridlevel,i,j) = residual(gridlevel-1,i,j); // Right hand side of current level is the restricted residual of previous level
		}

		gauss(n_iters); // Improves error on current level

		n+=1;
		vcycle(n_iters);
	}
	else if(this->gridlevel>0)
	{
		
//		print(left, MAXGRID, MAXGRID, gridlevel);
		prolong(left); // -1 gridlevel
		cout << ">>>>> GRIDLEVEL  = " << gridlevel << "\n"; 
		
		for(int i=0; i<MAXGRID;i+= currentstep)
		   for(int j=0; j<MAXGRID; j+= currentstep)
		   {

				left(gridlevel,i,j) += left(gridlevel+1,i,j); 	   

			}
		gauss(n_iters);
	n += 1;
		vcycle(n_iters);
	}

}
