#include "grid.h"
//TODO: particles.reserve(MAXNPART), donde?

Grid::Grid( int MAXGRID, const double L):
	MAXGRID(MAXGRID),
	L(L),
	h(L/MAXGRID), 
	G(1.),
	step(1),
	rhs(MAXGRID,MAXGRID),
	lhs(MAXGRID,MAXGRID),
	residual(MAXGRID,MAXGRID),
	error(MAXGRID,MAXGRID),
	fx(MAXGRID,MAXGRID),
	fy(MAXGRID,MAXGRID)
	{}

void Grid::add_particle(Particle part){
	/* Adds particle part to given grid */
	particles.push_back(part);
}

void Grid::compute_density(){
	/* Computes grid density using CIC mass assignment*/
	for(int p=0; p < particles.size(); ++p){	
		int i, ii, j, jj; 
		double xx = particles[p].position(0)/h;

		i = (int) xx;
		double yy = particles[p].position(1)/h;

		j = (int) yy; 
		double u = xx - i;
		double v = yy - j;


		ii = i + 1;
		if(ii >= MAXGRID)
		  ii = 0;
		jj = j + 1;
		if(jj >= MAXGRID)
		jj=0;

		this->rhs(i,j) += 4*M_PI*G*(1-u)*(1-v)/(h * h);
		this->rhs(ii,j) += 4*M_PI*G*(u)*(1-v)/(h*h);;
		this->rhs(i,jj) += 4*M_PI*G*v * (1-u)/(h*h); 
		this->rhs(ii,jj) += 4*M_PI*G*u*v/(h*h);
	}
}

double Grid::get_density(int i, int j){
	/* Obtains private variable density outside the class*/
	return this->rhs(i,j);
}

void Grid::compute_force(){
	/* Computes force on grid */
	vector<int> boundary;
	for( int i=0; i<MAXGRID;++i)
		for(int j=0; j<MAXGRID;++j)
		{
		 	boundary = applyBC(i,j,MAXGRID);	
			fx(i,j) = -(lhs(boundary[1],j) - lhs(boundary[0],j)) /(2.*h);
			fy(i,j) = - (lhs(i,boundary[3]) - lhs(i,boundary[2]))/(2.*h);
		}
}

vector2d<double> Grid::get_fx(){
	/* Obtains private variable fx outside the class*/
	return this->fx;
}
vector2d<double> Grid::get_fy(){
	/* Obtains private variable fy outside the class*/
	return this->fy;
}


