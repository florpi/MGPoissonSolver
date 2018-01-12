#include "grid.h"
//TODO: particles.reserve(MAXNPART), donde?

Grid::Grid( const int MAXGRID, const double L):
	MAXGRID(MAXGRID),
	L(L),
	h(L/MAXGRID), 
	G(1.),
	rho(MAXGRID,MAXGRID),
	phi(MAXGRID,MAXGRID),
	fx(MAXGRID,MAXGRID),
	fy(MAXGRID,MAXGRID)
	{}

void Grid::add_particle(Particle part){
	particles.push_back(part);
}

void Grid::compute_density(){
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

		this->rho(i,j) += 4*M_PI*G*(1-u)*(1-v)/(h * h);
		this->rho(ii,j) += 4*M_PI*G*(u)*(1-v)/(h*h);;
		this->rho(i,jj) += 4*M_PI*G*v * (1-u)/(h*h); 
		this->rho(ii,jj) += 4*M_PI*G*u*v/(h*h);
	}
}

double Grid::get_density(int i, int j){
	return this->rho(i,j);
}

void Grid::compute_force(){
	Multigrid mg(4,h, MAXGRID);
	mg.Initial_conditions(rho);
	mg.result(100);
	for( int i=0; i<MAXGRID;++i)
		for(int j=0; j<MAXGRID;++j){
			phi(i,j) = mg.left(0,i,j);
	}
	vector<int> boundary;
	for( int i=0; i<MAXGRID;++i)
		for(int j=0; j<MAXGRID;++j)
		{
		 	boundary = applyBC(i,j,1,MAXGRID);	
			fx(i,j) = -(phi(boundary[1],j) - phi(boundary[0],j)) /(2.*h);
			fy(i,j) = - (phi(i,boundary[3]) - phi(i,boundary[2]))/(2.*h);
		}
}

vector2d<double> Grid::get_fx(){
	return this->fx;
}
vector2d<double> Grid::get_fy(){
	return this->fy;
}


