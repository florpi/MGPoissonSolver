# MGPoissonSolver

An efficient 2d Multigrid gauss-seidel code to solve Poisson type equations. Gauss-seidel iterations are computationaly expensive, however using multiple grids with different resolutions helps convergence and brings execution times to reasonable levels. One of the main advantages to solve the gravitational force equation with a multigrid method is that it can be extended to non-linear equations.

Table of contents
=================

  * [Requirements](#Requirements)
  * [Build](#Build)
  * [Usage](#usage)
  * [Tests](#tests)
  
## Requirements
```shell
git (optional)
make
g++

## Build
```shell
$ git clone https://github.com/florpi/MGPoissonSolver.git 
$ cd MGPoissonSolver
$ make
```
## Usage
```console
Usage: bhg [command flags]

USAGE: 

   ./acceleration[-p <int>] 


Where: 
   -p <int>,  --parallel_flag <int>
	0 -> Executes sequential gauss seidel
	1 -> Executes parallel Black-red gauss seidel
 Default: 0
```
## Tests
For a point particle we obtain the following acceleration-distance relation in two dimensions:
![Alt text](results/results.png?raw=true "acceleration")

