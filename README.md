# MGPoissonSolver
An efficient 2d Multigrid gauss-seidel code to solve Poisson type equations. Gauss-seidel iterations are computationaly expensive, however using multiple grids with different resolutions helps convergence and brings execution times to reasonable levels. One of the main advantages to solve the gravitational force equation with a multigrid method is that it can be extended to non-linear equations.

For a point particle we obtain the following acceleration-distance relation in two dimensions:
![Alt text](results/results.png?raw=true "acceleration")
