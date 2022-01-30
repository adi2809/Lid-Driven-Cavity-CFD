# fortran solvers for lid driven cavity problem : 

we implemented the grid generation code in fortran ; the grid to be used is the cartesian non-uniform grid. we have incorporated the Rhie and Chow momentum interpolation scheme in the code to help avoid any kind of chequerboard oscillations in the pressure field (caused due to rapid switching between two parallel pressure field solutions) 
