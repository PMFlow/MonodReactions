This folder contains the following Matlab scrips and files:
#

- 'main_1D_ReactTransp_test.m' is the main program, based on a non-iterative scheme, which is used to test the 
	convergence of the BGRW solver through comparisons with an analytical solution of advective-diffusive 
	transport with nonlinear bimolecular reactions in case of constant flow velocity and \theta=1.

- 'BGRW_1D.m' is the function for the nonreactive advection-diffusion step computed with the BGRW scheme.

- 'reaction.m' is the function for nonlinear bimolecular reactions.

- 'plot_Figs_React.m' is a function used to plot comparisons of numerical and analytical solutions.

- 'main_1D_ReactTransp_test_L.m' is the main porgram used to test the convergence by an iterative L-scheme.

- 'BGRW_1D.m' is the BGRW function for the nonreactive advection-diffusion step in the iterative L-scheme.

- 'plot_conv_React.m' is a function which plots the norm of the convergence criterium of the L-scheme
	versus the number of iterations at successive time points.

- 'convc1' and 'convc2' are files containing convergence norms.
