This folder contains the following Matlab scrips and files:
#

- 'main_2D_Monod_Richards_L_test_Sat.m' is the main program used to test the convergence of the GRW flow L-scheme 
	coupled with an ITERATIVE BGRW L-scheme through comparisons with analytical solutions of flow and 
	reactive transport with nonlinear Monod reactions in case of UNSATURATED flow regime.

- 'F_unsaturated.m' is a function providing the source terms obtained by inserting the analytical solutions and their derivatives in the coupled 
	system of equations.

- 'velocity.m' is a function to compute velocity components according to Darcy's law.

- 'theta.m' is the function which provides the variable water content for unsaturated flow regime.

- 'plot_conv_Monod.m' is a function to plot norms of convergence criterium for flow and transport L-schemes.

- 'convfunsat', 'convc1unsat', and 'convc2unsat' are files containing convergence norms.

- 'BGRW_2D_Monod_Richards_L.m' is the BGRW function for the nonreactive advection-diffusion step used in the 
	iterative BGRW L-scheme.

- 'reaction_Monod2_L.m' is the function for double Monod reactions with constant biomass concentration c_3=1 
	used in the iterative BGRW L-scheme.
