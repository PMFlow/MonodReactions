This folder contains the following Matlab scrips and files:
#

- 'main_2D_Monod_Richards_test_Sat.m' is the main program used to test the convergence of the GRW flow L-scheme 
	coupled with a NONITERATIVE BGRW scheme through comparisons with an analytical solution of flow and 
	reactive transport with nonlinear Monod reactions in case of SATURATED flow regime.

- 'F_saturated.m' is a function providing the source terms derived from analytical solutions and the coupled 
	system of equations for saturated flow coupled with the Monod-reactive transport model.

- 'velocity.m' is a function to compute velocity components according to Darcy's law.

- 'theta.m' is the function which provides the water content (constant in the present case of saturated flow).

- 'plot_conv.m' is a function to plot norms of convergence criterium for the flow GRW L-scheme.

- 'convf' is the file containing convergence norms.

- 'BGRW_2D_Monod_Richards.m' is the BGRW function for the nonreactive advection-diffusion step.

- 'reaction_Monod2.m' is the function for double Monod reactions with constant biomass concentration c_3=1.
