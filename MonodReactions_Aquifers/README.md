This folder contains the following Matlab scrips and files:
#

- 'main_2D_Monod_Sat_Kraichnan_bgrw.m' is the main program for simulations of biodegradation in saturated flows
	through heterogeneous aquifers. The flow velocity is approximated by the Kraichnan procedure. The 
	transport step is computed with the non-iterative BGRW scheme.

- 'main_2D_Monod_Sat_Kraichnan_grw.m' is the main program for simulations of the same biodegradation process 
	by the non-iterative  unbiased GRW scheme.

- 'BGRW_2D_Monod_Richards.m' is the BGRW function for the nonreactive advection-diffusion step.

- 'GRW_2D_Monod_Richards.m' is the unbiased GRW function for the nonreactive advection-diffusion step.

- 'V_Kraichnan_Gauss_param.m' and 'V_Kraichnan_Gauss_func.m' are functions used to compute realizations of the
	Kraichnan velocities.

- 'initstate.m' is the function which initializes the random number generators used in the Kraichnan procedure.

- 'reaction_Monod3.m' is the function for double Monod reactions with two mobile and one immobile species.

- 'plot_fig_Monod_Sat.m' and 'subplot_surf.m' are functions used to plot results at the final time T.

- 'bgrwT30.mat', 'bgrwT300.mat', and 'grwT300.mat' are files containing results at the final time T.

- 'Comparison_Monod_Sat.m' is the function which compares the BGRW and GRW results at the final time T=300.

- 'main_2D_MoSaKr_plot.m' is the main program used to compare results provided by the GRW solver for different 
	variances and correlation lengths of the Kraichnan velocity field.
