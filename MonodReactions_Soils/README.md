This folder contains the following Matlab scrips and files:
#

- 'main_2D_Monod_Richards_L.m' is the main program for simulations of biodegradation in heterogeneous soils
	with transition from unsaturated to saturated flow regime modeled by a degenerate Richards equation,
	solved with an iterative GRW L-scheme. The transport step is solved with the iterative BGRW L-scheme.

- 'BGRW_2D_Monod_Richards_L.m' is the function for the nonreactive advection-diffusion step in the BGRW L-scheme.

- 'V_Kraichnan_Gauss_param.m' and 'K_r.m' are functions used to compute realizations of the saturated bydraulic
	conductivity modeled as a lognormal random function.

- 'initstate.m' is the function which initializes the random number generators used in the Kraichnan procedure.

- 'velocity.m' is a function to compute velocity components according to Darcy's law.

- 'theta_GM.m' is the function which provides the variable water content for saturated/unsaturated flow regimes
	according the the van Genuchten-Mualem model.

- 'reaction_Monod3.m' is the function for double Monod reactions with two mobile and one immobile species.

- 'plot_conv_Monod' is the function which plots norms of successive approximations used to assess the convergence 
	of the flow and transport L-schemes.

- 'convf', 'convc1', and 'convc2' are files containing convergence norms.

- 'subplot_surf_c.m' and 'subplot_surf_ptht.m' and are functions used to plot solutions at the final time T of the
	transport and flow equations, respectively. 

- 'subplot_solution.m' is the function which compares results for clay and loam soils.

- 'clayT1', 'clayT3', 'clayT5', 'loamT1', 'loamT3', and 'loamT5' are files containing solutions for clay and loam 
	soils at T=1,3, and 5 days.
