## Data repository accompanying the paper "Global random walk solvers for reactive transport and biodegradation processes in heterogeneous porous media"

This is a data repository accompanying the paper Global random walk solvers for reactive transport and biodegradation processes in heterogeneous porous media" 
by Nicolae Suciu and Florin A. Radu

The repository contains Matlab codes based on Global Random Walk (GRW) algorithms to solve coupled nonlinear problems of flow and reactive transport in porous media. 
Solutions of the Richards equation are obtained with the GRW L-scheme introduced in a previous publication (https://doi.org/10.1016/j.advwatres.2021.103935). 
The coupled system of multicomponent transport equations with nonlinear reactions of Monod type is solved with new numerical schemes based on both unbiased GRW and biased (BGRW) 
algorithms. The nonlinearity of the reaction system is solved with explicit non-iterative schemes, in case of saturated flow regime, and with linearization L-schemes if 
the flow domain contains both saturated and nonsaturated regions.

Convergence tests and applications of the GRW/BGRW schemes are included in the following folders:

#
## CodeVerification
with subfolders
### BimolecularReaction_ConstantVelocity
### Monod_SaturatedFlows
### Monod_UnsaturatedFlows
### Monod_VariablySaturatedFlows

- contains code verification and convergence tests preformed with analytical manufactured solutions.
#
## MonodReactions_Aquifers
- contains Matlab codes and functions, as well as simulation results for coupled saturated-flow and multicomponent Monod reactions, obtained with GRW and BGRW explicit 
    linearization schemes.
#
## MonodReactions_Soils 
- contains codes, functions, and results illustrating the use of the BGRW L-scheme to solve a coupled system of flow and multicomponent transport equations, with
    biodegradation governed by a double Monod model, in the general case with transition from unsaturated to saturated flow regime. 
