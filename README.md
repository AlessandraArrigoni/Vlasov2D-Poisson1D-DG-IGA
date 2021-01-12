# Vlasov2D-Poisson1D-DG-IGA
Repository with the MATLAB code used for my [master thesis](https://www.politesi.polimi.it/bitstream/10589/166380/1/2020_10_Arrigoni.pdf) *"Comparison of Eulerian and semi-Lagrangian discontinuous Galerkin methods for Vlasov-Poisson system"*, defended in October 2020 at EPFL (Lausanne, CH) and PoliMI (Milan, IT).

### Dependencies and references
The code needs the MATLAB/Octave library [GeoPDEs](http://rafavzqz.github.io/geopdes/) to solve the 1D Poisson problem with a classic isogeometric finite element method. In the code we refer to the following three research papers containing the description of the methods and the parameters for the tests:
* B. Ayuso De Dios and S. Hajian, 2012, *High order and energy preserving discontinuous Galerkin methods for the Vlasov-Poisson system*, (https://arxiv.org/abs/1209.4025)
* N. Crouseilles et al., 2011, *Discontinuous Galerkin semi-Lagrangian method for Vlasov-Poisson*, (https://doi.org/10.1051/proc/2011022)
* J. A. Rossmanith and D. C. Seal, 2011, *A positivity-preserving high-order semi-Lagrangian discontinuous Galerkin scheme for the Vlasov-Poisson equations*, (https://doi.org/10.1016/j.jcp.2011.04.018)

### Folders content
* `DG_QUAD_AD` : functions to build 2D cartesian grids and discontinuous FE spaces in one and two dimensions; functions to assemble the matrices for the 1D and 2D advection problem; functions to compute the errors and to visualize the DG numerical solution (create `.vtk` files).
* `EULERIAN-DG\timeschemes` : functions for the three temporal schemes used in the Eulerian approach (4th order Runge-Kutta and two different II order splitting schemes with Crank-Nicolson method).
* `Helpers` : functions based on GeoPDEs library for the solution of Poisson problem and the coupling of IGA and DG methods.
* `NeighboursStructures` : `.mat` files storing precomputed variables linked to the mesh structures.
* `SEMILG-DG` : function with the full semi-Lagrangian DG method for the Vlasov-Poisson system and basic tests to check the implementation.
* `Tests` : scripts with parameters for all the tests performed in the master thesis (convergence, weak Landau damping and two streams instability) with all the considered methods.
