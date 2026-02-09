# multiPhysicsSWInterFoam
This repository provides a set of custom [OpenFOAM®](https://www.openfoam.com/) solvers, the open-source CFD framework developed by OpenCFD Ltd., for coastal and ocean engineering applications.

It includes:
- **SWRiemannFoam**, a two-dimensional shallow water equation solver based on an approximate Riemann solver;
- **overInterDyMFoam+**, an extended three-dimensional Navier–Stokes solver with enhanced dynamic mesh and fluid–structure interaction capabilities;
- **multiPhysicsSWInterFoam**, a coupled multi-physics solver that integrates SWRiemannFoam and overInterDyMFoam+ to enable hybrid 2D–3D simulations.

Together, these solvers allow efficient modeling of large-scale wave propagation and localized, high-resolution flow interactions within a unified OpenFOAM framework.

<img src="docs/images/FOWTDinghy.png" width="100%">

# solvers
## SWRiemannFoam
SWRiemannFoam is a solver for the two-dimensional shallow water equations.

The solver is designed for simulating tsunami propagation and flood inundation problems,
and is intended for computing long-period waves such as tsunamis.

This solver was implemented based on the methodology presented in the following paper:

Murillo, J. and García-Navarro, P.,  
“Augmented versions of the HLL and HLLC Riemann solvers including source terms in one and two dimensions for shallow flow applications,”  
*Journal of Computational Physics*, Vol. 231, No. 20, pp. 6861–6906, 2012.  
https://www.sciencedirect.com/science/article/pii/S0021999112003464

<p align="center">
  <img src="docs/images/SWRiemann.png" width="35%">
  <img src="docs/images/tsunamiNankaiTrough.png" width="35%">
</p>

## overInterDyMFoam+
overInterDyMFoam+ is an extended three-dimensional Navier–Stokes solver for offshore and coastal flow simulations,
featuring enhanced dynamic mesh handling and fluid–structure interaction capabilities.

The solver is coupled with external physics engines to model rigid-body motion and mooring dynamics,
and is extended to account for internal fluid motion within floating or moving structures.

The integrated models include:
- **Bullet Physics** for rigid-body dynamics and collision handling  
  https://github.com/bulletphysics/bullet3
- **MoorDyn** for mooring system dynamics  
  https://github.com/FloatingArrayDesign/MoorDyn

<p align="center">
  <img src="docs/images/FOWT.png" width="65%">
</p>

## multiPhysicsSWInterFoam
multiPhysicsSWInterFoam is a solver for a wide range of coastal and ocean engineering problems.

The solver enables flexible coupling between a two-dimensional shallow water equation (SWE) solver
and a three-dimensional Navier–Stokes (NS) solver, allowing the simultaneous simulation of
large-scale, long-period wave propagation and localized high-resolution flow interactions.

The solver integrates the following numerical components:

- **Two-dimensional Shallow Water Equations (SWRiemannFoam)**
   - Computes large-scale wave propagation over the entire computational domain

- **Two-dimensional Laplace Equation**
   - Solves pressure distributions in narrow micro-scale gaps

- **overInterDyMFoam+**
   - Performs fluid–structure interaction (FSI) simulations using the overset mesh method
   - Models jointed and moored structures, including internal fluid motion

- **olaFlow (Higuera et al., 2017)**
   - Provides advanced wave generation capabilities
   - Simulates wave interaction with porous structures such as breakwaters
   - DOI: [![DOI](https://zenodo.org/badge/114743636.svg)](https://zenodo.org/badge/latestdoi/114743636)
   - Repository: https://github.com/phicau/olaFlow

<p align="center">
  <img src="docs/images/waveFlumes.png" width="100%">
</p>

The two-dimensional and three-dimensional domains can be coupled
via fixed connections at boundary or internal interfaces,
as well as through dynamically moving interfaces,
which allow the domains to be connected at arbitrary locations.

<p align="center">
  <img src="docs/images/arbitraryConnection.png" width="75%">
</p>

### note
In the current version of this solver, some memory deallocation procedures at the end of execution are not properly handled,
which may result in errors such as Bus error or Segmentation fault upon termination.

However, all output data and log files are generated correctly, and the results can be used without any problems.
This issue is planned to be addressed in a future update.

# compilation
To compile the solvers, first load the desired OpenFOAM environment,
then run the following script:

    ./Allwmake

Compilation has been tested and verified with OpenFOAM-v2512.
