# 2D-Navier-Stokes-Solver

## Purpose

As the field of `Computational Fluid Dynamics (CFD)` progresses, 
the fluid flows are more and more analysed by using simulations with the help of high speed computers. 
In order to solve and analyse these fluid flows we require intensive simulation involving mathematical equations which governs the fluid flow, 
these are **Navier Stokes (NS) equation**. Solving these equations has become a necessity as almost every problem which is related to fluid flow analysis call for solving of Navier Stokes equation. 
These NS equations are partial differential equations so different numerical methods are used to solve these equations. 
Solving these partial differential equations so different numerical methods requires large amount of computing power and huge amount of memory is in play. 
Only practical feasible way to solve these equation is write a parallel program to solve them, which can then be run on powerful hardware capable of parallel processing to get the desired results High speed supercomputer will provide us very good performance in terms of reduction in execution time. 
In paper focus will be on finite volume as a numerical method. 
We will also see what GPGPU `(General-Purpose computing on Graphics Processing Units)` is and how we are taking its advantages to solve CFD problems.

## Fluids Theory

Flow past two and three dimensional bodies of different geometries has been under
intense research since such geometries are encountered in many of the real life
problems and engineering applications. In many areas of engineering, circular cylinders
form the basic component of structures and machinery. For example; nuclear
power plant cooling systems, offshore structures, heat exchangers, cooling towers,
submarines, transmission cables, pipelines and structural components of bridges
under the effect of wind, etc. Separated flow behind such bodies leads to the formation
of wake region which encompasses complex phenomena like recirculation,
shear layers and vortex shedding. Beyond a certain Reynolds number, the wake
becomes asymmetric giving rise to alternate shedding of vortices forming the well
known von-Kármán vortex street. The surface pressure changes each time a vortex is
shed fromthe cylinder. Therefore the cylinder is subjected to fluctuating drag and lift
in addition to themean force components. Such fluctuating forces may cause structural
vibrations. Under certain circumstances this can trigger structural failure or
enhance undesirable flow mixing in the wake. Therefore it is very important to have
a thorough understanding of the flow phenomena in order to appropriately control
vortex shedding in practical engineering environments.

Two-dimensional analysis of flowpast any arbitrary body geometry allows a broad
spectrum of parameters to be considered and provides a baseline for more detailed
study. Depending on the particular application, three-dimensional flow effects may
be significant, but the increased difficulty of the analysis would limit the scope of detailed
parametric study. In the present research work, an incompressible unsteady
two-dimensional finite volume Navier-Stokes solver is developed to investigate viscous
flow past two-dimensional geometries. In this solver, the full Navier-Stokes
equations have been solved numerically in the physical plane itself without using
any transformation to the computational plane. A **Consistent Flux Reconstruction
(CFR)** technique has been implemented on a collocated unstructured mesh comprising
of triangular cells. The present CFR unstructured grid solver has been named
as CFRUNS. This solver is applied to simulate unconfined flow past two side by side
cylinders by varying the center to center distance between the cylinders. The flow
problems were analysed by keeping the cylinders stationary, rotating and rotationally
oscillating about their respective axes. Streamlines, vorticity and ¸2 contours
are plotted to visualize instantaneous flow field. POD is used to extract the coherent
structures in the flow, Fast Fourier Transform (FFT) is used to obtain the frequency
spectrumof the flow and the time variation of the flow is studied through the associated
stress fields. Important parameters characterizing the flow such as lift and drag
coefficients and Strouhal number are also computed and quantitatively compared
with results of other researchers. The results obtained for single circular and square
cylinder are primarily used for validating the present CFRUNS solver. Reasonably
good comparison is obtained between the present results and results available from
literature.

## Working

### Mesh (Cells used for Computation)

![Mesh](https://github.com/singhsidhukuldeep/2D-Navier-Stokes-Solver/raw/master/assets/mesh.png)

The governing equations are solved in the physical plane within the defined computational
domain. The computational domain is discretized into a set of small triangular
cells. The field variables u, v, and p are defined at their cell centers. The
governing equations are discretized by finite volume technique based on the integral
form of the equations to be solved. For any arbitrary cell which is completely
immersed in the flow domain (­), the stencil used to solve the discretized equation
is given. The grid related information required to performthe flow calculations
are nodal co-ordinates, neighbouring cells and neighbour to neighbouring
cells.

## CREDITS

>Kuldeep Singh Sidhu

Github: [github/singhsidhukuldeep](https://github.com/singhsidhukuldeep)
`https://github.com/singhsidhukuldeep`

Website: [Kuldeep Singh Sidhu (Website)](http://kuldeepsinghsidhu.com)
`http://kuldeepsinghsidhu.com`

LinkedIn: [Kuldeep Singh Sidhu (LinkedIn)](https://www.linkedin.com/in/singhsidhukuldeep/)
`https://www.linkedin.com/in/singhsidhukuldeep/`
