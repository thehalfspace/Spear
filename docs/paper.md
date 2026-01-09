---
title: 'SPEAR: A Julia package for 2D earthquake cycle simulations'
tags:
  - Julia
  - geophysics
  - seismology
  - earthquake cycle
  - numerical simulations
authors:
  - name: Prithvi Thakur
    orcid: 0000-0000-0000-0000
    equal-contrib: true
    affiliation: "1, 2" # (Multiple affiliations must be quoted)
  - name: Yihe Huang
    equal-contrib: true # (This is how you can denote equal contributions between multiple authors)
    affiliation: 2
  - name: Yoshihiro Kaneko
    equal-contrib: true # (This is how you can denote equal contributions between multiple authors)
    affiliation: 3
  - name: Peng Zhai
    equal-contrib: true # (This is how you can denote equal contributions between multiple authors)
    affiliation: 2
affiliations:
 - name: Center for Computation and Visualization, Brown University
   index: 1
 - name: Department of Earth and Environmental Sciences, Uniersity of Michigan
   index: 2
 - name: Department of Geophysics, Kyoto University
   index: 3
date: 09 January 2026
header-includes:
  - |
    \makeatletter
    \ifdefined\citeproc\else
      \newcommand{\citeproc}[2]{#2}
    \fi
    \makeatother
bibliography: refs.bib

---

# Summary
Simulating long-term earthquake sequences numerically is a topical problem in earth sciences. Such simulations encompass a wide range of spatial and temporal scales, with spatial scales ranging from millimeters along frictional contact to kilometers of fault-length, and temporal scales ranging from milliseconds during earthquake ruptures to decades during the aseismic period. Existing methods to simulate such multi-scale, nonlinear problem can be broadly classified into boundary-based methods [@Rice:1993; @Lapusta:2000; @Barbot:2019], volume-based methods [@Kaneko:2011; @Erickson:2014; @Thakur:2020], and hybrid methods [@Ma:2015; @Abdelmeguid:2019]. The boundary-based methods are very efficient in terms of time complexity and are suitable to study fault-friction and slip history in two and three dimensions. However, such methods cannot capture the heterogeneities in the bulk domain, and the effects of fault zone structures commonly observed in natural faults. Volume-based methods, like finite- or spectral-element are more suitable for such problems, at the cost of higher computation time.

This package contains a Julia implementation of spectral element method to solve two-dimensional fault-slip problem with linear elasticity, accounting for the dynamic inertial 
effects and bulk heterogeneities.


# Statement of need
Spear is a Julia package for 2D simulations of long-term fault-slip on strike-slip fault systems. Julia enables high-performance speeds while writing at a high-level programming interface, and is suitable for complex numerical simulations. Our package helps us perform complex numerical simulations with off-fault material heterogeneities and full inertial effects of seismic waves. Written in easy-to-use language Julia, and with minimal dependencies, our package is designed for users familiar with computational seismology, and users with minimal Julia or programming experience should be able to get started. 

The purpose of this package is to simulate long-term fault-slip dynamics with a focus on off-fault material heterogeneities and fault friction. We can study a wide variety of scientific problems like the effects of dynamic wave reflections due to the low-velocity fault zones, asymmetrical and complex fault zone structures, and time-dependent bulk property changes during the seismic cycle using our package. We have also performed benchmarks for our code in comparison to the other community softwares [@Erickson:2023].

# Features
- Antiplane strain deformation on planar fault with SH waves
- Fully dynamic treatment of seismic events
- Rate-and-state-dependent friction with aging laws
- Adaptive time-stepping to switch between interseismic and seismic events
- Customizable to include off-fault and on-fault heterogeneities
- Complex geometry of fault damage zones including gaussian and trapezoid fault damage zones
- Time-dependent healing of fault damage zones constrained by observations
- Precursory seismic velocity changes in the fault zone


# Methodology
We consider a two‐dimensional strike‐slip fault embedded in an elastic medium with Mode III rupture (Figure 1). This implies that the fault motion is out of the plane and only the depth variations of parameters are considered. The top boundary is stress‐free and represents Earth's free surface. The other three boundaries are absorbing boundaries that allow the waves to pass through. Since our model is symmetric across the fault, we restrict the computational domain to only one side of the fault. We model the fault damage zone as an elastic layer with lower seismic wave velocities compared to the host rock. On the fault boundary, we consider a rate- and state-dependent friction law that relates the shear strength on the fault to the
fault slip rate [@Dieterich:1979; @Ruina:1983; @Scholz:1998]. We use the regularized form for the shear strength interpreted as a thermally activated creep model [@Lapusta:2000; @Rice:1996]. 

We use a spectral element method to simulate dynamic ruptures and aseismic creep on the fault [@Kaneko:2011]. Full inertial effects are considered during earthquake rupture and an adaptive time stepping technique is used to switch from interseismic to seismic events based on a threshold maximum slip velocity of 5 mm s−1 on the fault. This fully dynamic modeling approach can simulate interseismic slip, earthquake nucleation, rupture propagation, and postseismic deformation during multiple seismic cycles in a single computational framework. We use a two-dimensional quadrilateral mesh with 5 Gauss-Lobatto-Legendre interpolation points inside each spectral-element. We implement @Kaneko:2011's algorithm in Julia [@Bezanson:2017] using a more efficient linear solver based on the Algebraic Multigrid scheme [@Ruge:1987] for the elliptic (interseismic) part of the earthquake sequence. We use the Algebraic Multigrid as a preconditioner while solving the sparse linear system using the conjugate gradient method. This combines the superior convergence properties of the Algebraic Multigrid with the stability of Krylov methods and is very well suited for symmetric, positive definite matrices. This iterative technique uses a fixed number of iterations independent of the mesh size. @Landry:2016 and @Landry:2019 have derived the equations to solve elliptic equations using the Geometric Multigrid in 2‐D and 3‐D. While the Geometric Multigrid has superior convergence properties, the Algebraic Multigrid is better suited for more complicated meshes and is scalable to a wide variety of problems as the solver works with the numerical coefficients of the linear system as opposed to the mesh structure. The detailed algorithm is described in @Tatebe:1993. In addition, we use the built‐in multithreading feature of Julia. Figure 2 shows the slip contours along the fault and the maximum slip rate along the fault for an example simulation. This example simulation can be run by referring to the USER_GUIDE.md Example 1. 

The software is open-source and available on Github under GNU license, and can be accessed at: https://github.com/thehalfspace/Spear. The instructions to run the simulations are provided in the [README.md](https://github.com/thehalfspace/Spear/blob/main/README.md) and the [USER_GUIDE.md](https://github.com/thehalfspace/Spear/blob/main/USER_GUIDE.md). Open source contributions are welcome from users and can be added under github issues and pull-requests. 


# Figures

![(a) A conceptual model of a strike-slip fault system with a fault damage zone and the seismogenic zone. This model is two-dimensional with along-strike dimension extending to infinity. (b) Two-dimensional modeling framework based on the conceptual model. The fault-zone half-width and the length scales can be adjusted.
\label{fig:example}](conceptual_model.png).

![Example simulation results showing: (a) The accumulated slip along the fault through multiple earthquake sequences. The orange line shows the slip-contours during an earthquake plotted every 1 sec, and the blue lines show the slip contours between earthquakes plotted every 2 years. (b) The maximum rate of slip along the fault through time. Any event above the threshold of 1 mm/s is interpreted as an earthquake, and the events below the threshold can be interpreted as slow-slip events.](results.png)


# List of Publications
Our package has been used in several publications, listed below:
Thakur, P., Huang, Y., & Kaneko, Y. (2020). Effects of low-velocity fault damage zones on long-term earthquake behaviors on mature strike-slip faults. Journal of Geophysical Research: Solid Earth, 125(8).
Thakur, P., & Huang, Y. (2021). Influence of fault zone maturity on fully dynamic earthquake cycles. Geophysical Res. Let.
Erickson, B. A., Jiang, J., Lambert, V., Abdelmeguid, M., Almquist, M., Ampuero, J. P., ..., Thakur, P., ..., & Yang, Y. (2022). Incorporating Full Elastodynamic Effects and Dipping Fault Geometries in Community Code Verification Exercises for Simulations of Earthquake Sequences and Aseismic Slip (SEAS). Seismological Res. Let..
Thakur, P., & Huang, Y. (2024) The effects of precursory velocity changes on earthquake nucleation and stress evolution in dynamic earthquake cycle simulations. Earth and Planetary Science Letters.
Zhai, P., & Huang, Y. (2024). The effects of characteristic weakening distance on earthquake nucleation styles in fully dynamic seismic cycle simulations. Journal of Geophysical Research: Solid Earth, 129(12).


# Acknowledgements
This study was supported by the National Science Foundation (Grant Award EAR-1943742) and JSPS KAKENHI (21H05206).


# References
`docs/refs.bib`
