# Instructions

Do this group assignment in groups of 3. You can make the groups yourself and/or use the Discussion board to find a group. This assignment counts for 15% of your grade.

Enroll for Group Assignments in the menu to create your group first!

# Overview

In this assignment you are asked to program a BEM (Blade Element Momentum) model for one of two rotors:

 - A wind turbine, in axial and yawed flow
 - A propeller, in axial and yawed flow, in cruise (optional: and energy harvesting)  

NOTE: choose one of the two rotor cases for the mandatory assignment, If you wish, you can do both rotors.

The BEM code should incorporate the Glauert correction for heavily loaded rotors for the wind turbine case, and the Prandtl tip and root corrections. Assume a constant airfoil along the span and use the provided polars of the:

 - DU airfoil for the wind turbine
 - ARA-D 8% airfoil for the propeller

# Tasks to be done

- Carlos
  - A flowchart of your code (1 page)
  - Plots with explanation of the influence of the tip correction
  - Plot the distribution of stagnation enthalpy as a function of radius at four locations: infinity upwind, at the rotor (upwind side), at the rotor (downwind side), infinity downwind.
  - Plot a representation of the system of circulation. Discuss the generation and release of vorticity in relation to the loading and circulation over the blade.
- Bernat
  - SHORT introduction (1 page)
  - Plots with explanation of results (alpha/inflow/a/a’/Ct/Cn/Cq vs r/R)
    - Spanwise distribution of angle of attack and inflow angle
    - Spanwise distribution of axial and azimuthal inductions
    - Spanwise distribution of thrust and azimuthal loading
    - Total thrust and torque
  - Plots with explanation of influence of number of annuli, spacing method (constant, cosine) and convergence history for total thrust.
  - Discuss the different inaccuracies introduced into the solution by using a single airfoil polar only.
- Niklas
  - Main assumptions with an explanation of their impact
  - For the cases of yawed HAWT, also plot the azimuthal variation (suggestion: polar contour plot)
  - Discuss the operational point of the airfoil in terms of lift and drag, and the relation between the distribution of lift coefficient and chord (chose one case, not necessary to do all cases)
  - SHORT discussion/conclusion, including the similarities and differences between the two rotor configurations, flow field and operation

## Pending
- (optional) Explanation of the design approach used for maximizing the Cp or efficiency
- (optional) Plots with explanation of the new designs
