# Improving Machining Yield in the Presence of Residual Stress

This repository contains software that makes it easy to:
1. Estimate the residual stress state of a workpiece by utilizing data about 
   its deformation.
2. Simulate sequences of potential tool passes. 

This software supports the larger, project-level goal of improving machining 
yield for blanks which have problematic residual stress fields (i.e. residual
stress fields that cause the workpiece to deform during and after machining
into a geometry that does not satisfy allowed tolerances).


---

## Simulating Sequences of Tool Paths ##

### Current Assumptions and Limitations:

1. The space of potential toolpaths is exactly the set of toolpaths which can be 
   taken by the machine. The set of toolpaths which can be taken by the machine 
   is governed by the G-code which the machine accepts. Some machines are only 
   able to perform, for example, linear interpolation and circular interpolation. 
   Others support toolpaths which are specified via cubic splines, quadratic splines, 
   and more exotic options.

   For full generality, the simulation engine must be able to represent and 
   simulate all toolpaths which the machine can perform. Furthermore, the toolpath 
   representations cannot be approximate representations representations. They 
   must be as **exact** as possible. 

   The simulation engine is currently able to represent toolpaths created by 
   cylindrical tools which follow cubic splines in 3D and are **not self-intersecting**. 
   Note that multiple toolpaths may intersect one another, but a single toolpath 
   should not intersect itself. Furthermore, the tool orientation with respect 
   to the toolpath is currently fixed and unchangeable. For all toolpaths, the 
   axis of rotational symmetry of the cylindrical tool is fixed and positioned 
   parallel to the y axis. See `notes/toolpath/toolpath_orientation_1.jpg` and 
   `notes/toolpath/toolpath_orientation_2.jpg` for more information.
   
2. Meshes are always constructed with tetrahedrons and the Free meshing technique. 
   Starting from some heuristically chosen density, the mesh is made more dense 
   until the part is successfully meshed. Currently, mesh quality is not considered.

3. TODO: 
   All simulations include an initial step which allows any residual stress 
   profile to relax to equilibrium. 

   If everything is being done in simulation (i.e. there is no real-world machining 
   process going on), there are exactly three sources of residual stress profiles: 
   the user-specified initial residual stress profile, the residual stress profile 
   which remains in a part after some chunk of the part has been removed during 
   simulation, and estimated residual stress profiles based on the deformations 
   which occur in simulation. If there is a real-world machining process going on, 
   the residual stress profile estimates come from deformations measured in 
   real-life, not in simulation. 
   
   In any case, it's possible that a user-specified residual stress profile or 
   an estimated residual stress profile does not satisfy mechanical equilibrium. 
   In general, a user will probably try to specify an initial residual stress 
   profile which is in equilibrium and any residual stress profile estimation 
   technique will probably attempt to produce an estimated residual stress profile 
   which does satisfy mechanical equilbirium. However, to mitigate against the 
   possibility of residual stress profiles which do not satisfy mechanical equilibrium, 
   a relaxation step is included as the first step in all simulations. During 
   this step, any residual stess profile which is not in equilibrium is allowed 
   to relax to equilibrium. If the residual stress profile already satisfies 
   mechanical equilibrium, nothing happens. If the residual stress profile does 
   not satify mechanical equilibrium, then some deformation will occur. For now, 
   we assume that such deformations are small and have negligible effect on the 
   fidelity of our simulations.
   
4. The user provides the initial part geometry via a .cae file. 

5. The clamping setup is represented as boundary conditions which exist in the 
   same coordinate system as the initial part geometry. For now, we only allow 
   displacement boundary conditions (restricting at most displacement in x, y, 
   and z directions and rotation in x, y, and z directions). We also expect that 
   the boundary conditions are applied to surfaces without curvature. This means 
   that we only allow clamping setups where the clamps are applied to surfaces 
   without curvature.

6. The only material properties which we care about are Young's Modulus and 
   Poisson's Ratio.

7. The user has the option to inject their own estimate of the entire residual 
   stress profile of a part after a tool pass plan has been committed (and potentially 
   conducted in real life). However, by default, the stress profile which resulted 
   from the simulation of the last committed tool pass plan is the stress profile 
   used as the starting stress profile of the part in the next commitment phase. 
---


---
## Estimating Residual Stress via Deformation Measurements

TODO....

---

