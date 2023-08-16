# Mitigating workpiece warpage caused by residual stress during machining!

---
### Current Assumptions / Limitations:
1. Toolpath geometries must be special right rectangular prisms.
2. Meshes are always constructed with tetrahedrons and the Free meshing technique. Starting from some heuristically chosen density, the mesh is made more dense until the part is successfully meshed. Currently, there is no concern about mesh quality. 
3. All simulations include an initial step which allows any residual stress profile to relax to equilibrium. 

   If everything is being done in simulation (i.e. there is no real-world machining process going on), there are exactly three sources of residual stress profiles: the user-specified initial residual stress profile, the residual stress profile which remains in a part after some chunk of the part has been removed during simulation, and estimated residual stress profiles based on the deformations which occur in simulation. If there is a real-world machining process going on, the residual stress profile estimates come from deformations measured in real-life, not in simulation. 
   
   In any case, it's possible that a user-specified residual stress profile or an estimated residual stress profile does not satisfy mechanical equilibrium. In general, a user will probably try to specify an initial residual stress profile which is in equilibrium and any residual stress profile estimation technique will probably attempt to produce an estimated residual stress profile which does satisfy mechanical equilbirium. However, to mitigate against the possibility of residual stress profiles which do not satisfy mechanical equilibrium, a relaxation step is included as the first step in all simulations. During this step, any residual stess profile which is not in equilibrium is allowed to relax to equilibrium. If the residual stress profile already satisfies mechanical equilibrium, nothing happens. If the residual stress profile does not satify mechanical equilibrium, then some deformation will occur. For now, we assume that such deformations are small and have negligible effect on the fidelity of our simulations.
4. The user provides the initial part geometry via a .cae file. 
5. The clamping setup is represented as boundary conditions which exist in the same coordinate system as the initial part geometry. For now, we only allow displacement boundary conditions (restricting at most displacement in x, y, and z directions and rotation in x, y, and z directions).

6. We assume that the clamping setup does not affect the way that the residual stress profile evolves in a part during the machining process. More concretely, it's easy to imagine a clamping setup which restricts the part in such a way that deformations due to residual stress stress don't happen until the clamps are released. For now, we assume that all deformations which occur due to residual stress are observable and not impacted by clamping. Under this assumption, it follows that we need only estimate the residual stress profile which existed in the part before any cuts occurred.
7. There is no functionality to incorporate real-life measurements.
8. For our purposes, the only material properties which we care about are Young's Modulus and Poisson's Ratio. 
9. When a full machining process is conducted in simulation, the user must only specify the part's initial residual stress profile. The stress profile which results from one tool pass simulation is always used as the initial stress profile for the next tool pass simulation. 
---
