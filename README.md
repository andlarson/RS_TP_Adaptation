# Mitigating workpiece warpage caused by residual stress during machining!

---
### Current Assumptions / Limitations:
1. Toolpath geometries must be special right rectangular prisms.
2. Mesh density is arbitrary.
3. No boundary conditions are imposed for any simulations.
4. Any residual stress profile specified by the user must be in equilibrium. 
5. The user provides the initial part geometry as a .cae file. 
6. There is no functionality to incorporate real-life measurements.
7. The only material properties that matter are Young's Modulus and Poisson's Ratio.
8. **There is no explicit functionality to support clamping.** 
9. When a full machining process is conducted in simulation, the user must only specify the part's initial residual stress profile. The stress profile which results from one tool pass simulation is always used as the initial stress profile for the next tool pass simulation. 
---
