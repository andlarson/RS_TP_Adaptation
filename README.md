# Mitigating workpiece warpage caused by residual stress during machining!

---
### Current Assumptions / Limitations:
1. Toolpath geometries must be special right rectangular prisms.
2. Mesh density is arbitrary.
3. No boundary conditions are imposed for any simulations.
4. Ordering between material removal and stress equilibration does not matter.
5. The user provides the initial part geometry as a .cae file. No functionality to incorporate real-life measurements.
6. The only material properties that matter are Young's modulus and Poisson's Ratio.
7. **No clamping is assumed.**
8. User must specify initial RS profile via user subroutine.
---
