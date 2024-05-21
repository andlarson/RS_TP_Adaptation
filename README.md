---
# !!!UPDATE!!! #

Significant changes to the software architecture are necessary. The functionality 
that we desire is not natively supported by Abaqus. In particular, Abaqus' meshing
routines are not robust, material removal is not natively supported in Abaqus
(and workarounds to do material removal in Abaqus do not fulfill our needs), and 
we need custom meshing algorithms to support canonicalization. As such, we are
in the process of rewriting a large portion of this repository.

---



---

# Improving Machining Yield in the Presence of Residual Stress #

This repository contains software that makes it easy to:
1. Estimate the residual stress state of a workpiece by utilizing real-world
   measurement data collected after the workpiece has deformed due to machining.
2. Simulate sequences of potential tool passes. 

This software supports the larger, project-level goal of improving machining 
yield for blanks which have problematic residual stress fields (i.e. residual
stress fields that cause the workpiece to deform during and after machining
into a geometry that does not satisfy allowed tolerances).

---



---

# Project File Structure # 

## TODO ##

---





