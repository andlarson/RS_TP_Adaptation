---

# UPDATE #

Significant changes to the software architecture and implementation are necessary. 
The functionality that we desire is not natively supported by Abaqus. In particular, 
Abaqus' meshing routines are not robust, material removal is not natively 
supported in Abaqus (and workarounds to do material removal in Abaqus do not 
fulfill our needs), and we need custom meshing algorithms to support 
canonicalization. As such, we are in the process of rewriting a large portion of 
this repository. This is primarily being done on the `migrate_to_cpp` branch.

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


---

# Dependencies #

This project assumes that its dependencies are available on the user's system.
Precisely, it expects that its dependencies are available 

The software depends on:
1. Header files for the Geogram, Threads, `TODO: Add more as necessary.` libraries.
2. Static libraries including:
   fTetWild. Precisely, the static library produced when fTetWild is built. This
   static library contains code from fTetWild, Geogram, GMP, IGL, Threads, etc. 

---



---

# Build Instructions #

1. Create a `build/` directory under `RS_TP_Adaptation/`.
2. `cd build/`
3. `cmake ..`
4. `cmake --build .`
5. Now the `build/` directory should contain the executable `rs_tp_adaptation`. Run the executable by issuing `./rs_tp_adaptation`.

Right now, automated installation is not provided. 

---


