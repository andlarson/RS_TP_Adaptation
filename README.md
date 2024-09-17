## Improving Machining Yield in the Presence of Residual Stress ##

This repository contains software that makes it easy to:
1. Use surface meshes, *created from point clouds captured while machining a 
   workpiece*, to estimate the residual stress field contained within the 
   workpiece. 

This software supports the larger, project-level goal of improving machining 
yield for blanks which have problematic residual stress fields (i.e. residual
stress fields that cause the workpiece to deform during and after machining
into a geometry that does not satisfy allowed tolerances).

## Dependencies ##

- TBD

## Build Instructions ##

```
mkdir build
cd build
cmake ..
cmake --build .
./rs_tp_adaptation --help 
```
