Significant changes to the software architecture are currently being made. These
are being carried out on the `rewrite` branch. 

It turns out that using Abaqus' preprocessor alone, via Abaqus' Python API, for
geometric manipulations is entirely insufficient for the goals of this project.
Other than being generally difficult to work with (e.g., Abaqus' Python API can
only be accessed through a custom Python interpreter, the only way to debug code
is a broken GUI-based debugger, etc.), Abaqus' preprocessor:
1. Doesn't offer meshing routines that are robust without humans in the
   loop.
2. Doesn't allow the user to control how surface meshes are generated on
   boundary representations (i.e. the .stl export functionality is a black
   box).
3. Requires expensive, unnecessary, and error-prone conversions between boundary
   representation of geometric objects and discretized representation of
   geometric objects.
4. Makes it hard to do canonicalization. 

The new software architecture, being developed on the `rewrite` branch, uses a
suite of specialized libraries:
1. OCCT, an open source geometry kernel.
2. CGAL, an open source computational geometry algorithms library.
3. fTetWild, an open source robust meshing algorithm.
4. Abaqus, a closed source FEA software. 
You won't find uses of OCCT and fTetWild in this repository, but they are used
in other repositories to prepare the inputs to this repository.
