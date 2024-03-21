This directory contains motivatation for the use of Blender.

The high level rationale is that Abaqus' preprocessor utilities aren't very robust.

More precisely, Abaqus is not good at:
    1) Converting a mesh to a geometry.
    2) Doing boolean operations between complex geometries.
    3) Representing all possible geometries of tool passes.



********************************************************************************
                        Converting Mesh to Geometry
********************************************************************************

At the time of writing this, I have sort of *observed* that a bottleneck in the
    software is converting a mesh to a geometry. It is necessary to go from a 
    mesh to a geometry in the following situations:
        1) When consecutive tool paths are simulated, a deformed mesh needs to
               be mapped to a geometry for a subsequent FEM analysis to
               occur.
        2) When data is captured from the real world, the conversion process
               goes like pointcloud -> .stl -> geometry to enable FEM analysis.
               Thus, a mesh to geometry conversion is necessary.
        3) When simulating real life, it's necessary to convert a deformed mesh
               to a geometry because the "real life data" is just a deformed mesh
               produced by a hidden FEM analysis.
I hypothesize that deformed meshes (i.e., meshes produced by a FEM simulation)
    can be non-manifold (i.e., not define a closed geometry), causing a failure in
    conversion to a geometry. I've also observed that when failure occurs in 
    converting from a mesh to a geometry, it usually occurs in regions where 
    there are lots of small features (i.e., small faces). My most precise hypothesis
    is that a deformed mesh can be non-manifold due to the rounding of floats.

At the time of writing this, I haven't convinced myself that my hypothesis is
    actually true. Thus, I've created this directory to convince myself. 

This directory contains a .odb which resulted from a simulation. When converting
    the deformed mesh in this .odb to a geometry using Abaqus' built in AddCells()
    method, the conversion fails! However, by default my implementation uses the 
    more robust Stitch() method. Using the Stitch() method, the conversion works. 
    However, I know that the Stitch() method can also fail....I just can't seem 
    to find a .odb which showcases the failure.



********************************************************************************
                 Boolean Operation Between Complex Geometries 
********************************************************************************

We sometimes want to do boolean operations between complex geometries. For example,
    when doing 3D shape comparison we want to find the symmetric difference
    between two arbitrary 3D shapes.
Abaqus has the ability to do boolean operations between geometries. However, as
    shown in failed_boolean_intersection.cae the boolean intersection operation
    is quite unreliable.



********************************************************************************
               Representing All Possible Geometries of Tool Passes 
********************************************************************************

Abaqus restricts the geometries which it can represent. For example, Abaqus restricts
    wires to be cubic splines with continuous first and second derivatives. This
    restriction means that we can't represent all possible geometries of tool
    passes that the machine is capable of.



