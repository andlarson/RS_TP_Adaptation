Is Blender really necessary?
The high level rationale is that Abaqus' preprocessor utilities aren't very 
    robust. This rationale is a bit fuzzy, so I want to be more precise.

More precisely, Abaqus does not seem to be good at: 
    1) Converting a mesh to a geometry.
    2) Doing boolean operations between complex geometries.
    3) Representing all possible geometries of tool passes.

There is also further motivation for remeshing, which Blender is capable of. Assume
    that a mesh -> geometry conversion were possible without remeshing. In general,
    the resulting geometry will have many small faces - the faces will match those
    of the mesh, and those faces are certainly nonuniform. These small faces on
    the geometry will make meshing the geometry quite difficult. By remeshing
    with a more uniform mesh and then converting to a geometry, the mesh
    generation process ought to be much easier. 



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
    there are lots of small features (i.e., small faces). After asking around
    online, other folks seem to think that my hypothesis are not true, but they
    haven't offered any better suggestions...

This directory contains a .odb which resulted from a simulation. When converting
    the deformed mesh in this .odb to a geometry using Abaqus' built in AddCells()
    method, the conversion fails! However, by default my implementation uses the 
    more robust Stitch() method. Using the Stitch() method, the conversion works. 
    However, I know that the Stitch() method can also fail....I just can't seem 
    to find a .odb which showcases the failure.

Debugging:
I think that the biggest hint is that the Stitch() method leads to success while
    the simple AddCell() method fails. The Stitch() method repairs gaps between
    edges. The upshot of this is that, maybe, a deformed mesh has tiny gaps
    between edges, which causes a conversion to a geometry to be impossible
    without stitching. To test this hypothesis, I'll do the following: 
        1) I'll import the undeformed mesh from the .odb into a fresh .cae, 
               create a geometric face for each face in the mesh, attempt to 
               create a geometry without using stitching (which should succeed 
               because this is the undeformed mesh after all), and if the geometry 
               cannot be created I'll export the collection of faces. 
        2) Next, I'll do the same thing as (1) but this time, I'll use the
               deformed geometry. In this case, I know that it won't be possible
               to create a geometry without using stitching.
    If, as I expect, test (1) results in a successful conversion from mesh to
        geometry and test (2) results in an unsuccessful conversion from mesh
        to geometry, I will import the .stl of the deformed geometry into
        Blender and look for gaps between edges defining the faces. 

Result From Test (1):
Surprisingly, converting an undeformed mesh to a geometry also fails without
    stitching! This is unexpected. In particular, when I query the part (which
    is just composed of a bunch of faces) geometry, it tells me that all the
    vertices are invalid entites and some vertices are imprecise. Note that
    there is no real good defintion of invalid entities - the documentation
    loosley defines an invalid entity as something which makes it impossible
    to create a closed volume. An imprecise entity refers to some feature which
    forces Abaqus to use some tolerance around the entity to create a closed
    volume. Two possibilites occur to me: my routine that maps from mesh -> 
    geometry is incorrect or the floating point data in the output database is not 
    sufficiently precise. 
I went back to the model that produced the offending .odb, wrote out the job
    instructions to the .inp file, and then invoked the job from the command line
    using the "abaqus job=*.inp user=*.o output_precision=full interactive"
    incantation. The necessary files for doing this are in the
    floating_point_hypthesis directory. I then re-ran test (1) using the resulting
    .odb file, producing the file called full_precision.cae. Sadly, this .cae
    file still contains malformed geometry - all the vertices are still invalid
    entities, many of the vertices are imprecise, and all of the edges are
    free edges.
The next possibility to consider is that my routine for doing the mesh ->
    geometry conversion is somehow incorrect. All material for this hypothesis
    is located in the implementation_error directory. To test this, I'll first 
    look at the plugin provided by Abaqus. The compiled and decompiled versions
    of this plugin are located in the abaqus_mesh_to_geom_plugin directory.
    Sadly, the decompiled code seems to use some encryption mechanism to hide
    the source....What's frustrating is that there are many different ways, each
    full of different heuristics, to do the mesh -> geometry conversion. Here
    are the different things I have tried:
    
        1) Create a face for each element face, stitch faces as they are
               created, don't do overall stitching. Attempt shell -> solid
               conversion. VERY slow! 
           Result: stitching_every_face.cae
           Shell -> solid conversion works. There are no invalid entities
               and there are some imprecise entities.
        2) Create a face for each element face, don't stitch faces as they
               are created, do overall stitching. Attempt shell -> solid
               conversion.
           Result: overall_stitch.cae
           Shell -> solid conversion works. There are no invalid entities
               and there are some imprecise entities.
        3) Create a face for each element face, don't do any stitching at
               all. Attempt shell -> solid conversion.
           Result: no_stitch.cae
           Shell -> solid conversion fails. All vertices are invalid entites
               and there are some imprecise entities.
        4) Use the block box plug-in provided by Abaqus. 
           Result: via_abaqus_plugin.cae
           Shell -> solid conversion works. Succeeds and was actually quite
               fast. There are no imprecise or invalid entities.

Conclusion:
I still can't precisely say why (3) fails while (1), (2), and (4) succeed. I also
    can't precisely describe why (1) and (2) result in imprecise entities
    while (3) results in no imprecise entities.
Nonetheless, it seems like the plug-in is the way to go...



********************************************************************************
                 Boolean Operation Between Complex Geometries 
********************************************************************************

We sometimes want to do boolean operations between complex geometries. For example,
    when doing 3D shape comparison we want to find the symmetric difference
    between two arbitrary 3D shapes.
Abaqus has the ability to do boolean operations between geometries. failed_boolean_
    intersection.cae showcases a failure due to an attempted boolean opertion
    between geometries which resulted from mesh -> geometry conversions via
    a overall stitch. boolean_intersection_geoms_from_plugin.cae showcases a
    successful boolean operation between geometries which resulted from mesh ->
    geometry conversions via the plug-in. Evidently, the internal representation
    of the geometry affects the success of boolean operations.

Conclusion:
When the mesh -> geometry plugin provided by Abaqus is used, the boolean capabilites
    in Abaqus seem sufficient.



********************************************************************************
               Representing All Possible Geometries of Tool Passes 
********************************************************************************

Abaqus restricts the geometries which it can represent. For example, Abaqus restricts
    wires to be cubic splines with continuous first and second derivatives. This
    restriction means that we can't represent all possible geometries of tool
    passes that the machine is capable of.

Conclusion:
Eventaully we will want to represent all possible geometries of tool passes, but
    for now, the functionality in Abaqus is sufficient.



********************************************************************************
               Remeshing Complex Geometries Defined by Many Faces
********************************************************************************

When simulating sequential tool passes, converting in-machine scans to geometries,
    and converting simulated in-machine scans to geometries, it's necessary to
    mesh a geometric representation that is highly faceted.
All experimental results are in the meshing_highly_faceted_geometries/ subdirectory.
    The .odb files in this directory contain the example meshes which were used.

The experiments in meshing_highly_faceted_geometries/simple_approach/ were conducted 
    by using a simple procedure:
    1) Import a .odb containing a deformed mesh.
    2) Convert the mesh to a geometry with the mesh -> geometry plugin provided
           by Abaqus.
    3) Mesh the resulting geometry. 

Alternative approaches might include using the virtual topology toolset to simplify
    the geometry and/or doing remeshing to make the mesh more uniform.

The experiments in meshing_highly_faceted_geometries/alternative_approach/ were
    conducted by using the alternative procedure:
    1) Import a .odb containing a deformed mesh.
    2) Convert the mesh to a geometry with the mesh -> geometry plugin provided
           by Abaqus.
    3) Use the default virtual topology settings to simplify the geometry. Do
           this recursively until no further simplification is possible.
    3) Mesh the resulting geometry. 

Conclusion:
A major downside of the simple procedure is that it produces a mesh with a very
    large number of faces. The virtual topology utility helps to mitigate this
    problem, although it does not provide any guarantees.
