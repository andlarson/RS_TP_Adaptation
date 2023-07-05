import util.geom as geom

""" 
There needs to be a nice relationship between the residual stress field
  definition and the tool pass definition. There are a couple ways we can
  go about this:
  1) Partition the part to define residual stress. Then add more partitions
     for the tool passes. 
     The downside with this approach is that we will have overlapping 
     partitions. The volume of a cut might pass through multiple regions 
     which have different residual stress tensors.
     This is the way I did things by hand. When doing things by hand it
     has the upside that changes to the mesh do not affect part partitons.
  2) Partition the part to define residual stress. Then mesh the part. Then
     collect elements from the mesh which need to be removed to do the
     part cutting.
     We might eventually want to mesh the part in a way which depends on
     the cuts. This approach might complicate that.
     This approach would allow us to potentially use sets in a clever way
     since a collection of mesh elements can be bundled into a set.
  3) Mesh the part. Then define both residual stress and the part cutting
     in terms of the mesh elements.
     Some potential problem due to dependency between the mesh and the
     cuts.
  4) Create a single part for the initial geometry and a part for each cut.
     Create instances for each of these in the Assembly module. Use the
     cut functionality in the Assembly module to create the post-cut geometry,
     which still has a residual stress pattern in its remaining volume.
     This approach deals with arbitrary cuts in space pretty well. It
     natively removes only overlapping volume. 

A fundumental realization:
Cuts are specified as arbitrary paths in space. However, after the first cut
is made the geometry of the part might change due to deformation. This
means that the portion of material which is removed for a single cut depends
on the deformation which occurred during the previous cut. Therefore, under
the assumption of non-negligible deformation, each cut requires a single 
simulation and the simulations must be chained together.
""" 


class ToolPass:
    
    def __init__(self, tool_pass):
    # type: (geom.SpecRightRectPrism) -> None

        self.tool_pass = tool_pass



class ToolPassPlan:
    
    def __init__(self, tool_passes):
    # type: (list[ToolPass]) -> None

        self.passes_done = []
        self.passes_todo = tool_passes


    def add(self, tool_pass):
    # type: (ToolPass) -> None

        self.passes_todo.append(tool_pass)

    
    def pop(self):
    # type: (None) -> ToolPass 
        
        tool_pass = self.passes_todo.pop(0)
        self.passes_done.append(tool_pass)

        return tool_pass

