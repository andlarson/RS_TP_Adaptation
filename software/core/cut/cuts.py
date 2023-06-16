import util.geom as geom



class ToolPass:

    # TODO: 
    # For now, the implementation of a ToolPass is quite brittle. Each
    #   tool pass is made up of exactly 8 vertices. It is assumed that these 8
    #   vertices live in the same coordinate system as the part and that they
    #   make up a rectangular prism. It is also assumed that the rectangular
    #   prism lives entirely within the part and that one face of rectangular
    #   prism is coincident with a face of the part (this essentially means that
    #   each chunk of the part's surface can only be claimed once). Furthermore, 
    #   all passes are assumed to be specified in the coordinate system of the part's 
    #   initial geometry (i.e. before any deformation occurs). This differs
    #   from the way that machine cuts are specified in G-Code ("move the tool
    #   to this location, lower it, ...").
    # 
    # There needs to be a nice relationship between the residual stress field
    #   definition and the tool pass definition. There are a couple ways we can
    #   go about this:
    #   1) Partition the part to define residual stress. Then add more partitions
    #      for the tool passes. 
    #      The downside with this approach is that we will have overlapping 
    #      partitions. The volume of a cut might pass through multiple regions 
    #      which have different residual stress tensors.
    #      This is the way I did things by hand. When doing things by hand it
    #      has the upside that changes to the mesh do not affect part partitons.
    #   2) Partition the part to define residual stress. Then mesh the part. Then
    #      collect elements from the mesh which need to be removed to do the
    #      part cutting.
    #      We might eventually want to mesh the part in a way which depends on
    #      the cuts. This approach might complicate that.
    #      This approach would allow us to potentially use sets in a clever way
    #      since a collection of mesh elements can be bundled into a set.
    #   3) Mesh the part. Then define both residual stress and the part cutting
    #      in terms of the mesh elements.
    #      Some potential problem due to dependency between the mesh and the
    #      cuts.
    def __init__(self, tool_pass):
    # type: (geom.RectPrism) -> None

        self.tool_pass = tool_pass



class ToolPasses:
    
    def __init__(self):
    # type: (None) -> None

        # Keep track of what's been done and what needs to be done. 
        self.passes_done = []
        self.passes_todo = []


    def add_pass(self, tool_pass):
    # type: (ToolPass) -> None

        self.cuts_todo.append(tool_pass)




