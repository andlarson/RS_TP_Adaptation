import util.geom as geom



class ToolPass:

    # TODO: For now, the implementation of a ToolPass is quite brittle. Each
    #   tool pass is made up of exactly 8 vertices. It is assumed that these 8
    #   vertices live in the same coordinate system as the part and that they
    #   make up a rectangular prism. Any portion of the part which exists 
    #   inside the rectangular prism will be removed from the part via element
    #   deletion. It is neither assumed nor required that the rectangular prism
    #   exists entirely within the part. Only the portion of the rectangular
    #   prism intersecting the part will be removed. Furthermore, all passes 
    #   are assumed to be specified in the coordinate system of the part's 
    #   initial geometry (i.e. before any deformation occurs). This may be
    #   problematic if the part geometry changes drastically due to deformation
    #   after only a portion of the passes have been performed, but it 
    #   simplifies things for now.
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




