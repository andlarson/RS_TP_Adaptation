"""
This file contains functionality for making calls into Blender's Python APIs.
"""



def remesh(message: str) -> str:
    """Remeshes a mesh with some heuristically chosen new attributes (i.e.
           number of vertices in new mesh, number of faces in new mesh,
           etc.). 
                   
       Args:
           message: Unmodified message data.
    
       Returns:
           Absolute path to .stl file containing the remeshed mesh.

       Raises:
           None.
    """
    
    return "Received instruction to remesh!"



def compute_volume_symmetric_difference(message: str) -> int:
    """Computes the volume of the symmetric difference between two 3D shapes.
                   
       Args:
           message: Unmodified message data.
    
       Returns:
           Volume of the symmetric difference.
    
       Raises:
           None.
    """

    return -1



def remove_material(message: str) -> str:
    """Removes a volume of material. 

       TODO: To utilize full capability of Blender, this function should accept
           a shape and a tool pass path (i.e. a spline or the like). Then Blender
           native operations could be used to do the removal (constructing the
           shape of the tool pass path only if necessary).
                   
       Args:
           message: Unmodified message data.

       Returns:
           Absolute path to .stl file containing the shape resulting from the
               removal.
    
       Raises:
           None.
    """

    return "Received instruction to remove material."




