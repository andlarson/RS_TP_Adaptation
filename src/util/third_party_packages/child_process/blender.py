"""
This file contains functionality for making calls into Blender's Python APIs.
"""



def remesh(stl_to_remesh: str, save_location) -> None:
    """Remeshes a mesh with some heuristically chosen new attributes (i.e.
           number of vertices in new mesh, number of faces in new mesh,
           etc.). 
                   
       Args:
           stl_to_remesh: Absolute path to .stl file to remesh.
           save_location: Absolute path to desired save location. Should end
                              with .stl.
    
       Returns:
           None.

       Raises:
           None.
    """
    
    return



def compute_volume_symmetric_difference(stl1: str, stl2: str) -> float:
    """Computes the volume of the symmetric difference between two 3D shapes.
                   
       Args:
           stl1: Absolute path to .stl file containing a shape.
           stl2: Absolute path to .stl file containing a shape.
    
       Returns:
           Volume of symmetric difference.
    
       Raises:
           None.
    """

    return -1



def remove_material(base_geom: str, to_remove: str, save_loc: str) -> None:
    """Removes a volume of material. 

       TODO: To utilize full capability of Blender, this function should accept
           a shape and a tool pass path (i.e. a spline or the like). Then Blender
           native operations could be used to do the removal (constructing the
           shape of the tool pass path only if necessary).
                   
       Args:
           base_geom: Absolute path to .stl file containing the base geometry from
                          which material will be removed.
           to_remove: Absolute path to .stl file containing the geometry to remove
                          from the base geometry.
           save_loc:  Absolute path to desired save location of .stl which contains
                          the geometry post-material removal. Should end with .stl.

       Returns:
           None.
    
       Raises:
           None.
    """

    return 




