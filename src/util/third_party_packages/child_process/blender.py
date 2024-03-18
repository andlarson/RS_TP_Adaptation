"""
This file contains functionality for making calls into Blender's Python APIs.
"""



def remesh(path: str, save_dir_path: str) -> str:
    """Remeshes a mesh with some heuristically chosen new attributes (i.e.
           number of vertices in new mesh, number of faces in new mesh,
           etc.). 
                   
       Args:
           path:          Absolute path to .stl file containing mesh to remesh.
           save_dir_path: Absolute path to directory in which intermediate/result
                              files should be saved.
    
       Returns:
           Absolute path to .stl file containing the remeshed mesh.

       Raises:
           None.
    """




def compute_volume_symmetric_difference(path1: str, path2: str, 
                                        save_dir_path: str) -> int:
    """Computes the volume of the symmetric difference between two 3D shapes
           described by .stl file.
                   
       Args:
           path1:         Absolute path to .stl file describing a shape.
           path2:         Absolute path to .stl file describing a shape.
           save_dir_path: Absolute path to directory in which intermediate/result
                              files should be saved.
    
       Returns:
           Volume of the symmetric difference.
    
       Raises:
           None.
    """



def remove_material(path1: str, path2: str, save_dir_path: str) -> str:
    """Removes the material volume described by the 3D shape in the .stl file
           at path2 from the material volume described by the 3D shape in the
           .stl file at path1.

       TODO: To utilize full capability of Blender, this function should accept
           a shape and a tool pass path (i.e. a spline or the like). Then Blender
           native operations could be used to do the removal (constructing the
           shape of the tool pass path only if necessary).
                   
       Args:
           path1:         Absolute path to .stl file describing shape to remove 
                              volume from.
           path2:         Absolute path to .stl file describing shape of volume to
                              remove.
           save_dir_path: Absolute path to directory in which intermediate/result
                              files should be saved.
    
       Returns:
           Absolute path to .stl file containing the shape resulting from the
               removal.
    
       Raises:
           None.
    """
