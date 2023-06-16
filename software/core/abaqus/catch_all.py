
'''
For now, this is a catch-all file which interfaces with the Abaqus scripting
  interface.
This file will probably need to be broken up by functionality.
'''

from abaqus import *
from abaqusConstants import * 


# Build a part in the current Abaqus model database.
# TODO: For now, we assume the part is a rectangular prism. This certainly won't
#   be the case in the future for either initial geometries or tool paths.
def build_part(rect_prism):
# type: (RectPrism) -> None


