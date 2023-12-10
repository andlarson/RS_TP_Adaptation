"""
Testbed for small components of the implementation.
"""

import sys

# DEBUG
# For Abaqus PDB and allows main.py to be run from any directory.
# Imperfect solution, but good enough for now. Every time this file is imported,
#    this line runs. This pollutes the sys.path list.
sys.path.append("/home/andlars/Desktop/RS_TP_Adaptation/src")

import numpy as np

import util.geom as geom

l1 = [(40.0, 10.0, 40.0), (0.0, 10.0, 40.0), (0.0, 10.0, 0.0), (40.0, 10.0, 0.0)]
l2 = [(0.0, 10.0, 0.0), (40.0, 10.0, 0.0), (40.0, 10.0, 40.0), (0.0, 10.0, 40.0)]

s1 = geom.seq_points(l1)
s2 = geom.seq_points(l2)

print(geom.seq_points_equal(s1, s2))
