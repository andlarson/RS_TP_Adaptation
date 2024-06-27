"""
Building a solid sweep via datums -> wire -> solid sweep.
"""

import sys

from abaqus import *
from abaqusConstants import *


def dp(message):
# type: (str) -> None

    sys.__stderr__.write(message + "\n")



def sketch_rectangle(length, width, part_name, model_name, mdb):
# type: (float, float, str, str, Any) -> Any

    sketch = build_constrained_sketch(None, part_name, model_name, mdb)

    # Draw the rectangle.
    sketch.rectangle((width, 0), (-width, length))

    return sketch



def build_constrained_sketch(transform, sketch_name, model_name, mdb):
# type: (Any, str, str, Any) -> None

    if transform is not None:
        sketch = mdb.models[model_name].ConstrainedSketch(sketch_name, sheetSize=1, transform=transform)
    else:
        sketch = mdb.models[model_name].ConstrainedSketch(sketch_name, sheetSize=1)

    if sketch == None:
        raise RuntimeError("Failed to build sketch!!")

    return sketch



if __name__ == "__main__":

    mdb = openMdb("test.cae")

    part = mdb.models["Model-1"].parts["Part-1"]

    datums = []

    feature = part.DatumPointByCoordinate((1, 2, 3))
    datums.append(part.datums[feature.id])
    first_datum = part.datums[feature.id]

    feature = part.DatumPointByCoordinate((4, 5, 3))
    datums.append(part.datums[feature.id])

    feature = part.DatumPointByCoordinate((3, 2, 3))
    datums.append(part.datums[feature.id])
    last_datum = part.datums[feature.id]

    feature = part.WireSpline(datums, mergeType=SEPARATE)

    edge = part.edges[len(part.edges)-1:]

    # Length = .05
    # Radius = .03
    sketch = sketch_rectangle(.05, .03, "Part-1", "Model-1", mdb)

    df1 = part.DatumPointByCoordinate((0, 0, 0))
    df2 = part.DatumPointByCoordinate((0, 0, 1))
    datum_axis_feature = part.DatumAxisByTwoPoint(part.datums[df1.id], part.datums[df2.id])
    datum_axis = part.datums[datum_axis_feature.id]

    mdb.saveAs("post_test.cae")

    sweep = part.SolidSweep(path=edge, profile=sketch, sketchUpEdge=datum_axis, sketchOrientation=RIGHT)

    # At the first datum and last datum, we need to cap off the toolpath.

    """
    face_at_first_datum = part.faces.findAt(first_datum.pointOn)
    face_at_last_datum = part.faces.findAt(last_datum.pointOn)

    edge_ids = face_at_first_datum.getEdges()
    edge_on_first_face = part.edges[edge_ids[0]]
    """

    mdb.saveAs("post_test.cae")
    mdb.close()