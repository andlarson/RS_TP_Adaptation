"""
Let's see if an Abaqus Feature object is truly useless...
"""

import sys

from abaqus import *
from abaqusConstants import *


def dp(message):
# type: (str) -> None

    sys.__stderr__.write(message + "\n")



if __name__ == "__main__":

    mdb = openMdb("test.cae")

    part = mdb.models["Model-1"].parts["Part-1"]

    dp("There are currently " + str(len(part.features)) + " features associated with the part.")
    dp("There are currently " + str(len(part.datums)) + " datums associated with the part.")
    dp("The type of the container for a part's datums is " + str(type(part.datums)))

    feature = part.DatumPointByCoordinate((1, 2, 3))

    dp("There are currently " + str(len(part.features)) + " features associated with the part.")
    dp("There are currently " + str(len(part.datums)) + " datums associated with the part.")

    dp("The id of the new feature is " + str(feature.id))
    dp("The name of the new feature is " + str(feature.name))

    feature = part.DatumPointByCoordinate((4, 5, 6))

    dp("There are currently " + str(len(part.features)) + " features associated with the part.")
    dp("There are currently " + str(len(part.datums)) + " datums associated with the part.")

    dp("The id of the new feature is " + str(feature.id))
    dp("The name of the new feature is " + str(feature.name))

    dp("The content of the feature rep for the part is " + str(part.features))
    dp("The content of the datum repository for the part is " + str(part.datums))

    dp("Now, let's try to access the Datum repo using the feature is....")
    datum = part.datums[feature.id]

    mdb.saveAs("post_test.cae")
    mdb.close()