"""
Figuring out what linear relationships exist for a linear elastic FEA setup. 
In particular, I want to figure out if there is a linear relationship between
   an applied force vector and a displacement vector. It's pretty obvious to
   me that the linear relationship must depend on position.

Conclusion:
When purely elastic behavior is modeled, there is a linear relationship between
   applied traction and displacement on a per-node basis. That is, for each node,
   there is a 3x3 matrix which relates an applied traction vector to the
   displacement of the node. The matrix differs for each node!
"""

# Resolve imports.
sys.path.append("/home/andlars/Desktop/RS_TP_Adaptation/src")

import sys
import random

import numpy as np
import numpy.linalg as linalg
from abaqus import *
from abaqusConstants import *
import odbAccess
import regionToolset

import core.abaqus.abaqus_shim as shim
from util.debug import *


# Modify an MDB so that a traction of a particular magnitude in a particular 
#    direction is applied. Then run the job.
# The direction is expected to be a 1D numpy array of length 3 which is normalized.
def apply_force(mdb, mag, direction, step_name, job_name):

    # Delete all pre-existing loads.
    dp("Deleting all loads in the model...")
    for name in mdb.models["Model-1"].loads.keys():
        del mdb.models["Model-1"].loads[name]

    # Make sure there is only a single datum associated with the part.
    assert(len(mdb.models["Model-1"].parts["Part-1"].datums) == 1)

    # Extract the id of the datum.
    datum_id = mdb.models["Model-1"].parts["Part-1"].datums.keys()[0]

    # Use this datum to locate the target face.
    point = mdb.models["Model-1"].parts["Part-1"].datums[datum_id].pointOn

    # Note that this is a FaceArray object i.e. a sequence of faces.
    faces = mdb.models["Model-1"].rootAssembly.instances["Part-1-1"].faces.findAt(((point, )))

    # Use the target face to build a surface.
    surface = mdb.models["Model-1"].rootAssembly.Surface(side1Faces=faces, name="Surf")

    # Create the surface traction.
    load_name = "surface_traction"
    mdb.models["Model-1"].SurfaceTraction(name=load_name, createStepName=step_name, 
                                          magnitude=float(mag), region=surface, field='',
                                          distributionType=UNIFORM, directionVector=((0, 0, 0), tuple(direction)), localCsys=None, traction=GENERAL)

    dp("Running the simulation with surface traction with direction " + str(mdb.models["Model-1"].loads[load_name].directionVector))

    # Delete all the jobs.
    dp("Deleting all jobs in the model...")
    for key in mdb.jobs.keys():
        del mdb.jobs[key]

    # Create a new job and run it.
    dp("Creating job with name " + job_name)
    mdb.Job(job_name, "Model-1", resultsFormat=BOTH)
    dp("Running job " + job_name)
    job = mdb.jobs[job_name]
    job.submit()
    job.waitForCompletion()



# Extract the displacements of some nodes from an ODB.
# If the nodes numbers are specified, use those. Otherwise choose some nodes
#    randomly.
def extract_nodal_displacements(name, step_name, cnt, nodes=None):

    # Open the ODB.
    odb = odbAccess.openOdb(name + ".odb")

    """
    dp("The names of the parts in the ODB are " + str(odb.parts.keys()))
    dp("The number of nodes in the part in the odb is " + str(len(odb.parts["PART-1"].nodes)))
    dp("The number of nodes in the assembly in the odb is " + str(len(odb.rootAssembly.nodes)))
    dp("The names of the steps are " + str(odb.steps.keys()))
    dp("There are " + str(len(odb.steps[force_app_step_name].frames)) + " frames")
    dp("The names of the field outputs for the last frame are " + str(odb.steps[force_app_step_name].frames[-1].fieldOutputs.keys()))
    disp_field_output = odb.steps[force_app_step_name].frames[-1].fieldOutputs["U"]
    dp("For the displacement field output object, the type is " + str(disp_field_output.type) + " and the description is " + disp_field_output.description)
    dp("The length of the FieldLocationArray in the displacement FieldOutput object is " + str(len(disp_field_output.locations)))
    dp("The length of the values in the displacement FieldOutput object is " + str(len(disp_field_output.values)))
    for val in disp_field_output.values:
        dp("The position of this FieldValue object is " + str(val.position))
        dp("The node label of this FieldValue object is " + str(val.nodeLabel))
        dp("The data associated with this value is " + str(val.data))
    init_step_disp_field_output = odb.steps["Initial"].frames[-1].fieldOutputs["U"]
    dp("For the initial step displacement field output object, the type is " + str(init_step_disp_field_output.type) + " and the description is " + init_step_disp_field_output.description)
    """

    disp_output = odb.steps[step_name].frames[-1].fieldOutputs["U"]
    if nodes is None:
        # If no nodes are specified, pick some randomly.
        total_node_cnt = len(disp_output.values)
        dp("There were " + str(total_node_cnt) + " nodes detected")
        nodes = random.sample(range(total_node_cnt), cnt)
        dp("The nodes randomly selected are " + str(nodes))
    else:
        # The nodes are specified.
        dp("The nodes being used are " + str(nodes))

    # Record the displacements of the nodes.
    ret = {}
    for val in disp_output.values:
        if val.nodeLabel in nodes:
            dp("The displacement of node " + str(val.nodeLabel) + " is " + str(val.data)) 
            displacement = np.array(val.data)
            ret[val.nodeLabel] = displacement 

    return ret



if __name__ == "__main__":

    mdb = openMdb("testing_linear_relationships.cae")

    # Make sure that the step exists.
    force_app_step_name = "Apply_Force"
    assert(force_app_step_name in mdb.models["Model-1"].steps.keys())

    # Run #1.
    direction1 = np.array([0, -1, 0])
    mag1 = 500000
    job1 = "job1"
    apply_force(mdb, mag1, direction1, force_app_step_name, job1)
    disps1 = extract_nodal_displacements(job1, force_app_step_name, 5)
    nodes = disps1.keys()

    # Run #2.
    direction = np.array([1, -1, 0])
    direction2 = (1 / linalg.norm(direction)) * direction
    mag2 = 5000000
    job2 = "job2"
    apply_force(mdb, mag2, direction2, force_app_step_name, job2)
    disps2 = extract_nodal_displacements(job2, force_app_step_name, 0, nodes=nodes)

    # Run #3.
    direction = np.array([1, -1, 1])
    direction3 = (1 / linalg.norm(direction)) * direction
    mag3 = 2500000
    job3 = "job3"
    apply_force(mdb, mag3, direction3, force_app_step_name, job3)
    disps3 = extract_nodal_displacements(job3, force_app_step_name, 0, nodes=nodes)

    # Run #4.
    direction = np.array([.333, 0, -1])
    direction4 = (1 / linalg.norm(direction)) * direction
    mag4 = 1230000
    job4 = "job4"
    apply_force(mdb, mag4, direction4, force_app_step_name, job4)
    disps4 = extract_nodal_displacements(job4, force_app_step_name, 0, nodes=nodes)


    # *****************************************************
    #               Computing Linear Relationship
    # *****************************************************

    # The applied force vectors.
    fv1 = mag1 * direction1
    fv2 = mag2 * direction2
    fv3 = mag3 * direction3
    fv4 = mag4 * direction4

    # Build the force matrix.
    force = np.stack((fv1, fv2, fv3), axis=0).T
    dp("The force matrix is: ")
    dp(str(force))

    # Build the displacement matrix for a particular node!
    node = nodes[0]
    dp("Computing the linear relationship for node " + str(node) + "!")
    displacement = np.stack((disps1[node], disps2[node], disps3[node]), axis=0).T
    dp("The displacement matrix is: ")
    dp(str(displacement))

    # Solve for the linear relationship for this node!
    res = np.dot(force, linalg.inv(displacement))
    dp("The linear relationship for this node is described by: ")
    dp(str(res))

    # Use the linear relationship for this node to predict the displacement 
    #    for the traction which was not used to compute it. 
    actual_disp = disps4[node]
    dp("The actual displacement of node " + str(node) + " is " + str(actual_disp))
    predicted_displacement = np.dot(linalg.inv(res), fv4)
    dp("The predicted displacement is " + str(predicted_displacement))
