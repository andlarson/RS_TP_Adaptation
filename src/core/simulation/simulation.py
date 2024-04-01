"""
Provides top-level functionality to do simulations of interest in Abaqus.
"""

from __future__ import annotations
import os
from typing import Optional, Any, TYPE_CHECKING

from abaqus import *
from abaqusConstants import * 
import numpy as np

# Imports only used for static analysis / type checking.
# This prevents cyclic imports, which are sometimes necessary for type annotations,
#     from causing runtime errors. 
if TYPE_CHECKING:
    import src.core.machining.machining as mach
    import src.core.tool_pass.tool_pass as tp
    import src.core.metadata.metadata as md
    import src.core.residual_stress.residual_stress as rs

import src.core.abaqus.abaqus_shim as shim
import src.core.boundary_conditions.boundary_conditions as bc
import src.core.metadata.naming as naming
import src.core.metadata.abaqus_metadata as abq_md
import src.util.geom as geom

from src.util.debug import *        



def sim_tool_pass_plan(tool_pass_plan: tp.ToolPassPlan, name: str, target_dir: str, 
                       machining_invariants: mach._MachiningInvariants,
                       mdb_md: abq_md.AbaqusMdbMetadata, mdb: Any,
                       stress_state: str) -> None:
    """Simulates consecutive tool passes and saves off the results.
    
       Args:
           tool_pass_plan:      The tool pass plan to simulate. 
           name:                Desired name of the MDB and the files which result
                                    from the tool path simulations. 
           target_dir:          Absolute path to directory. This directory should
                                    have a subdirectory with name name created in
                                    it already.
           machining_invariants: Invariants for this machining process.
           mdb_md:               Metadata for the MDB.
           mdb:                  Abaqus MDB object. The MDB in which to do the
                                     simulations. Note that this MDB is saved as.
                                     This MDB need not be empty, but the tool passes
                                     are simulated on top of whatever content it
                                     has.
           stress_state:         Either an absolute path to an object file containing
                                     a stress subroutine or an absolute path to a .sim
                                     file containing the stress state which resulted
                                     from a simulation. 
       
       Returns:
           None.

       Raises:
           None.
    """
    
    # Change the CWD to the directory which should already exist for simulation
    #     results.
    cwd = os.getcwd()
    sim_result_dir = os.path.join(target_dir, name)
    os.chdir(sim_result_dir)
    
    for tool_pass in tool_pass_plan.plan:
        _sim_single_tool_pass(name, tool_pass, machining_invariants, mdb_md, mdb, stress_state)    
        tool_pass_plan.pop()
    
    # Don't permanently change the CWD.
    os.chdir(cwd)

    # Save the MDB which has all the simulated tool passes in it.
    new_mdb_path = os.path.join(sim_result_dir, name)
    shim.save_mdb_as(new_mdb_path, mdb)



def _sim_single_tool_pass(job_name: str, tool_pass: tp.ToolPass, 
                          machining_invariants: mach._MachiningInvariants,
                          mdb_md: abq_md.AbaqusMdbMetadata, mdb: Any,
                          stress_state: str) -> None:
    """Simulates a single tool pass. This tool pass may be the first tool pass
           in a tool pass plan, or the nth tool pass in a tool pass plan.

       Args:
           job_name:             Desired name of the job.
           tool_pass:            The tool path to simulate.
           machining_invariants: Invariants for this machining process.
           mdb_md:               Metadata for the MDB.
           mdb:                  Abaqus MDB object. The MDB in which the tool pass
                                     will be simulated. The MDB may already contain
                                     simulated tool passes.
           stress_state:         Either an absolute path to an object file containing
                                     a stress subroutine or an absolute path to a .sim
                                     file containing the stress state which resulted
                                     from a simulation. 

       Returns:
           None.

       Raises:
           AssertionError: The MDB is not in a recognized state.
    """
    
    # The way that the MDB is simulated depends on its content.
    if shim.check_simple_standard_mdb(False, mdb):
        _sim_first_tool_pass(job_name, tool_pass, machining_invariants, mdb_md, mdb, stress_state)
    elif shim.check_job_submissions(False, mdb):
        _sim_nth_tool_pass(job_name, tool_pass, machining_invariants, mdb_md, mdb)
    else:
        assert False, "Can't figure out how to do next tool pass..."



def _sim_first_tool_pass(job_name: str, tool_pass: tp.ToolPass, 
                         machining_invariants: mach._MachiningInvariants,
                         mdb_md: abq_md.AbaqusMdbMetadata, mdb: Any, 
                         stress_state: str) -> None:
    """Simulates the first tool pass in an MDB. The MDB is assumed to not contain
           any tool pass simulations.
       
       Args:
           job_name:             The desired name of the job.
           tool_pass:            The tool path to simulate.
           machining_invariants: Invariants for this machining process.
           mdb:                  Abaqus MDB object. The MDB in which the tool pass
                                     will be simulated. 
           mdb_md:               The metadata for the MDB.
           stress_state:         Either an absolute path to an object file containing
                                     a stress subroutine or an absolute path to a .sim
                                     file containing the stress state which resulted
                                     from a simulation. 

       Returns:
           None.

       Raises:
           RuntimeError: The MDB is not in the expected state (e.g. it does not
                             have the expected content).
    """

    if not shim.check_simple_standard_mdb(False, mdb):
        raise RuntimeError("The MDB is not in the expected state.")

    names = naming.ModelNames(naming.ModelTypes.FIRST_TOOL_PASS, mdb_md)

    shim.create_material(machining_invariants.material, names.new_model_name, mdb)
    shim.create_section(names.new_model_name, mdb)
    shim.assign_section_to_whole_part(names.pre_tool_pass_part_name, names.new_model_name, mdb)
    
    # Since this modified the input file directly, this should be the last thing
    #     which happens before job creation and submission.
    _do_boilerplate_tp_sim_ops(tool_pass, stress_state, names, mdb_md, machining_invariants, mdb)

    job = shim.create_job(job_name, names.new_model_name, mdb_md, mdb) 

    if stress_state.endswith(".o"):
        # Associate the user subroutine with the job.
        shim.add_user_subroutine(job, stress_state)

    shim.run_job(job)



def _sim_nth_tool_pass(job_name: str, tool_pass: tp.ToolPass, 
                       machining_invariants: mach._MachiningInvariants, 
                       mdb_md: abq_md.AbaqusMdbMetadata, mdb: Any) -> None:
    """Simulates the nth tool pass in an MDB. The MDB is assumed to already
           contain at least one tool pass simulation.

       This function does not accept a stress subroutine because it uses the
           stress state which resulted from the last tool pass simulation as
           the initial stress state for the current simulation.

       Args:
           job_name:             Desired name of the job.
           tool_pass:            The tool pass to simulate.
           machining_invariants: Invariants associated with this machining process.
           mdb_md:               The metadata for the MDB.
           mdb:                  Abaqus MDB object. The MDB in which the tool pass
                                     will be simulated. 

       Returns:
           None.

       Raises:
           None.
    """

    names = naming.ModelNames(naming.ModelTypes.NTH_TOOL_PASS, mdb_md)
    last_odb_file_name = naming.last_odb_file_name(mdb_md)

    # The .sim file contains the stress information which resulted from the last
    #     simulation.
    # Using the name instead of the absolute path assumes something about the
    #     CWD. 
    last_sim_file_name = naming.last_sim_file_name(mdb_md)

    # Create the new model with the deformed part in it from the ODB. 
    shim.create_model_and_part_from_mesh(last_odb_file_name, names.new_model_name, mdb_md, mdb)

    # Do material and section creation.
    shim.create_material(machining_invariants.material, names.new_model_name, mdb)
    shim.create_section(names.new_model_name, mdb)

    # Do the section assignment only after the orphan mesh has been mapped to a
    #    geometry. This cannot be done before mapping to a geometry. 
    shim.assign_section_to_whole_part(names.pre_tool_pass_part_name, names.new_model_name, mdb)

    _do_boilerplate_tp_sim_ops(tool_pass, last_sim_file_name, names, mdb_md, machining_invariants, mdb)

    # Map the stress profile which existed at the end of the last simulation
    #    onto the new geometry. 
    # Since this modifies the input file directly, this should be the last
    #    thing that happens before the job is submitted and runs.
    shim.inp_map_stress(last_sim_file_name, names.new_model_name, mdb)

    # Create the job and submit it.
    job = shim.create_job(job_name, names.new_model_name, mdb_md, mdb) 
    shim.run_job(job)



def _do_boilerplate_tp_sim_ops(tool_pass: tp.ToolPass, 
                               stress_state: str,
                               names: naming.ModelNames, 
                               mdb_metadata: abq_md.AbaqusMdbMetadata,
                               machining_invariants: mach._MachiningInvariants,
                               mdb: Any
                              ) -> None:
    """Does operations which need to be done for just about any tool pass
           simulation.
       
       This modifies the input file directly, so it should be the last thing
           which is done before a job is created and submitted.
       
       Args:
           tool_pass:            The tool pass to simulate.
           stress_state:         Either an absolute path to an object file containing
                                     a stress subroutine or an absolute path to a .sim
                                     file containing the stress state which resulted
                                     from a simulation. 
           names:                The names associated with the model which contains
                                     the tool pass material removal in the MDB.
           mdb_metadata:         Metadata associated with the the MDB.
           machining_invariants: Invariants for this machining process.
           mdb:                  Abaqus MDB object. The MDB in which the tool pass
                                     will be simulated. 

       Returns:
           None.

       Raises:
           None.
    """
    
    post_tool_pass_instance = _do_cut(tool_pass, names, mdb_metadata, mdb)

    bc.apply_BCs(machining_invariants.boundary_conditions, shim.STANDARD_INITIAL_STEP_NAME, post_tool_pass_instance, names.new_model_name, mdb)

    shim.mesh(post_tool_pass_instance, 20, names.new_model_name, mdb)

    last_step_name = mdb_metadata.models_metadata[names.new_model_name].step_names[-1]

    if stress_state.endswith(".o"):
        # A stress state in an object file may not achieve mechanical equilibrium.
        shim.create_step(names.equilibrium_step_name, last_step_name, names.new_model_name, mdb_metadata, mdb)
        shim.create_step(names.deformation_step_name, names.equilibrium_step_name, names.new_model_name, mdb_metadata, mdb)
        shim.inp_add_stress_subroutine(names.new_model_name, mdb)
    elif stress_state.endswith(".sim"):
        # A stress state in an .sim file was produced by the Abaqus kernel - 
        #     it is guaranteed to satisfy mechanical equilibrium.
        shim.create_step(names.deformation_step_name, last_step_name, names.new_model_name, mdb_metadata, mdb)
        shim.inp_map_stress(stress_state, names.new_model_name, mdb)
    else:
        assert False, "Unknown file type."



def _do_cut(tool_pass: tp.ToolPass, 
            names: naming.ModelNames, 
            mdb_metadata: abq_md.AbaqusMdbMetadata,
            mdb: Any
           ) -> Any:
    """Cuts out the tool pass geometry from a part.

       Deletes the part instances used to do the cutting.

       Args:
           tool_pass:           The tool pass to simulate.
           names:               The names associated with the model which contains
                                    the tool pass material removal in the MDB.
           mdb_metadata:        Metadata associated with the the MDB.
           mdb:                 Abaqus MDB object. The MDB in which the tool pass
                                    will be simulated. 
    
       Returns:
           Abaqus PartInstance object. An instance of the part which has material
               removed in the region of the tool pass.
    
       Raises:
           None.
    """

    model = mdb.models[names.new_model_name]
    initial_geom_part = model.parts[names.pre_tool_pass_part_name]
    initial_geom_instance = shim.instance_part_into_assembly(names.pre_tool_pass_part_name, initial_geom_part, False, names.new_model_name, mdb)  

    tool_pass_part = shim.build_tool_pass_part(names.tool_pass_part_name, tool_pass, names.new_model_name, mdb_metadata, mdb)
    tool_pass_part_instance = shim.instance_part_into_assembly(names.tool_pass_part_name, tool_pass_part, False, names.new_model_name, mdb)

    cutting_instances = (tool_pass_part_instance, )
    post_tool_pass_part = shim.cut_instances_in_assembly(names.post_tool_pass_part_name, initial_geom_instance, cutting_instances, names.new_model_name, mdb_metadata, mdb)

    shim.delete_assembly_feature(names.pre_tool_pass_part_name, model)    
    shim.delete_assembly_feature(names.tool_pass_part_name, model)

    shim.assign_only_section_to_part(post_tool_pass_part, names.new_model_name, mdb)

    # Simplify the part geometry using some heuristics.
    # This cannot be done earlier because parts with virtual topologies cannot be used in boolean operations...
    shim.add_virtual_topology(post_tool_pass_part)

    # Instance the post cut geometry.
    # The instance is independent so that meshing can be done in the Assembly
    #    module.
    post_tool_pass_instance = shim.instance_part_into_assembly(names.post_tool_pass_part_name, post_tool_pass_part, False, names.new_model_name, mdb)

    return post_tool_pass_instance



def _simulate_traction_app(traction: geom.Vec3D,
                           point_near_face: geom.Point3D, 
                           model_to_copy: str,
                           machining_invariants: mach._MachiningInvariants,
                           mdb_md: abq_md.AbaqusMdbMetadata,
                           mdb: Any) -> str:
    """Simulates applying a traction to a face.

       Args:
           traction:             The traction vector to apply.
           point_near_face:      A point which lies on, or very near, the face on which
                                     to apply the traction.
           model_to_copy:        The name of the model to copy. Assumed that this
                                    model contains a single part instance. This part
                                    instance should have the face the point is near.
                                    The face on this part instance is the face to
                                    which the traction is applied.
           machining_invariants: Invariants for this machining process.
           mdb_md:               Metadata associated with the MDB. 
           mdb:                  Abaqus MDB object. The MDB to use.
    
       Returns:
           The name of the job which was run. This is not an absolute path.
    
       Raises:
           None.
    """

    names = naming.ModelNames(naming.ModelTypes.TRACTION_APP, mdb_md)
    
    shim.copy_model(model_to_copy, names.new_model_name, mdb_md, mdb)
    
    instance_name = shim.get_only_instance_name(mdb.models[names.new_model_name])
    instance = mdb.models[names.new_model_name].rootAssembly.instances[instance_name]

    bc.apply_BCs(machining_invariants.boundary_conditions, 
                 shim.STANDARD_INITIAL_STEP_NAME, 
                 instance, names.new_model_name, mdb)
    
    # Add a step for traction application. 
    model_metadata = mdb_md.models_metadata[names.new_model_name]
    last_step_name = model_metadata.step_names[-1]
    shim.create_step(names.traction_step_name, last_step_name, names.new_model_name, 
                     mdb_md, mdb)
    
    face = shim.get_closest_face([point_near_face], instance, True)

    shim.create_traction(shim.STANDARD_SURFACE_TRACTION_NAME, names.traction_step_name,
                         traction, face, names.new_model_name, mdb) 

    shim.mesh(instance, 20, names.new_model_name, mdb)

    job = shim.create_job(names.traction_job_name, names.new_model_name, mdb_md, mdb)

    shim.run_job(job)

    job_name = mdb_md.models_metadata[names.new_model_name].job_name

    return job_name 



def estimate_residual_stresses(deformed_geometry: str, target_geometry: str,
                               tool_pass: tp.ToolPass, 
                               machining_invariants: mach._MachiningInvariants,
                               path: str
                              ) -> rs.ConstantResidualStressField:
    """Estimates the residual stresses which existed in a region of material
           which was removed and caused a deformation of the workpiece.
       
       Args:
           deformed_geometry:    Absolute path to .odb or .stl file containing the
                                     deformed geoemetry from real life.
           target_geometry:      Absolute path to .odb or .stl file containing the
                                     target geometry (aka the pre-cut, pre-deform
                                     geometry) from real life.
           tool_pass:            The tool pass that removed material.
           machining_invariants: Invariants for this machining process.
           path:                 Absolute path to directory. Any MDBs which 
                                     need to be created in order to recover 
                                     the residual stresses are placed in this 
                                     directory.
    
       Returns:
           The residual stress field in the material that was removed by the
               committed tool pass plan.
    
       Raises:
           None.
    """

    mdb_md, mdb = _prepare_for_stress_estimation(deformed_geometry, target_geometry,
                                                 tool_pass, machining_invariants, path)

    return _recover_constant_residual_stress(tool_pass, machining_invariants, path, mdb_md, mdb)



def _prepare_for_stress_estimation(deformed_geometry: str, target_geometry: str,
                                   tool_pass: tp.ToolPass, 
                                   machining_invariants: mach._MachiningInvariants,
                                   path: str
                                  ) -> tuple[abq_md.AbaqusMdbMetadata, Any]:
    """This function prepares for the stress estimation step by consolidating
           the real world data into a single MDB.

       Args:
           deformed_geometry:    Absolute path to .odb or .stl file containing the
                                     deformed geoemetry from real life.
           target_geometry:      Absolute path to .odb or .stl file containing the
                                     target geometry (aka the pre-cut, pre-deform
                                     geometry) from real life.
           tool_pass:            The tool pass that removed material.
           machining_invariants: Invariants associated with the machining process.
           path:                 Absolute path to directory. Any MDBs which 
                                     need to be created in order to recover 
                                     the residual stresses are placed in this 
                                     directory.

       Returns:
           An MDB and associated metadata. The MDB contains a model for the
               deformed geometry and a model for the post-cut, pre-deform
               geometry. Each model contains exactly one part. Additionally,
               the model which contains the post-cut, pre-deform geometry
               has a material, a section definition, and a single part instance.

       Raises:
           None.
    """

    # The .stl/.odb formats for the deformed and target geometries aren't very
    #     useful since we want to apply a traction to the post-cut, pre-deform
    #     geometry. Thus, we need to create an MDB, put the deformed geometry
    #     and target geoemtry into it, and then construct the post-cut,
    #     pre-deform geometry.

    new_mdb_name = shim.STANDARD_STRESS_ESTIMATION_MDB_NAME
    mdb = shim.create_mdb(new_mdb_name, path)
    full_path = os.path.join(path, new_mdb_name)
    mdb_md = abq_md.AbaqusMdbMetadata(full_path) 
    
    # TODO: Simplify geometry?

    names = naming.ModelNames(naming.ModelTypes.DEFORMED_GEOM, mdb_md) 
    shim.import_mesh(deformed_geometry, names.deformed_geom_part_name, names.new_model_name, mdb_md, mdb)

    names = naming.ModelNames(naming.ModelTypes.TARGET_GEOM, mdb_md)
    shim.import_mesh(target_geometry, names.target_geom_part_name, names.new_model_name, mdb_md, mdb)

    # Since tractions are going to be applied to the target geometry, it must
    #     have a material, etc.
    shim.create_material(machining_invariants.material, names.new_model_name, mdb)
    shim.create_section(names.new_model_name, mdb)
    shim.assign_section_to_whole_part(names.target_geom_part_name, names.new_model_name, mdb)

    # Construct the post-cut, pre-deform geometry.
    # Note that this also does geometry simplification!
    _do_cut(tool_pass, names, mdb_md, mdb)

    # Get rid of the model created by default.
    shim.delete_model(shim.STANDARD_MODEL_NAME, mdb_md, mdb)

    return mdb_md, mdb



def _recover_constant_residual_stress(tool_pass: tp.ToolPass,
                                      machining_invariants: mach._MachiningInvariants,
                                      path: str,
                                      mdb_md: abq_md.AbaqusMdbMetadata,
                                      mdb: Any
                                     ) -> rs.ConstantResidualStressField:
    """Recovers the residual stresses which existed in the region of material 
           which was removed by a committed tool pass. 

       This technique assumes that the residual stress field which existed in
           the region of material was constant before it was removed.
        
       Args:
           tool_pass:            The tool pass that removed material.
           machining_invariants: Invariants for this machining process.
           path:                 Absolute path to directory. Any MDBs which 
                                     need to be created in order to recover 
                                     the residual stresses are placed in this 
                                     directory.
           mdb_md:               Metadata for MDB.
           mdb:                  Abaqus MDB object. This MDB is expected to contain 
                                     two models. The first model should contain a 
                                     single part, the deformed geometry. The second 
                                     model should also contain a single part, the 
                                     post-cut, pre-deformed geometry. 

       Returns:
           The residual stress field in the material that was removed by the
               committed tool pass plan.
    
       Raises:
           None.
    """
   
    # Find a point on the trench for the tool pass.
    trench_point = _find_trench_point(tool_pass) 
    
    # Some heustics for gradient ascent.
    cur = geom.Vec3D(np.array([10**9, 0, 0]))
    iteration_cnt = 3 
    step_size = 10**8

    for _ in range(iteration_cnt):
        grad = _compute_gradient(cur, trench_point, machining_invariants, path, mdb_md, mdb)
        cur = cur + geom.Vec3D(step_size * grad.rep)

    # Map the "good traction vector" to the components of the stress tensor
    #     assuming a constant state of stress.
    
    return None


     
def _find_trench_point(toolpass: tp.ToolPass) -> geom.Point3D:
    """Finds a point on the trench of a tool pass.

       Args:
           toolpass: The tool pass for which to find the trench.
    
       Returns:
           Point which is on the trench of a tool pass.
    
       Raises:
           None.
    """
    
    # Assumes an orientation of the tool with respect to the path that the tool
    #     follows.
    median_index = int(len(toolpass.path.v_list) / 2)
    point_on_trench = toolpass.path.v_list[median_index]

    return point_on_trench 



def _compute_gradient(at: geom.Vec3D, 
                      point_near_face: geom.Point3D,
                      machining_invariants: mach._MachiningInvariants,
                      path: str,
                      mdb_md: abq_md.AbaqusMdbMetadata,
                      mdb: Any
                     ) -> geom.Vec3D:
    """Computes an estimate of the gradient of the volume-intersection function
           at a particular point. 
       
       Implementation Details:
       Computes the value of the volume-intersection function at small offsets
           (in orthogonal directions) from the point of offset. Formally, this
           technique is a finite-differences style gradient estimation. 

       Args:
           at:                   The traction vector to apply. This is the point
                                     at which the gradient is estimated. 
           point_on_face:        Point on the face on which to apply the tractions.
           machining_invariants: Invariants for this machining process.
           path:                 Absolute path to a directory where simulation
                                     results will be dumped.
           mdb_md:               Metadata associated with the MDB.
           mdb:                  Abaqus MDB object. This MDB is expected to contain 
                                     at least two models. The first model should 
                                     contain a single part, the deformed geometry. 
                                     The second model should also contain a single 
                                     part, the post-cut, pre-deformed geometry. 
    
       Returns:
           An estimate of the gradient.
    
       Raises:
           None.
    """
    
    # Offset relative to length. 
    length = at.len()
    offset = length * .05

    offset_in_x = at + geom.Vec3D(np.array([offset, 0, 0]))
    offset_in_y = at + geom.Vec3D(np.array([0, offset, 0]))
    offset_in_z = at + geom.Vec3D(np.array([0, 0, offset]))
    
    vol_diff = _evaluate_vol_intersec(at, point_near_face, machining_invariants, path, mdb_md, mdb)
    vol_diff_x = _evaluate_vol_intersec(offset_in_x, point_near_face, machining_invariants, path, mdb_md, mdb)
    vol_diff_y = _evaluate_vol_intersec(offset_in_y, point_near_face, machining_invariants, path, mdb_md, mdb)
    vol_diff_z = _evaluate_vol_intersec(offset_in_z, point_near_face, machining_invariants, path, mdb_md, mdb)
    
    # Moving in the direction of the slope is equivalent to moving in the direction
    #     of greatest ascent.
    x_slope = (vol_diff_x - vol_diff) / offset
    y_slope = (vol_diff_y - vol_diff) / offset
    z_slope = (vol_diff_z - vol_diff) / offset
    
    return geom.Vec3D(np.array([x_slope, y_slope, z_slope]))
    


def _evaluate_vol_intersec(at: geom.Vec3D, 
                           point_near_face: geom.Point3D,
                           machining_invariants: mach._MachiningInvariants,
                           path: str,
                           mdb_md: abq_md.AbaqusMdbMetadata,
                           mdb: Any
                          ) -> int :
    """Evaluates the volume intersection function at some point. 

       The volume-intersection function is a function with domain which is the
           set of all possible traction vectors and range which is all positive 
           reals. It is defined for a particular target geometry, a particular
           deformed geometry, and a face. To evaluate the function at a point, 
           you simply apply the traction (the point) to the target geometry on 
           the specified face, overlap the resulting geometry in space with the 
           deformed geometry, and compute the volume of the pieces which are
           in the intersection of the volumes.

       Implementation Details:
       By applying traction to the target geometry, and not the deformed geometry,
           there is no uncertainty in the location of the trench (which may be
           in a different location due to deformation).
                   
       Args:
           at:                   The traction vector to apply. This is the point
                                     at which the function is evaluated. 
           point_on_face:        Point on the face on which to apply the tractions.
           machining_invariants: Invariants for this machining process.
           path:                 Absolute path to a directory where simulation
                                     results will be dumped.
           mdb_md:               Metadata associated with the MDB.
           mdb:                  Abaqus MDB object. This MDB is expected to contain 
                                     at least two models. The first model should 
                                     contain a single part, the deformed geometry. 
                                     The second model should also contain a single 
                                     part, the post-cut, pre-deformed geometry.
    
       Returns:
           The value of the volume intersection function (a scalar representing the
               volume) when evaluated at the point of interest.
    
       Raises:
           None.
    """
    
    # Step 1:
    # Simulate the traction application.

    names = naming.ModelNames(naming.ModelTypes.INTERSECTION, mdb_md)

    # Change the CWD to fix location of simulation results.
    orig_cwd = os.getcwd()
    os.chdir(path)
    job_name = _simulate_traction_app(at, point_near_face, names.target_geom_model_name,
                                      machining_invariants, mdb_md, mdb)
    path_to_odb = os.path.join(path, job_name) + ".odb"
    os.chdir(orig_cwd)

    # Step 2:
    # Create the model in the MDB for computing symmetric difference. 

    # Create the model by importing the result of traction application.
    shim.import_mesh(path_to_odb, names.post_traction_target_geom_part_name, names.new_model_name, mdb_md, mdb)
    
    # Copy the deformed geometry part into the model.
    shim.copy_part(names.deformed_geom_part_name, names.existing_deformed_geom_part_name,
                   names.new_model_name, names.deformed_geom_model_name,
                   mdb_md, mdb) 

    # Step 3:
    # Compute the volume of the intersection.
    volume = _compute_vol_intersection(names, mdb_md, mdb)

    return volume 



def _compute_vol_intersection(names: naming.ModelNames, 
                              mdb_md: abq_md.AbaqusMdbMetadata, 
                              mdb: Any) -> int:
    """Computes the volume of the intersection between two parts. 

       Args:
           names:  Names for a intersection operation. 
           mdb_md: Metadata associated with the MDB.
           mdb:    Abaqus MDB object. Assumed to contain a model ready for the
                       intersection operation. 
    
       Returns:
           Volume of the symmetric difference.
    
       Raises:
           None.
    """

    model_name = names.new_model_name
    part1_name = names.post_traction_target_geom_part_name
    part2_name = names.deformed_geom_model_name
    intersection_part_name = names.intersection_part_name
    part1_remainder_name = names.remainder_part_name

    # Step 1: 
    # Instance the parts into the assembly.
    part1 = mdb.models[model_name].parts[part1_name]
    part1_instance = shim.instance_part_into_assembly(part1_name, part1, False, model_name, mdb)

    part2 = mdb.models[model_name].parts[part2_name]
    part2_instance = shim.instance_part_into_assembly(part2_name, part2, False, model_name, mdb)

    # Step 2:
    # Find the intersection of the two parts.
    # Cut, out of first part, the overlapping portion of the second part. Some
    #     portion of the first part remains. Cut this remaining portion out
    #     of the first part. This yields the intersection of the two parts.

    # Cut out the overlapping portion of the second part.
    part1_remainder = shim.cut_instances_in_assembly(part1_remainder_name, part1_instance,
                                                     (part2_instance,), model_name, mdb_md,
                                                     mdb)
    
    remainder_part_instance = shim.instance_part_into_assembly(part1_remainder_name, 
                                                               part1_remainder, False, 
                                                               model_name, mdb)
    
    # Unsuppress the instance of part1 so it can be used in the subsequent cutting
    #     operation.
    shim.resume_assembly_feature(part1_name, mdb.models[model_name])
    part1_instance = mdb.models[model_name].rootAssembly.instances[part1_name]

    # Do the cut to get the intersection.
    intersection = shim.cut_instances_in_assembly(intersection_part_name, part1_instance, 
                                                  (remainder_part_instance,), model_name,
                                                  mdb_md, mdb)

    # Step 3:
    # Compute the volume of the intersection.
    return shim.compute_part_volume(intersection)



# *****************************************************************************
#                                 DEPRECATED 
# *****************************************************************************


def _compute_symmetric_diff(names: naming.ModelNames, 
                            mdb_md: abq_md.AbaqusMdbMetadata, 
                            mdb: Any) -> int:
    """DEBRECATED. Replaced by computing the intserction.
    
       Computes the symmetric difference in volume between two part geometries. 

       In set theory, the symmetric difference between two sets A, B is the set
           of elements which are in either A, B and are not in the intersection
           of A, B. This function computes the volume of the symmetric difference. 
                   
       Args:
           names:  Names for a symmetric difference operation. 
           mdb_md: Metadata associated with the MDB.
           mdb:    Abaqus MDB object. Assumed to contain a model ready for the
                       symmetric difference operation.
    
       Returns:
           Volume of the symmetric difference.
    
       Raises:
           None.
    """

    model_name = names.new_model_name
    part1_name = names.post_traction_target_geom_part_name
    part2_name = names.deformed_geom_model_name
    part1_remainder_name = names.part1_remainder_part_name
    intersection_part_name = names.intersection_part_name
    union_part_name = names.union_part_name
    symm_diff_part_name = names.symm_diff_part_name

    # Step 1: 
    # Instance the parts into the assembly.
    part1 = mdb.models[model_name].parts[part1_name]
    part1_instance = shim.instance_part_into_assembly(part1_name, part1, False, model_name, mdb)

    part2 = mdb.models[model_name].parts[part2_name]
    part2_instance = shim.instance_part_into_assembly(part2_name, part2, False, model_name, mdb)

    # Step 2:
    # Find the intersection of the two parts.
    # Cut, out of first part, the overlapping portion of the second part. Some
    #     portion of the first part remains. Cut this remaining portion out
    #     of the first part. This yields the intersection of the two parts.

    # Cut out the overlapping portion of the second part.
    part1_remainder = shim.cut_instances_in_assembly(part1_remainder_name, part1_instance,
                                                     (part2_instance,), model_name, mdb_md,
                                                     mdb)
    
    remainder_part_instance = shim.instance_part_into_assembly(part1_remainder_name, 
                                                               part1_remainder, False, 
                                                               model_name, mdb)
    
    # Unsuppress the instance of part1 so it can be used in the subsequent cutting
    #     operation.
    shim.resume_assembly_feature(part1_name, mdb.models[model_name])
    part1_instance = mdb.models[model_name].rootAssembly.instances[part1_name]

    # Do the cut to get the intersection.
    intersection = shim.cut_instances_in_assembly(intersection_part_name, part1_instance, 
                                                  (remainder_part_instance,), model_name,
                                                  mdb_md, mdb)

    # Step 3:
    # Find the symmetric difference by using the intersection.
    intersection_instance = shim.instance_part_into_assembly(intersection_part_name, 
                                                             intersection, False, 
                                                             model_name, mdb)

    # Resume the two part instances and merge them. The resulting part is the
    #     union.
    shim.resume_assembly_feature(part1_name, mdb.models[model_name])
    part1_instance = mdb.models[model_name].rootAssembly.instances[part1_name]
    shim.resume_assembly_feature(part2_name, mdb.models[model_name])
    part1_instance = mdb.models[model_name].rootAssembly.instances[part2_name]

    union = shim.merge_instances_in_assembly(union_part_name, (part1_instance, part2_instance), 
                                             False, model_name, mdb_md, mdb)
    union_instance = shim.instance_part_into_assembly(union_part_name, union, False, model_name, mdb)
    
    # Remove the intersection from the union to yield the symmetric difference.
    symmetric_difference = shim.cut_instances_in_assembly(symm_diff_part_name, union_instance,
                                                          (intersection_instance,), model_name,
                                                          mdb_md, mdb)

    # Step 4:
    # Compute the volume of the symmetric difference. 
    return shim.compute_part_volume(symmetric_difference)



def _recover_linear_relationships(points: list[geom.Point3D], 
                                  face_point: geom.Point3D,
                                  mdb_metadata: abq_md.AbaqusMdbMetadata,
                                  commit_phase_md: md.CommitmentPhaseMetadata, 
                                  mdb: Any,
                                  path: str
                                 ) -> list[Any]:
    """DEPRECATED. Was used to recover linear relationships between applied
          tractions and nodal displacements in an effort to find good tractions. 
    
       Applies test tractions to some face of a workpiece to recover the linear 
           relationships which exist between applied traction and the displacement 
           of some points.

       Args:
           points:          The points at which to recover the linear relationships.
           face_point:      A point which is on the face on which tractions will
                                be applied. You might ask: Why not just pass in
                                the face itself? The reason is that an Abaqus Face
                                object is associated with a particular part/assembly
                                in a model. Since this function creates a model
                                for each traction which is applied, the Abaqus Face
                                object would be invalid for these new models.
           mdb_metadata:    Metadata associated with the MDB.
           commit_phase_md: Metadata associated with the commitment phase.
           mdb:             Abaqus MDB object. The MDB containing the part which has the
                                face. This MDB should contain only a basic part geometry
                                (i.e. it should contain a single model, a single part,
                                no non-default steps, etc.)
           path:            Absolute path to the directory where simulation results 
                                (.odb, etc.) will be dumped.
    
       Returns:
           List of 3 x 3 numpy arrays which describe the linear relationships
               between the applied tractions and the displacements of the
               points. 
           For each matrix M, the matrix M relates the quantities: M * t = d. 
               In this equation, t is a traction vector and d is a displacement vector.
           The number of matrices equals the number of points passed to this
               function.
           Note that these linear relationships only hold when tractions are
               applied to the face passed to this function.
    
       Raises:
           None.
    """

    if not shim.check_simple_standard_mdb(True, mdb):
        raise RuntimeError("The MDB should contain a basic geometry.")

    # Set the CWD to dictate location of simulation results.
    orig_cwd = os.getcwd()
    os.chdir(path)

    x = geom.Vec3D(np.array([1, 0, 0]))
    y = geom.Vec3D(np.array([0, 1, 0]))
    z = geom.Vec3D(np.array([0, 0, 1]))
    traction_vecs = (x, y, z)
   
    job_paths = []
    for vec in traction_vecs:
        job_name = _simulate_traction_app(vec, face_point, shim.STANDARD_MODEL_NAME, mdb_metadata, commit_phase_md, mdb)
        job_paths.append(os.path.join(path, job_name))

    os.chdir(orig_cwd)
    
    displacements = []
    for p in job_paths:
        odb_path = p + ".odb"
        points_to_nodes = shim.read_displacements(points, odb_path)
        displacements.append(points_to_nodes)

    # Useful for understanding difference between points passed in and the nodes
    #     at which the linear relationships are computed.
    for d in displacements:
        dp("For a dictionary of displacements: ")
        for key in d:
            dp("At point " + str(key) + " the closest node is " + str(d[key][0]) + " and the displacement of this node is " + str(d[key][1])) 
    
    # And compute the linear relationship for each point. 
    linear_relationships = []
    for point in points: 
        # Get the displacement of the point for each traction vector application. 
        disps = []
        for d in displacements:
            disps.append(d[point][1])

        # Use the displacements to compute the linear relationship.
        M = _compute_linear_relationship(traction_vecs, tuple(disps))
        linear_relationships.append(M)

    return linear_relationships 


    
def _compute_linear_relationship(tractions: tuple[geom.Vec3D, geom.Vec3D, geom.Vec3D],
                                 displacements: tuple[geom.Vec3D, geom.Vec3D, geom.Vec3D]
                                ) -> Any:
    """DEPRECATED. Was used to compute the linear relationships between a traction
           applied on a surface and the displacement of some points.
           
       Computes the linear relationship which relates an applied traction to
           the displacement of a point.
        
       Args:
           tractions:     Three linearly independent traction vectors.
           displacements: Three displacments of a single point. The displacements
                              produced by the applied tractions.
    
       Returns:
           Numpy Array object of size 3x3. If this matrix is called M, then it
               can be used to do M * t = d, where t is a traction vector and
               d is a displacement vector.
    
       Raises:
           None.
    """

    # Normalize the tractions.
    tractions = (tractions[0].normalize(), tractions[1].normalize(), tractions[2].normalize())

    # Check that the tractions are linearly independent.
    tractions = np.stack((tractions[0].components(), tractions[1].components(), tractions[2].components()), axis=-1)
    try:
        np.linalg.inv(tractions)
    except BaseException as e:
        raise RuntimeError("Tractions are linearly dependent!") 

    # Solve for the matrix.
    displacements = np.stack((displacements[0].components(), displacements[1].components(), displacements[2].components()), axis=-1)
    res = np.linalg.solve(tractions.T, displacements.T).T

    return res
