"""
Provides top-level functionality to do simulations of interest in Abaqus.
"""

import os
from typing import Optional, Any
import pathlib
import random
import shutil
import copy
import sys
import enum

from abaqus import *
from abaqusConstants import * 
import numpy as np
# The two imports below are special. They import functionality provided by
#     Abaqus plug-ins. In order for these imports to be resolved correctly,
#     it's necessary to add installation-dependent directories to the path.
PATH_STL_IMPORT = r'/opt/caen/abaqus/abaqus-2024/SIMULIA/EstProducts/2024/'\
                  r'linux_a64/code/python3.10/lib/abaqus_plugins/stlImport'
if PATH_STL_IMPORT not in sys.path:
    sys.path.append(PATH_STL_IMPORT)
PATH_STL_EXPORT = r'/opt/caen/abaqus/abaqus-2024/SIMULIA/EstProducts/2024/'\
                  r'linux_a64/code/python3.10/lib/abaqus_plugins/stlExport'
if PATH_STL_EXPORT not in sys.path:
    sys.path.append(PATH_STL_EXPORT)
import stl2inp
import stlExport_kernel

import src.core.abaqus.abaqus_shim as shim
import src.core.tool_pass.tool_pass as tp
import src.core.boundary_conditions.boundary_conditions as bc
import src.core.metadata.metadata as md
import src.core.metadata.naming as naming
import src.core.residual_stress.residual_stress as rs
import src.core.metadata.abaqus_metadata as abq_md
import src.core.real_world_data.real_world_data as rwd
import src.util.geom as geom

from src.util.debug import *        



def sim_tool_pass_plan(tool_pass_plan: tp.ToolPassPlan, name: str, target_dir: str, 
                       commit_phase_md: md.CommitmentPhaseMetadata,
                       mdb_md: abq_md.AbaqusMdbMetadata, mdb: Any,
                       stress_subroutine: Optional[str] = None) -> None:
    """Simulates consecutive tool passes and saves off the results.
    
       Args:
           tool_pass_plan:    The tool pass plan to simulate. 
           name:              Desired name of the MDB and the files which result
                                  from the tool path simulations. 
           target_dir:        Absolute path to directory. This directory should
                                  have a subdirectory with name name created in
                                  it already.
           commit_phase_md:   The metadata associated with the commitment phase
                                  in which to simulate the tool pass plan.
           mdb_md:            Metadata for the MDB.
           mdb:               Abaqus MDB object. The MDB in which to do the
                                  simulations. Note that this MDB is saved as.
                                  This MDB need not be empty, but the tool passes
                                  are simulated on top of whatever content it
                                  has.
           stress_subroutine: Absolute path to stress subroutine. If this argument is not 
                                  specified, then the initial stress state of the 
                                  part will be sourced from the last simulation in 
                                  the previous commitment phase. If this argument 
                                  is specified, it is used instead.
       
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
        _sim_single_tool_pass(name, tool_pass, commit_phase_md, mdb_md, mdb, stress_subroutine)    
        tool_pass_plan.pop()
    
    # Don't permanently change the CWD.
    os.chdir(cwd)

    # Save the MDB which has all the simulated tool passes in it.
    new_mdb_path = os.path.join(sim_result_dir, name)
    shim.save_mdb_as(new_mdb_path, mdb)



def _sim_single_tool_pass(job_name: str, tool_pass: tp.ToolPass, 
                          commit_metadata: md.CommitmentPhaseMetadata, 
                          mdb_md: abq_md.AbaqusMdbMetadata, mdb: Any,
                          stress_subroutine: Optional[str] = None) -> None:
    """Simulates a single tool pass. This tool pass may be the first tool pass
           in a tool pass plan, or the nth tool pass in a tool pass plan.

       Args:
           job_name:          Desired name of the job.
           tool_pass:         The tool path to simulate.
           commit_metadata:   The record associated with the search for a good next 
                                  tool pass to do.
           mdb_md:            Metadata for the MDB.
           mdb:               Abaqus MDB object. The MDB in which the tool pass
                                  will be simulated. The MDB may already contain
                                  simulated tool passes.
           stress_subroutine: Path to stress subroutine. If this argument is not 
                                  specified, then the initial stress state of the 
                                  part will be sourced from the last simulation in 
                                  the previous commitment phase. If this argument 
                                  is specified, it is used instead.

       Returns:
           None.

       Raises:
           AssertionError: The MDB is not in a recognized state.
    """
    
    # The way that the MDB is simulated depends on its content.
    if shim.check_simple_standard_mdb(False, mdb):
        if stress_subroutine is not None:
            _sim_first_tool_pass(job_name, tool_pass, commit_metadata, mdb_md, mdb, stress_subroutine)
        else:
            _sim_first_tool_pass(job_name, tool_pass, commit_metadata, mdb_md, mdb)  
    elif shim.check_job_submissions(False, mdb):
        _sim_nth_tool_pass(job_name, tool_pass, commit_metadata, mdb_md, mdb)
    else:
        raise AssertionError("Can't figure out how to do next tool pass...")



def _sim_first_tool_pass(job_name: str, tool_pass: tp.ToolPass, 
                         commit_metadata: md.CommitmentPhaseMetadata, 
                         mdb_md: abq_md.AbaqusMdbMetadata, mdb: Any, 
                         stress_subroutine: Optional[str] = None) -> None:
    """Simulates the first tool pass in an MDB. The MDB is assumed to not contain
           any tool pass simulations.
       
       Args:
           job_name:          The desired name of the job.
           tool_pass:         The tool path to simulate.
           commit_metadata:   The record associated with the search for a good next 
                                  tool pass to do.
           mdb:               Abaqus MDB object. The MDB in which the tool pass
                                  will be simulated. 
           mdb_md:            The metadata for the MDB.
           stress_subroutine: Path to stress subroutine. If this argument is not 
                                  specified, then the initial stress state of the 
                                  part will be sourced from the last simulation in 
                                  the previous commitment phase. If this argument 
                                  is specified, it is used instead.

       Returns:
           None.

       Raises:
           RuntimeError: The MDB is not in the expected state (e.g. it does not
                             have the expected content).
    """

    if not shim.check_simple_standard_mdb(False, mdb):
        raise RuntimeError("The MDB is not in the expected state.")

    names = naming.ModelNames(naming.ModelTypes.FIRST_TOOL_PASS, mdb_md)

    shim.create_material(commit_metadata.init_part.material, names.new_model_name, mdb)
    shim.create_section(names.new_model_name, mdb)
    shim.assign_section_to_whole_part(names.pre_tool_pass_part_name, names.new_model_name, mdb)
    
    if stress_subroutine is not None:
        _do_boilerplate_tp_sim_ops(tool_pass, True, names, mdb_md, commit_metadata, mdb)
    else:
        _do_boilerplate_tp_sim_ops(tool_pass, False, names, mdb_md, commit_metadata, mdb)

    if stress_subroutine is not None:
        # Since this modifies the input file directly, this should be the last
        #    thing that happens before the job is submitted and runs.
        shim.inp_add_stress_subroutine(names.new_model_name, mdb)
    else:
        # Use the stress state from the last toolpath in the previous commit.
        path_sim_file = commit_metadata.path_last_commit_sim_file

        # Since this modifies the input file directly, this should be the last
        #    thing that happens before the job is submitted and runs.
        shim.inp_map_stress(path_sim_file, names.new_model_name, mdb)

    job = shim.create_job(job_name, names.new_model_name, mdb_md, mdb) 

    if stress_subroutine is not None:
        # Associate the user subroutine with the job.
        shim.add_user_subroutine(job, stress_subroutine)

    shim.run_job(job)



def _sim_nth_tool_pass(job_name: str, tool_pass: tp.ToolPass, 
                       commit_metadata: md.CommitmentPhaseMetadata, 
                       mdb_md: abq_md.AbaqusMdbMetadata, mdb: Any) -> None:
    """Simulates the nth tool pass in an MDB. The MDB is assumed to already
           contain at least one tool pass simulation.

       This function does not accept a stress subroutine because it uses the
           stress state which resulted from the last tool pass simulation as
           the initial stress state for the current simulation.

       Args:
           job_name:        Desired name of the job.
           tool_pass:       The tool pass to simulate.
           commit_metadata: The record associated with the search for a good next 
                                tool pass to do.
           mdb_md:          The metadata for the MDB.
           mdb:             Abaqus MDB object. The MDB in which the tool pass
                                will be simulated. 

       Returns:
           None.

       Raises:
           None.
    """

    names = naming.ModelNames(naming.ModelTypes.NTH_TOOL_PASS, mdb_md)
    last_odb_file_name = naming.last_odb_file_name(mdb_md)

    # The .sim file contains the stress information which resulted from the last
    #    simulation.
    last_sim_file_name = naming.last_sim_file_name(mdb_md)

    # Create the new model with the deformed part in it from the ODB. 
    shim.create_model_and_part_from_odb(names.pre_tool_pass_part_name, names.new_model_name, last_odb_file_name, mdb_md, mdb)

    # Do material and section creation.
    shim.create_material(commit_metadata.init_part.material, names.new_model_name, mdb)
    shim.create_section(names.new_model_name, mdb)

    _orphan_mesh_to_geometry(names.pre_tool_pass_part_name, names.new_model_name, mdb)

    # Do the section assignment only after the orphan mesh has been mapped to a
    #    geometry. This cannot be done before mapping to a geometry. 
    shim.assign_section_to_whole_part(names.pre_tool_pass_part_name, names.new_model_name, mdb)

    _do_boilerplate_tp_sim_ops(tool_pass, False, names, mdb_md, commit_metadata, mdb)

    # Map the stress profile which existed at the end of the last simulation
    #    onto the new geometry. 
    # Since this modifies the input file directly, this should be the last
    #    thing that happens before the job is submitted and runs.
    shim.inp_map_stress(last_sim_file_name, names.new_model_name, mdb)

    # Create the job and submit it.
    job = shim.create_job(job_name, names.new_model_name, mdb_md, mdb) 
    shim.run_job(job)



def _do_boilerplate_tp_sim_ops(tool_pass: tp.ToolPass, 
                               user_defined_stress: bool,
                               names: naming.ModelNames, 
                               mdb_metadata: abq_md.AbaqusMdbMetadata,
                               commit_metadata: md.CommitmentPhaseMetadata, 
                               mdb: Any
                              ) -> None:
    """Does operations which need to be done for just about any tool pass
           simulation.
       
       Args:
           tool_pass:           The tool pass to simulate.
           user_defined_stress: Flag indicating if a user defined stress profile
                                    (which may not satisfy mechanical equilibrium)
                                    is going to be used for this simulation.
           names:               The names associated with the model which contains
                                    the tool pass material removal in the MDB.
           mdb_metadata:        Metadata associated with the the MDB.
           commit_metadata:     Metadata associated with the commitment phase.
           mdb:                 Abaqus MDB object. The MDB in which the tool pass
                                    will be simulated. 

       Returns:
           None.

       Raises:
           None.
    """
    
    post_tool_pass_instance = _do_cut(tool_pass, names, mdb_metadata, mdb)

    bc.apply_BCs(commit_metadata.BCs, shim.STANDARD_INITIAL_STEP_NAME, post_tool_pass_instance, names.new_model_name, mdb)

    shim.mesh(post_tool_pass_instance, 20, names.new_model_name, mdb)

    last_step_name = mdb_metadata.models_metadata[names.new_model_name].step_names[-1]
    if user_defined_stress:
        shim.create_step(names.equilibrium_step_name, last_step_name, names.new_model_name, mdb_metadata, mdb)
        shim.create_step(names.deformation_step_name, names.equilibrium_step_name, names.new_model_name, mdb_metadata, mdb)
    else:
        shim.create_step(names.deformation_step_name, last_step_name, names.new_model_name, mdb_metadata, mdb)



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
                           mdb_metadata: abq_md.AbaqusMdbMetadata,
                           commit_phase_md: md.CommitmentPhaseMetadata,
                           mdb: Any) -> str:
    """Simulates applying a traction to a face.

       Args:
           traction:        The traction vector to apply.
           point_near_face: A point which lies on, or very near, the face on which
                                to apply the traction.
           model_to_copy:   The name of the model to copy. Assumed that this
                                model contains a single part instance. This part
                                instance should have the face the point is near.
                                The face on this part instance is the face to
                                which the traction is applied.
           mdb_metadata:    Metadata associated with the MDB. 
           commit_phase_md: Metadata associated with the commitment phase.
           mdb:             Abaqus MDB object. The MDB to use.
    
       Returns:
           The name of the job which was run. This is not an absolute path.
    
       Raises:
           None.
    """

    if not shim.check_first_model_simple_assembly(True, mdb_metadata, mdb):
        raise RuntimeError("The MDB is in an unexpected state.")

    names = naming.ModelNames(naming.ModelTypes.TRACTION_APP, mdb_metadata)
    
    shim.copy_model(model_to_copy, names.new_model_name, mdb_metadata, mdb)
    
    instance_name = shim.get_only_instance_name(mdb.models[names.new_model_name])
    instance = mdb.models[names.new_model_name].rootAssembly.instances[instance_name]

    bc.apply_BCs(commit_phase_md.BCs, shim.STANDARD_INITIAL_STEP_NAME, 
                 instance, names.new_model_name, mdb)
    
    # Add a step for traction application. 
    model_metadata = mdb_metadata.models_metadata[names.new_model_name]
    last_step_name = model_metadata.step_names[-1]
    shim.create_step(names.traction_step_name, last_step_name, names.new_model_name, 
                     mdb_metadata, mdb)
    
    face = shim.get_closest_face([point_near_face], instance, True)

    shim.create_traction(shim.STANDARD_SURFACE_TRACTION_NAME, names.traction_step_name,
                         traction, face, names.new_model_name, mdb) 

    shim.mesh(instance, 20, names.new_model_name, mdb)

    job = shim.create_job(names.traction_job_name, names.new_model_name, mdb_metadata, mdb)

    shim.run_job(job)

    job_name = mdb_metadata.models_metadata[names.new_model_name].job_name

    return job_name 



def create_mdb_from_odb(new_mdb_name: str, new_mdb_path: str, odb_path: str
                       ) -> tuple[Any, abq_md.AbaqusMdbMetadata, naming.ModelNames]:
    """Creates an MDB with a part in it by sourcing the content of an ODB.
           Saves the MDB before returning.
      
       Implementation Detail:
       The MDB is populated with a single part and that part is named in a
           generic way.
        
       Args:
           new_mdb_name: The desired name of the new MDB.
           new_mdb_path: Absolute path to directory. The MDB (.cae file) will be
                             created in this directory.
           odb_path:     Absolute path to ODB (.odb file).
     
       Returns:
           Tuple containing the new Abaqus MDB object, the metadata that accompanies 
               it, and the names associated with it.
    
       Raises:
           None.
    """

    # Create the new MDB.
    mdb = shim.create_mdb(new_mdb_name, new_mdb_path)
    full_path = os.path.join(new_mdb_path, new_mdb_name)
    mdb_metadata = abq_md.AbaqusMdbMetadata(full_path) 
    
    # Generate the names for the stuff in the MDB.
    names = naming.ModelNames(naming.ModelTypes.DEFAULT_MODEL, mdb_metadata)

    # Fetch the result of the committed tool pass and use it as the initial state
    #     of the MDB.
    shim.create_part_from_odb(names.part_from_odb_name, names.new_model_name, odb_path, mdb_metadata, mdb)
    _orphan_mesh_to_geometry(names.part_from_odb_name, names.new_model_name, mdb)
    
    # Save the MDB so it is visible in the file system.
    shim.save_mdb(mdb)

    return mdb, mdb_metadata, names



def _stl_to_mdb(stl_path: str, mdb_name: str, save_dir: str) -> Any:
    """Converts the representation of a workpiece from .stl (an open standard)
           format to .cae format (i.e. an Abaqus model database).

       Args:
           stl_path: Absolute path to .stl file to convert.
           mdb_name: Name of the MDB which will be created. 
           save_dir: Absolute path to directory in which to create the MDB.
                         Note that the MDB is slated to be saved to this directory,
                         but isn't actually saved to it.
    
       Returns:
           Abaqus MDB object with a single model and a single part. Note that
               this MDB is open and it has not been saved.
    
       Raises:
           None.
    """
    
    mdb = shim.create_mdb(mdb_name, save_dir)

    # It isn't documented how Abaqus plug-in's collect and and use information
    #     from an MDB. The STL import function, for example, does not take a
    #     MDB as one of its arguments. I hypothesize that the STL import plugin
    #     uses the current MDB which is open by default.
    
    # This function isn't documented in Abaqus' documentation. I found it by
    #     using the plugin and then looking at the .rpy file. Note that the
    #     mergeNodesTolerance is set to same value used by default in the
    #     .rpy file.
    stl2inp.STL2inp(stlfile=stl_path, modelName=shim.STANDARD_MODEL_NAME, 
                    mergeNodesTolerance=1E-06)

    return mdb 



class ModuleNames(enum.Enum):
    """Module names from which .stl files can be exported."""
    PART = "Part"
    MESH = "Mesh"



def _export_stl_from_mdb(module_name: ModuleNames, model_name: str,
                         stl_path: str, mdb: Any, part: Any = None
                        ) -> None:
    """Exports the data from a module to .stl format.

       Args:
           module_name: The name of the module from which to export .stl
                            data.
           model_name:  The name of the model containing the module to export.
           stl_path:    Absolute path to desired save location for .stl file. 
           mdb:         Abaqus MDB object which is open. 
       
       Optional Args:
           part:        Abaqus Part object. Note that there may be many parts 
                            in the Part module. In this case, the STL export 
                            plugin also uses which part is currently being displayed 
                            in the viewport to determine which part should be 
                            exported in .stl format. Thus, if an export from 
                            the Part module is requested, it's necessary to also
                            pass this argument.
    
       Returns:
           None.
    
       Raises:
           None.
    """
    



    





def _orphan_mesh_to_geometry(part_name: str, model_name: str, mdb: Any) -> None:
    """Converts an orphan mesh to a solid geometry by converting all the faces
           of the orphan mesh to regions/surfaces, and then melding all the
           regions/surfaces together to form a solid geometry.
           
       For some context, the result of a simulation is an orphan mesh (i.e. 
           just a bunch of vertices and elements connecting the vertices 
           together). An orphan mesh by itself is not very useful. For example, 
           to do boolean operations between parts, the parts need to have 
           geometries associated with them. This function adds geometric features 
           to a part, giving it a geometry.

       Args:
           part_name: The name of the part which contains the orphan mesh. This
                          part will have the solid geometry associated with it.
           model_name: The name of the model the part is in.
           mdb:        Abaqus MDB object.

       Returns:
           None.

       Raises:
           None.
    """

    part = mdb.models[model_name].parts[part_name]

    unique_elem_faces = shim.get_unique_element_faces(part)

    # For each unique face in the mesh, we know the face is on the surface if it is
    #    associated with exactly one element.
    # We need to construct a region for each face on the surface of the part. 
    # Each region is then used to build a geometric face feature associated
    #    with the part.

    cnt = 0
    for elem_face in unique_elem_faces:
        if len(shim.get_mesh_face_elements(elem_face)) == 1:

            face_reg = shim.build_region_with_elem_face(elem_face, part)

            shim.add_face_from_region(face_reg, part)

            # Clear the cache every once in a while for a potential speed up. 
            cnt += 1
            if cnt % 50 == 0:
                part.clearGeometryCache()

    # Note: The convert_shell_to_solid() function is significantly less failure 
    #     prone than add_solid_from_faces().
    shim.convert_shell_to_solid(part)
    
    # Remove any dependency on the orphan mesh.
    # Not doing this causes the orphan mesh to still appear in the Assembly
    #    module, which can make appearence confusing. 
    shim.suppress_part_feature(shim.STANDARD_ORPHAN_MESH_FEATURE_NAME, part)



def estimate_residual_stresses(deformed_geometry_data: rwd.RealWorldData | rwd.RealWorldDataFromSim,
                               target_geometry_data: rwd.RealWorldData | rwd.RealWorldDataFromSim,
                               tool_pass: tp.ToolPass,
                               commit_phase_md: md.CommitmentPhaseMetadata,
                               path: str
                              ) -> rs.ConstantResidualStressField:
    """Estimates the residual stresses which existed in a region of material
           which was removed and caused a deformation of the workpiece.
       
       Args:
           deformed_geometry_data: The real world data collected by scanning
                                       the workpiece after it underwent some
                                       deformation due to material removal.
           target_geometry_data:   The real world data collected by scanning
                                       the workpiece before it underwent some
                                       deformation due to material removal.
           tool_pass:                The tool pass that removed material.
           commit_phase_md:          The commitment phase during which the tool 
                                         pass was committed to.
           path:                     Absolute path to directory. Any MDBs which 
                                         need to be created in order to recover 
                                         the residual stresses are placed in this 
                                         directory.
    
       Returns:
           The residual stress field in the material that was removed by the
               committed tool pass plan.
    
       Raises:
           None.
    """

    return _recover_constant_residual_stress(deformed_geometry_data, target_geometry_data,
                                             tool_pass, commit_phase_md, path)



def _recover_constant_residual_stress(deformed_geometry_data: rwd.RealWorldData | rwd.RealWorldDataFromSim,
                                      target_geometry_data: rwd.RealWorldData | rwd.RealWorldDataFromSim, 
                                      tool_pass: tp.ToolPass,
                                      commit_phase_md: md.CommitmentPhaseMetadata,
                                      path: str
                                     ) -> rs.ConstantResidualStressField:
    """Recovers the residual stresses which existed in the region of material 
           which was removed by a committed tool pass. 

       This technique assumes that the residual stress field which existed in
           the region of material was constant before it was removed.
        
       Args:
           deformed_geometry_data: The real world data collected by scanning
                                       the workpiece after it underwent some
                                       deformation due to material removal.
           target_geometry_data:   The real world data collected by scanning
                                       the workpiece before it underwent some
                                       deformation due to material removal.
           tool_pass:              The tool pass that removed material.
           commit_phase_md:        The commitment phase during which the tool 
                                       pass was committed to.
           path:                   Absolute path to directory. Any MDBs which 
                                       need to be created in order to recover 
                                       the residual stresses are placed in this 
                                       directory.

       Returns:
           The residual stress field in the material that was removed by the
               committed tool pass plan.
    
       Raises:
           None.
    """
   
    # Step 1:
    # Removing the chunk of material due to the tool pass from the target geometry
    #     results in the geometry that would result if no residual stresses were
    #     present.
    # TODO: Do this in Blender. Note that result will be .stl file. 

    # Step 2:
    # Find a point on the trench for the tool pass.
    # This point will be used to locate the trench in the real-world part data. 
    tool_pass = commit_phase_md.committed_tpp[1].plan[0]
    trench_point = _find_trench_point(tool_pass) 

    # Step 3:
    # Do gradient descent on the volume difference function. 
    cur = geom.Vec3D(np.array([100000000, 0, 0]))
    iteration_cnt = 3 
    step_size = 1

    for _ in range(iteration_cnt):
        grad = _compute_gradient(cur, trench_point, deformed_geometry_data, target_geometry_data, commit_phase_md, path)
        cur = cur - geom.Vec3D(step_size * grad.rep)

    # Step 4:
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
                      deformed_geometry_data: rwd.RealWorldData | rwd.RealWorldDataFromSim,
                      target_geometry_data: rwd.RealWorldData | rwd.RealWorldDataFromSim, 
                      commit_phase_md: md.CommitmentPhaseMetadata,
                      path: str,
                     ) -> geom.Vec3D:
    """Computes an estimate of the gradient of the volume-difference function
           at a particular point. 
       
       Implementation Details:
       Computes the value of the volume-difference function at small offsets
           (in orthogonal directions) from the point of offset. Formally, this
           technique is a finite-differences style gradient estimation. 

       Args:
           at:                     The traction vector to apply. This is the point
                                       at which the gradient is estimated. 
           point_on_face:          Point on the face on which to apply the tractions.
           deformed_geometry_data: The real world data collected by scanning
                                       the workpiece after it underwent some
                                       deformation due to material removal.
           target_geometry_data:   The real world data collected by scanning
                                       the workpiece before it underwent some
                                       deformation due to material removal.
           commit_phase_md:        Metadata associated with the commitment phase
                                       to which the deformed MDB belongs. 
           path:                   Absolute path to a directory where simulation
                                       results will be dumped.
    
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
    
    vol_diff = _evaluate_vol_diff(at, point_near_face, deformed_geometry_data, target_geometry_data, commit_phase_md, path)
    vol_diff_x = _evaluate_vol_diff(offset_in_x, point_near_face, deformed_geometry_data, target_geometry_data, commit_phase_md, path)
    vol_diff_y = _evaluate_vol_diff(offset_in_y, point_near_face, deformed_geometry_data, target_geometry_data, commit_phase_md, path)
    vol_diff_z = _evaluate_vol_diff(offset_in_z, point_near_face, deformed_geometry_data, target_geometry_data, commit_phase_md, path)
    
    # Moving in the direction of the slope is equivalent to moving in the direction
    #     of greatest ascent.
    x_slope = (vol_diff_x - vol_diff) / offset
    y_slope = (vol_diff_y - vol_diff) / offset
    z_slope = (vol_diff_z - vol_diff) / offset
    
    return geom.Vec3D(np.array([x_slope, y_slope, z_slope]))
    


def _evaluate_vol_diff(at: geom.Vec3D, 
                       point_near_face: geom.Point3D,
                       deformed_geometry_data: rwd.RealWorldData | rwd.RealWorldDataFromSim,
                       target_geometry_data: rwd.RealWorldData | rwd.RealWorldDataFromSim, 
                       commit_phase_md: md.CommitmentPhaseMetadata,
                       path: str,
                      ) -> int :
    """Evaluates the volume difference function at some point. 

       The volume-difference function is a function with domain which is the
           set of all possible traction vectors and range which is all positive 
           reals. It is defined for a particular target geometry, a particular
           deformed geometry, and a face. To evaluate the function at a point, 
           you simply apply the traction (the point) to the target geometry on 
           the specified face, overlap the resulting geometry in space with the 
           deformed geometry, and compute the volume of the pieces which are not
           in the intersection of the volumes (i.e. the symmetric difference).

       Implementation Details:
       By applying traction to the target geometry, and not the deformed geometry,
           there is no uncertainty in the location of the trench (which may be
           in a different location due to deformation).
                   
       Args:
           at:                     The traction vector to apply. This is the point
                                       at which the function is evaluated. 
           point_on_face:          Point on the face on which to apply the tractions.
           deformed_geometry_data: The real world data collected by scanning
                                       the workpiece after it underwent some
                                       deformation due to material removal.
           target_geometry_data:   The real world data collected by scanning
                                       the workpiece before it underwent some
                                       deformation due to material removal.
           commit_phase_md:        Metadata associated with the commitment phase
                                       to which the deformed MDB belongs. 
           path:                   Absolute path to a directory where simulation
                                       results will be dumped.
    
       Returns:
           The value of the volume difference function (a scalar representing the
               volume) when evaluated at the point of interest.
    
       Raises:
           None.
    """
    
    # Step 1:
    # Use the data for the target geometry to produce a MDB in which the
    #     traction application can happen.
    
    # Change the CWD to fix location of simulation results.
    orig_cwd = os.getcwd()
    os.chdir(path)

    # Step 2:
    # Simulate the traction application.
    job_name = _simulate_traction_app(at, point_near_face, shim.STANDARD_MODEL_NAME,
                                      target_geometry_mdb_md, commit_phase_md, 
                                      target_geometry_mdb)
    # Restore CWD.
    os.chdir(orig_cwd)

    # Step 3:
    # Clean up. Now the ODB produced by the job can be used.
    shim.save_mdb(target_geometry_mdb)
    shim.close_mdb(target_geometry_mdb)

    post_traction_odb = os.path.join(path, job_name + ".odb")

    # Step 4:
    # Map the results from the ODB (produced by applying a traction to the target
    #     geometry) to a .stl. 

    # Step 5:
    # Compute the symmetric difference in volume between the two parts. 

    return symm_diff_vol 



# *****************************************************************************
#                                 DEPRECATED 
# *****************************************************************************



def _compute_symmetric_diff(part1_name: str, part2_name: str, model_name: str,
                            mdb_md: abq_md.AbaqusMdbMetadata, mdb: Any) -> int:
    """DEPRECATED. This is now being done in Blender due to Abaqus' inflexible
           front end.

       Computes the symmetric difference in volume between two part geometries. 

       In set theory, the symmetric difference between two sets A, B is the set
           of elements which are in either A, B and are not in the intersection
           of A, B. This function computes the volume of the symmetric difference. 
           
                   
       Args:
           part1_name: Name of the first part.
           part2_name: Name of the second part.
           model_name: Name of the model containing both parts. Assumed that the
                           model is in its default state, aside from containing
                           two parts.
           mdb_md:     Metadata associated with the MDB.
           mdb:        Abaqus MDB object. Assumed to be open.
    
       Returns:
           Volume of the symmetric difference.
    
       Raises:
           None.
    """

    # Check that the model is in the expected state.
    if not shim.check_two_part_model(True, model_name, mdb):
        raise RuntimeError("Model in unexpected state.")
    
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
    PART1_REMAINDER_PART_NAME = part1_name + "_Remainder"
    part1_remainder = shim.cut_instances_in_assembly(PART1_REMAINDER_PART_NAME, part1_instance,
                                                     (part2_instance,), model_name, mdb_md,
                                                     mdb)
    
    remainder_part_instance = shim.instance_part_into_assembly(PART1_REMAINDER_PART_NAME, 
                                                               part1_remainder, False, 
                                                               model_name, mdb)
    
    # Unsuppress the instance of part1 so it can be used in the subsequent cutting
    #     operation.
    shim.resume_assembly_feature(part1_name, mdb.models[model_name])
    part1_instance = mdb.models[model_name].rootAssembly.instances[part1_name]

    # Do the cut to get the intersection.
    INTERSECTION_PART_NAME = "Intersection"
    intersection = shim.cut_instances_in_assembly(INTERSECTION_PART_NAME, part1_instance, 
                                                  (remainder_part_instance,), model_name,
                                                  mdb_md, mdb)

    # Step 3:
    # Find the symmetric difference by using the intersection.

    intersection_instance = shim.instance_part_into_assembly(INTERSECTION_PART_NAME, intersection, False, model_name, mdb)

    # Resume the two part instances and merge them. The resulting part is the
    #     union.
    shim.resume_assembly_feature(part1_name, mdb.models[model_name])
    part1_instance = mdb.models[model_name].rootAssembly.instances[part1_name]
    shim.resume_assembly_feature(part2_name, mdb.models[model_name])
    part1_instance = mdb.models[model_name].rootAssembly.instances[part2_name]

    UNION_PART_NAME = "Union"
    union = shim.merge_instances_in_assembly(UNION_PART_NAME, (part1_instance, part2_instance), 
                                             False, model_name, mdb_md, mdb)
    union_instance = shim.instance_part_into_assembly(UNION_PART_NAME, union, False, model_name, mdb)
    
    # Remove the intersection from the union to yield the symmetric difference.
    SYMM_DIFF_PART_NAME = "Symmetric_Difference"
    symmetric_difference = shim.cut_instances_in_assembly(SYMM_DIFF_PART_NAME, union_instance,
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
