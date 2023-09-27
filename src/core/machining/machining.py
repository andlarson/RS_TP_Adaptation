import shutil
import os

import core.part.part as part
import core.simulation.simulation as sim
import core.boundary_conditions.boundary_conditions as bc
import core.metadata.metadata as md 
import core.metadata.abaqus_metadata as abq_md
import core.metadata.naming as naming
import core.tool_pass.tool_pass as tp
import core.abaqus.abaqus_shim as shim
from util.debug import *


STANDARD_POST_COMMIT_FILE_NAME_PREFIX = "Post_Committed_Tool_Pass_"


class MachiningProcess:

    # The object which defines the whole machining process for a single part.
    # 
    # Notes:
    #    Does a one-time switch of the CWD. If init_part already has an underlying
    #       MDB, then the CWD is switched to the directory that the MDB already
    #       lives in. If init_part does not have an underlying MDB, an MDB is
    #       created at some known location (TODO: Fill in this behavior) and the
    #       CWD is switched to the directory of that this new MDB lives in. This
    #       is the only place that the CWD is permanently changed.
    #
    # Arguments:
    #    totally_in_simulation - Boolean.
    #                            Flag which indicates if the whole machining process
    #                               is being conducted in simulation i.e. there is
    #                               no real-world machining process going on.
    #    init_part             - Part object.
    #                            The initial part geometry before any machining has
    #                               taken place i.e. the block of material from which
    #                               the part will be machined.
    #    boundary_conditions   - List of BC objects.
    #                            The boundary conditions which should be used for
    #                               all simulations in this machining process. Usually,
    #                               these reflect the clamping conditions of the
    #                               part. Without any boundary conditions, a finite
    #                               element simulation is unconstrained. 
    def __init__(self, init_part, boundary_conditions):
    # type: (bool, part.Part, List[bc.BC]) -> None
  
        assert(len(boundary_conditions) > 0)

        os.chdir(os.path.dirname(init_part.path_to_mdb))

        self.metadata_committed_tool_pass_plans = []
        self.stress_profile_estimates = []
        self.boundary_conditions = boundary_conditions 

        if isinstance(init_part, part.UserDefinedPart):
            # When this is supported, it's expected that a UserDefinedPart is
            #    first constructed in an MDB so the method of setting up the
            #    metadata is the same as the AbaqusDefinedPart.
            raise RuntimeError("Not yet supported.")

        elif isinstance(init_part, part.AbaqusDefinedPart):
            first_tp_metadata = md.CommittedToolPassPlanMetadata(init_part, init_part.path_to_mdb, boundary_conditions)
            self.metadata_committed_tool_pass_plans.append(first_tp_metadata)

        else:
            raise RuntimeError("Unrecognized part type.")



    # Commit to a tool pass plan. If this exact tool pass plan was already simulated
    #    in this commitment phase, no additional simulations are done. If this
    #    tool pass plan has not yet been simulated, it is simulated. 
    #
    # Notes:
    #    This function has two effects on the file system:
    #       This function saves the simulation artifacts (the .odb, .sim, .inp, 
    #          etc. and the final .cae file) in a subdirectory of the directory 
    #          that the MDB for this commitment phase lives in.
    #       This function saves the .cae file which contains the part geometry
    #          resulting from the committed tool passes in the directory that 
    #          the MDB for this commitment phase lives in.
    #
    # Arguments:
    #    tool_pass_plan - ToolPassPlan object.
    #    save_name      - String.
    #                     Name of the subdirectory, the .cae file in the subdirectory.
    #
    # Returns:
    #    None. 
    def commit_tool_passes(self, tool_pass_plan, save_name):
    # type: (tp.ToolPassPlan, str) -> None

        # For the first commitment phase, the user must supply a stress profile.
        if len(self.metadata_committed_tool_pass_plans) == 0 and \
           len(self.stress_profile_estimates) == 0:
            raise AssertionError("A user-specified stress profile must be supplied \
                                  for the first commitment phase!")

        # If the user did not supply a stress profile estimate for this commitment 
        #    phase, record as such.
        if len(self.stress_profile_estimates) < len(self.metadata_committed_tool_pass_plans):
            self.stress_profile_estimates.append(None)

        self.metadata_committed_tool_pass_plans[-1].committed_tool_pass_plan = tool_pass_plan 

        num_commits = len(self.metadata_committed_tool_pass_plans)

        # If the tool pass matches one which was already simulated, just duplicate
        #    the results of the previous simulation into a new directory. 
        match = False
        for name, already_simulated_plan in self.metadata_committed_tool_pass_plans[-1].simulated_tool_pass_plans:
            if tp.compare_tool_pass_plans(tool_pass_plan, already_simulated_plan):
                shutil.copytree(name, save_name)
                
                # Rename the .cae and .jnl files in this new subdirectory to conform.
                shutil.move(os.path.join(save_name, name + ".cae"), os.path.join(save_name, save_name + ".cae"))
                shutil.move(os.path.join(save_name, name + ".jnl"), os.path.join(save_name, save_name + ".jnl"))

                match = True

                dp("For committed tool pass number " + str(num_commits) + ", the tool pass plan was already simulated so no additional simulation was needed!")

        # Otherwise, simulate the plan.
        if not match:
            self.sim_potential_tool_passes(tool_pass_plan, save_name)

        # Create the MDB for the next commitment phase and the metadata that
        #    accompanies it.
        new_mdb_name = STANDARD_POST_COMMIT_FILE_NAME_PREFIX + str(num_commits)
        mdb = shim.create_mdb(new_mdb_name, os.getcwd()) 
        mdb_metadata = abq_md.AbaqusMdbMetadata(new_mdb_name) 

        # Generate the names for the stuff in the MDB.
        names = naming.new_model_names(mdb_metadata, True)

        # Create a part in the pre-existing lone model which comes with a new
        #    MDB.
        odb_name = shim.STANDARD_JOB_PREFIX + str(len(tool_pass_plan.plan)) + ".odb"
        odb_path = os.path.join(save_name, odb_name)
        shim.create_part_from_odb(names["pre_tool_pass_part_name"], names["new_model_name"], odb_path, mdb_metadata, mdb)

        # Map the orphan mesh to a part geometry.
        sim.orphan_mesh_to_geometry(names["pre_tool_pass_part_name"], names["new_model_name"], mdb)

        # Remove redundant features.
        new_part = shim.get_part(names["pre_tool_pass_part_name"], names["new_model_name"], mdb)
        shim.add_virtual_topology(new_part)

        # Propagate the material definitions and sections from the ODB.
        shim.create_material_from_odb(odb_path, names["new_model_name"], mdb)
        shim.create_section_from_odb(odb_path, names["new_model_name"], mdb)

        # Do the necessary section assignment.
        shim.assign_section_to_whole_part(names["pre_tool_pass_part_name"], names["new_model_name"], mdb)

        # Save the MDB so it is visible.
        save_path = os.path.join(os.getcwd(), new_mdb_name)
        shim.save_mdb_as(save_path, mdb)

        # The starting point for the next commitment phase.
        abaqus_part = part.AbaqusDefinedPart(save_name, save_path)

        # The initial state of the next commitment phase depends on the result
        #    of the simulation. 
        self.metadata_committed_tool_pass_plans.append(md.CommittedToolPassPlanMetadata(abaqus_part, abaqus_part.path_to_mdb, self.boundary_conditions))



    # Simulate some potential tool passes and save off the results. 
    #
    # Notes:
    #    This function has a single effect on the file system:
    #       This function saves the simulation artifacts (the .odb, .sim, .inp, 
    #          etc. and the final .cae file) in a subdirectory of the directory 
    #          that the MDB for this commitment phase lives in.
    #
    # Arguments:
    #    tool_pass_plan - ToolPassPlan object.
    #    name           - String.
    #                     The name of the subdirectory in which the simulation artifacts
    #                        will be placed. Also, the name of the MDB which results
    #                        from these simulations and lives in the subdirectory. 
    #
    # Returns:
    #    None. 
    def sim_potential_tool_passes(self, tool_pass_plan, save_name):
    # type: (tp.ToolPassPlan, str) -> None

        # For the first commitment phase, the user must supply a stress profile.
        if len(self.metadata_committed_tool_pass_plans) == 0 and \
           len(self.stress_profile_estimates) == 0:
            raise AssertionError("A user-specified stress profile must be supplied \
                                  for the first commitment phase!")

        # The start point for each sequence of simulations is always the MDB which
        #    resulted from the last committed tool pass (or the very initial MDB).
        path_to_mdb = self.metadata_committed_tool_pass_plans[-1].path_initial_mdb

        # A new MDB is created for this sequence of simulations. Therefore, a
        #    new metadata data structure must exist and accompany this new MDB.
        self.metadata_committed_tool_pass_plans[-1].per_mdb_metadata.append(abq_md.AbaqusMdbMetadata(path_to_mdb))

        sim.sim_tool_pass_plan(tool_pass_plan, save_name, self.metadata_committed_tool_pass_plans[-1])



    # Provide a stress profile as the start point for the current commitment
    #    phase.
    # By default, the stress profile for the start point of the current commitment
    #    phase comes from the output of last simulation in the previous commitment
    #    phase. This function overrides that default behavior!!
    #
    # Notes:
    #    At the beginning of each commitment phase, the part needs to have some stress
    #       profile associated with it. By default, the part's stress profile is sourced
    #       from the output of the last simulation in the previous commitment phase.
    #       This doesn't work in the first commitment phase. Also, it's possible to
    #       glean information from deformations observed in real-life to improve stress
    #       profile estimates.
    #    
    # Arguments:
    #    path - String.
    #           Path to the stress subroutine which represents the stress profile. 
    #    
    # Return:
    #    None.
    def record_estimated_stress_profile(self, path):

        self.stress_profile_estimates.append(path)



    def estimate_stress_via_last_tool_pass(self):
    # type: (None) -> None
        
        raise RuntimeError("Not yet supported.")