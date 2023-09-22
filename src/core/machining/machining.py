import shutil
import os

import core.part.part as part
import core.simulation.simulation as sim
import core.boundary_conditions.boundary_conditions as bc
import core.metadata.metadata as md 
import core.metadata.abaqus_metadata as abq_md
import core.tool_pass.tool_pass as tp
import core.abaqus.abaqus_shim as shim


STANDARD_POST_COMMIT_FILE_NAME_PREFIX = "Post_Committed_Tool_Pass_"


class MachiningProcess:
   
    def __init__(self, init_part, boundary_conditions):
    # type: (part.Part, List[bc.BC]) -> None
  
        self.metadata = []
        if isinstance(init_part, part.UserDefinedPart):
            # When this is supported, it's expected that a UserDefinedPart is
            #    first constructed in an MDB so the method of setting up the
            #    metadata is the same as the AbaqusDefinedPart.
            raise RuntimeError("Not yet supported.")

        elif isinstance(init_part, part.AbaqusDefinedPart):
            first_tp_metadata = md.CommittedToolPassPlanMetadata(init_part, init_part.path_to_mdb, boundary_conditions)
            self.metadata.append(first_tp_metadata)

        # The boundary conditions are assumed to be fixed for the whole machining process. 
        self.boundary_conditions = boundary_conditions 



    # Commit to a tool pass plan. If this exact tool pass plan was already simulated
    #    in this commitment phase, no additional simulations are done. If this
    #    tool pass plan has not yet been simulated, it is simulated. 
    #
    # Notes:
    #    This function has two effects on the file system:
    #       This function saves the simulation artifacts (the .odb, .sim, .inp, 
    #          etc. and the final .cae file) in a subdirectory of the directory 
    #          that the MDB for this commitment phase lives in.
    #       This function also saves the .cae file which represents the MDB for
    #          the next commitment phase in the directory that the MDB for this
    #          commitment phase lives in.
    #
    # Arguments:
    #    tool_pass_plan - ToolPassPlan object.
    #    save_name      - String.
    #
    # Returns:
    #    None. 
    def commit_tool_passes(self, tool_pass_plan, save_name):
    # type: (tp.ToolPassPlan, str) -> None

        self.metadata[-1].committed_tool_pass_plan = tool_pass_plan 

        # If the tool pass matches one which was already simulated, just duplicate
        #    the results of the previous simulation into a new directory. 
        match = False
        for name, simd_tpp in self.metadata[-1].simulated_tool_pass_plans:
            if tp.compare_tool_pass_plans(tool_pass_plan, simd_tpp):
                shutil.copytree(name, save_name)
                match = True

        # Otherwise, simulate the plan.
        if not match:
            self.sim_potential_tool_passes(tool_pass_plan, save_name)

        # Create the MDB for the next commitment phase and the metadata that
        #    accompanies it.
        num_commits = len(self.metadata)
        new_mdb_name = STANDARD_POST_COMMIT_FILE_NAME_PREFIX + str(num_commits)
        mdb = shim.create_mdb(new_mdb_name, os.getcwd()) 
        mdb_metadata = abq_md.AbaqusMdbMetadata(new_mdb_name) 

        # The last tool pass in the tool pass plan generated a .odb file which 
        #    needs to be mapped to a part geometry and placed into the MDB
        #    that was just created.
        odb_path = save_name + "/" + shim.STANDARD_JOB_PREFIX + str(len(tool_pass_plan.plan)) 
        shim.create_model_from_odb(odb_path, "Model-1", mdb_metadata, mdb)

        # The starting point for the next commitment phase.
        abaqus_part = part.AbaqusDefinedPart(save_name, new_mdb_name)

        # The initial state of the next commitment phase depends on the result
        #    of the simulation. 
        self.metadata.append(md.CommittedToolPassPlanMetadata(abaqus_part, abaqus_part.path_to_mdb, self.boundary_conditions))

        # Save the MDB so it is visible.
        shim.save_mdb(mdb)



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

        # The start point for each sequence of simulations is always the MDB which
        #    resulted from the last committed tool pass (or the very initial MDB).
        path_to_mdb = self.metadata[-1].path_initial_mdb

        # A new MDB is created for this sequence of simulations. Therefore, a
        #    new metadata data structure must exist and accompany this new MDB.
        self.metadata[-1].per_mdb_metadata.append(abq_md.AbaqusMdbMetadata(path_to_mdb))

        # Set the CWD to the directory where the MDB for this commitment phase
        #    lives in.
        os.chdir(self.metadata[-1].path_initial_mdb)

        sim.sim_consecutive_tool_passes(tool_pass_plan, save_name, self.metadata[-1])



    def estimate_stress_via_last_tool_pass(self):
    # type: (None) -> None
        
        raise RuntimeError("Not yet supported.")