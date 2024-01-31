"""
This module contains top-level, user-facing functionality for:
    1) Performing simulations of potential tool passes.
    2) Recovering information about the residual stress field.

The MachiningProcess object is the top-level object that a user of this library
    should create. This object's methods implement the functionality listed
    above.
"""

import shutil
import os

import src.core.part.part as part
import src.core.simulation.simulation as sim
import src.core.boundary_conditions.boundary_conditions as bc
import src.core.metadata.metadata as md 
import src.core.metadata.abaqus_metadata as abq_md
import src.core.metadata.naming as naming
import src.core.tool_pass.tool_pass as tp
import src.core.abaqus.abaqus_shim as shim
import src.core.material_properties.material_properties as mp

from src.util.debug import *


STANDARD_POST_COMMIT_FILE_NAME_PREFIX = "Post_Commit_"


class MachiningProcess:

    def __init__(self, init_part: part.Part, boundary_conditions: list[bc.BC]) -> None:
        """Defines the whole machining process for a single part.
   
           Implementation Detail:
           Does a one-time switch of the CWD. If the initial part geometry already 
               has an underlying MDB, then the CWD is switched to the directory 
               that the MDB already lives in. If the initial part geometry does not 
               have an underlying MDB, an MDB is created at some known location 
               (TODO: Fill in this behavior) and the CWD is switched to the directory 
               that this new MDB lives in. This is the only place that the CWD 
               is permanently changed.
           
           Args:
               init_part:           The initial part geometry before any machining 
                                        has taken place (e.g. the blank). 
               boundary_conditions: The boundary conditions which should be used for all 
                                        simulations in this machining process. Usually, 
                                        these reflect the clamping conditions of the part. 
                                        Without any boundary conditions, a finite element 
                                        simulation is unconstrained. Thus, there should
                                        always be at least one boundary condition.
           
           Returns:
               None.
           
           Raises:
               RuntimeError: Fewer than one boundary conditions were specified.
        """
     
        if (len(boundary_conditions) == 0):
            raise RuntimeError("At least one boundary condition must be supplied.")

        self.metadata_committed_tool_pass_plans = []
        self.stress_profile_estimates = []
        self.boundary_conditions = boundary_conditions 

        if isinstance(init_part, part.UserDefinedPart):
            # When this is supported, it's expected that a UserDefinedPart is
            #    first constructed in an MDB so the method of setting up the
            #    metadata is the same as the AbaqusDefinedPart.
            raise AssertionError("Not yet supported.")

        elif isinstance(init_part, part.AbaqusDefinedPart):
            first_tp_metadata = md.CommittedToolPassPlanMetadata(init_part, init_part.path_to_mdb, boundary_conditions)
            self.metadata_committed_tool_pass_plans.append(first_tp_metadata)

            os.chdir(os.path.dirname(init_part.path_to_mdb))

        else:
            raise AssertionError("Unrecognized part type.")



    def commit_tool_passes(self, tool_pass_plan: tp.ToolPassPlan, save_name: str) -> None:
        """Commits to a tool pass plan. 
           
           It is assumed that a user of this library does some simulations of
               potential tool passes and then determines which tool passes
               that they are content with. In this case, the user will commit
               to these tool passes, and if there is a real machining process
               going on, do these tool passes in real life. Committing to a
               tool pass plan is just the way that the user tells this library
               that they are content with this set of tool passes. 
           
           Implementation Details:
           If this exact tool pass plan was already simulated in this commitment 
               phase, no additional simulations are done. If this tool pass plan 
               has not yet been simulated, it is simulated. 
    
           This function has two effects on the file system:
               This function saves the simulation artifacts (the .odb, .sim, .inp, 
                   etc. and the final .cae file) in a subdirectory of the directory 
                   that the MDB for this commitment phase lives in.
               This function saves the .cae file which contains the part geometry
                   resulting from the committed tool passes in the directory that 
                   the MDB for this commitment phase lives in.
           
           Args:
               tool_pass_plan: The tool pass plan to commit to. 
               save_name:      Name of the subdirectory that this function creates 
                                   and the .cae file which is created in the 
                                   subdirectory.
           
           Returns:
               None. 
           
           Raises:
               RuntimeError: No stress profile was supplied for the first commitment
                                 phase. If no stress profile is specified for the
                                 first commitment phase, no deformation will occur
                                 when tool paths are simulated.
        """

        # For the first commitment phase, the user must supply a stress profile.
        if len(self.metadata_committed_tool_pass_plans) == 0 and \
           len(self.stress_profile_estimates) == 0:
            raise RuntimeError("A user-specified stress profile must be supplied \
                                  for the first commitment phase!")

        # If the user did not supply a stress profile estimate for this commitment 
        #    phase, record as such.
        if len(self.stress_profile_estimates) < len(self.metadata_committed_tool_pass_plans):
            self.stress_profile_estimates.append(None)

        self.metadata_committed_tool_pass_plans[-1].committed_tool_pass_plan = tool_pass_plan 

        self._lazy_sim_potential_tool_passes(tool_pass_plan, save_name)
        self._prepare_next_commit(tool_pass_plan, save_name)



    def sim_potential_tool_passes(self, tool_pass_plan: tp.ToolPassPlan, 
                                  save_name: str, save_dir: str) -> None:
        """Simulates potential tool passes and saves off the results. 
        
           Args:
               tool_pass_plan: The tool pass plan to simulate.
               save_dir:       Absolute path to a directory. In this directory,
                                   a subdirectory with name save_name will be 
                                   created.
               save_name:      The name of the subdirectory in which the simulation 
                                   artifacts (.odb file, .sim file, etc. and the 
                                   final .cae file) will be placed. Also, the name 
                                   of the artifacts themselves.
                               If a file or directory exists at the path specified
                                   by save_dir and save_name (concatenated), it will
                                   be permanently deleted to make room for the new
                                   subdirectory.
        
           Returns:
               None. 
           
           Raises:
               RuntimeError: No stress profile was supplied for the first commitment
                                 phase. If no stress profile is specified for the
                                 first commitment phase, no deformation will occur
                                 when tool paths are simulated.
        """

        # For the first commitment phase, the user must supply a stress profile.
        if len(self.metadata_committed_tool_pass_plans) == 0 and \
           len(self.stress_profile_estimates) == 0:
            raise RuntimeError("A user-specified stress profile must be supplied \
                                  for the first commitment phase!")

        # The start point for each sequence of simulations is always the MDB which
        #    resulted from the last commitment phase (or the initial MDB in the
        #    very first commitment phase).
        path_to_mdb = self.metadata_committed_tool_pass_plans[-1].path_initial_mdb

        # A new MDB is created for this sequence of simulations. Therefore, a
        #    new metadata data structure must exist and accompany this new MDB.
        mdb_metadata = abq_md.AbaqusMdbMetadata(path_to_mdb)
        self.metadata_committed_tool_pass_plans[-1].per_mdb_metadata.append(mdb_metadata)

        # If there is an estimated stress profile for this commitment phase, use it.
        if len(self.metadata_committed_tool_pass_plans) == len(self.stress_profile_estimates):
            sim.sim_tool_pass_plan(tool_pass_plan, save_name, save_dir, self.metadata_committed_tool_pass_plans[-1], 
                                   self.stress_profile_estimates[-1])
        else:
            sim.sim_tool_pass_plan(tool_pass_plan, save_name, save_dir, self.metadata_committed_tool_pass_plans[-1])



    def record_estimated_stress_profile(self, path: str):
        """Records a stress profile to be used as the initial stress state for
               the current commitment phase.

           By default, the stress profile for the start point of the current 
               commitment phase comes from the output of last simulation in 
               the previous commitment phase. This function overrides that 
               default behavior!!
           
           At the beginning of each commitment phase, the part needs to have some 
               stress profile associated with it. By default, the part's stress 
               profile is sourced from the output of the last simulation in the 
               previous commitment phase. This doesn't work in the first 
               commitment phase. Also, it's possible to glean information from 
               deformations observed in real-life to improve stress profile 
               estimates.
           
           Each commitment phase should have, at most, one initial stress state.
              
           Args:
               path: Path to the object file produced by compiling the stress
                         subroutine which defines the desired stress profile.
              
           Return:
               None.

           Raises:
               None.
        """

        if len(self.stress_profile_estimates) >= len(self.metadata_committed_tool_pass_plans):
            raise RuntimeError("Trying to pass too many estimated stress \
                                profiles. Each commitment phase should have, at \
                                most, one accompanying estimated stress profile!")

        self.stress_profile_estimates.append(path)



    def estimate_stress_via_last_tool_pass(self):
        """Will implement the technique to recover stress information due to
              deformations which result from a tool pass."""
        
        raise RuntimeError("Not yet supported.")



    def _lazy_sim_potential_tool_passes(self, tool_pass_plan: tp.ToolPassPlan, save_name: str) -> None:
        """Simulates a tool pass plan if it hasn't already been simulated. If the 
               toolpass plan was already simulated, copies the results into a new 
               directory and the names are updated accordingly. 
           
           Implementation Details:
           Assumes that CWD contains the directories of simulation results.
           
           Args:
              tool_pass_plan: The tool pass plan which may or may not need to be 
                                  re-simulated.
              save_name:      The name of the simulation results.
           
           Returns:
               None. 

           Raises:
               None.
        """

        match = False
        for name, already_simulated_plan in self.metadata_committed_tool_pass_plans[-1].simulated_tool_pass_plans:
            if tp.compare_tool_pass_plans(tool_pass_plan, already_simulated_plan):
                shutil.copytree(name, save_name)
                
                # Rename the .cae and .jnl files in this new subdirectory.
                shutil.move(os.path.join(save_name, name + ".cae"), os.path.join(save_name, save_name + ".cae"))
                shutil.move(os.path.join(save_name, name + ".jnl"), os.path.join(save_name, save_name + ".jnl"))

                match = True

                num_commits = len(self.metadata_committed_tool_pass_plans)
                dp("For commit " + str(num_commits) + ", the tool pass plan was already simulated so no additional simulation was needed!")
                break

        # Simulate the plan if necessary.
        if not match:
            self.sim_potential_tool_passes(tool_pass_plan, save_name)



    def _prepare_next_commit(self, tool_pass_plan: tp.ToolPassPlan, save_name: str) -> None:
        """Prepares for the next commitment phase by propagating the results
               of the last commitment phase.

           In particular, this means taking the geometry which resulted from
               the last tool pass in the previous commitment phase and using
               that as the initial geometry for the new, next commitment
               phase.
           
           Implementation Details:
           Assumes that CWD contains the directories of simulation results.
        
           Args:
              tool_pass_plan: The tool pass plan which may or may not need to be re-simulated.
              save_name: The desired name of the simulation results.
        
           Returns:
               None.

           Raises:
               None.
        """

        # Record the path of the .sim file. The .sim file may be used in the
        #     next commitment phase to set the initial stress state.
        num_jobs = len(tool_pass_plan.plan)
        sim_file_name = shim.STANDARD_JOB_PREFIX + str(num_jobs) + ".sim"
        sim_file_path = os.path.join(os.getcwd(), save_name, sim_file_name) 

        # Create the MDB for the next commitment phase and the metadata that
        #     accompanies it.
        num_commits = len(self.metadata_committed_tool_pass_plans)
        new_mdb_name = STANDARD_POST_COMMIT_FILE_NAME_PREFIX + str(num_commits)
        mdb = shim.create_mdb(new_mdb_name, os.getcwd()) 
        mdb_metadata = abq_md.AbaqusMdbMetadata(new_mdb_name) 

        # Generate the names for the stuff in the MDB.
        names = naming.ModelNames(mdb_metadata, True)

        # Create a part in the pre-existing lone model which comes with a new
        #     MDB.
        odb_name = shim.STANDARD_JOB_PREFIX + str(len(tool_pass_plan.plan)) + ".odb"
        odb_path = os.path.join(save_name, odb_name)
        shim.create_part_from_odb(names.pre_tool_pass_part_name, names.new_model_name, odb_path, mdb_metadata, mdb)

        # Map the orphan mesh to a part geometry.
        sim.orphan_mesh_to_geometry(names.pre_tool_pass_part_name, names.new_model_name, mdb)

        # Save the MDB so it is visible.
        save_path = os.path.join(os.getcwd(), new_mdb_name)
        shim.save_mdb_as(save_path, mdb)

        # Create the material based on the material of the very first part in
        #     the machining process. 
        first_commit_phase_metadata = self.metadata_committed_tool_pass_plans[0]
        very_first_part = first_commit_phase_metadata.init_part
        material = very_first_part.material

        # The starting point for the next commitment phase.
        abaqus_part = part.AbaqusDefinedPart(save_name, save_path, material)

        # The initial state of the next commitment phase depends on the result
        #     of the simulation. 
        new_commit_metadata = md.CommittedToolPassPlanMetadata(abaqus_part, abaqus_part.path_to_mdb, self.boundary_conditions)
        self.metadata_committed_tool_pass_plans.append(new_commit_metadata)

        # Record the stress state after the last tool pass in the previous commit.
        self.metadata_committed_tool_pass_plans[-1].path_last_commit_sim_file = sim_file_path

