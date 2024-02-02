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
import pathlib

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



class MachiningProcess:

    def __init__(self, init_part: part.Part, boundary_conditions: list[bc.BC]) -> None:
        """Defines the whole machining process for a single part.
   
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
               RuntimeError: At least one boundary condition must be specified.
        """
     
        if len(boundary_conditions) == 0:
            raise RuntimeError("At least one boundary condition must be supplied.")

        self.commitment_phase_metadata = []
        self.stress_profile_estimates = []
        self.boundary_conditions = boundary_conditions 

        if isinstance(init_part, part.UserDefinedPart):
            # When this is supported, it's expected that a UserDefinedPart is
            #    first constructed in an MDB so the method of setting up the
            #    metadata is the same as the AbaqusDefinedPart.
            raise AssertionError("Not yet supported.")
        elif isinstance(init_part, part.AbaqusDefinedPart):
            first_tp_metadata = md.CommittedToolPassPlanMetadata(init_part, init_part.path_to_mdb, boundary_conditions)
            first_tp_metadata.first_commitment_phase = True
            self.commitment_phase_metadata.append(first_tp_metadata)
        else:
            raise AssertionError("Unrecognized part type.")



    def commit_tool_passes(self, tool_pass_plan: tp.ToolPassPlan, save_name: str,
                           save_dir: str) -> None:
        """Commits to a tool pass plan. 
           
           It is assumed that a user of this library does some simulations of
               potential tool passes and then determines which tool passes
               that they are content with. In this case, the user will commit
               to these tool passes, and if there is a real machining process
               going on, do these tool passes in real life. Committing to a
               tool pass plan is just the way that the user tells this library
               that they are content with this set of tool passes. 

           TODO: When real-life data is available, the user just needs to supply 
               the geometries of the tool passes which were performed and the
               real-life shape of the part after the tool pass was performed.
           
           Implementation Details:
           If this exact tool pass plan was already simulated in this commitment 
               phase, no additional simulations are done. If this tool pass plan 
               has not yet been simulated, it is simulated. 
    
           Args:
               tool_pass_plan: The tool pass plan to commit to. 
               save_name:      The name of the subdirectory in which the simulation 
                                   artifacts (.odb file, .sim file, etc. and the 
                                   final .cae file) will be placed. Also, the name 
                                   of the artifacts themselves.
               save_dir:       Absolute path to a directory. In this directory,
                                   a subdirectory with name save_name will be 
                                   created.
                               If this path is already occupied (i.e. something
                                   exists there) the behavior of this function is
                                   undefined. 
           
           Returns:
               None. 
           
           Raises:
               RuntimeError: No stress profile was supplied for the first commitment
                                 phase. If no stress profile is specified for the
                                 first commitment phase, no deformation will occur
                                 when tool paths are simulated.
        """

        # For the first commitment phase, the user must supply a stress profile.
        if len(self.commitment_phase_metadata) == 1 and \
           len(self.stress_profile_estimates) == 0:
            raise RuntimeError("A user-specified stress profile must be supplied \
                                  for the first commitment phase!")

        # If necessary, start a new commitment phase.
        if self._is_new_commitment_phase():
            self._start_new_commitment_phase(save_name, save_dir) 

        # If the user did not supply a stress profile estimate for this commitment 
        #    phase, record as such.
        if len(self.stress_profile_estimates) < len(self.commitment_phase_metadata):
            self.stress_profile_estimates.append(None)
        
        last_commit_metadata = self.commitment_phase_metadata[-1]
        last_commit_metadata.committed_tool_pass_plan = tool_pass_plan 
        last_commit_metadata.committed_tool_pass_plan_path = os.path.join(save_name, save_dir)

        self._lazy_sim_potential_tool_passes(tool_pass_plan, save_name, save_dir)



    def sim_potential_tool_passes(self, tool_pass_plan: tp.ToolPassPlan, 
                                  save_name: str, save_dir: str) -> None:
        """Simulates potential tool passes and saves off the results. 
        
           Args:
               tool_pass_plan: The tool pass plan to simulate.
               save_name:      The name of the subdirectory in which the simulation 
                                   artifacts (.odb file, .sim file, etc. and the 
                                   final .cae file) will be placed. Also, the name 
                                   of the artifacts themselves.
               save_dir:       Absolute path to a directory. In this directory,
                                   a subdirectory with name save_name will be 
                                   created.
                               If this path is already occupied (i.e. something
                                   exists there) the behavior of this function is
                                   undefined. 
        
           Returns:
               None. 
           
           Raises:
               RuntimeError: No stress profile was supplied for the first commitment
                                 phase. If no stress profile is specified for the
                                 first commitment phase, no deformation will occur
                                 when tool paths are simulated.
        """

        # For the first commitment phase, the user must supply a stress profile.
        if len(self.commitment_phase_metadata) == 1 and \
           len(self.stress_profile_estimates) == 0:
            raise RuntimeError("A user-specified stress profile must be supplied \
                                  for the first commitment phase!")

        # If necessary, start a new commitment phase.
        if self._is_new_commitment_phase():
            self._start_new_commitment_phase(save_name, save_dir) 

        # The start point for each sequence of simulations is always the initial
        #     state of the current commitment phase. 
        commit_metadata = self.commitment_phase_metadata[-1]
        path_to_mdb = commit_metadata.path_initial_mdb

        # A new MDB is created for this sequence of simulations. Therefore, a
        #    new metadata data structure must exist and accompany this new MDB.
        mdb_metadata = abq_md.AbaqusMdbMetadata(path_to_mdb)
        commit_metadata.per_mdb_metadata.append(mdb_metadata)

        # If there is an estimated stress profile for this commitment phase, use it.
        if len(self.commitment_phase_metadata) == len(self.stress_profile_estimates):
            sim.sim_tool_pass_plan(tool_pass_plan, save_name, save_dir, 
                                   commit_metadata, self.stress_profile_estimates[-1])
        else:
            sim.sim_tool_pass_plan(tool_pass_plan, save_name, save_dir, commit_metadata)



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

        if len(self.stress_profile_estimates) >= len(self.commitment_phase_metadata):
            raise RuntimeError("Trying to pass too many estimated stress \
                                profiles. Each commitment phase should have, at \
                                most, one accompanying estimated stress profile!")

        self.stress_profile_estimates.append(path)



    def estimate_stress_via_last_tool_pass(self):
        """
            
           Args:
               None.

           Returns:
               None.

           Raises:
               None.
        """
        
        raise RuntimeError("Not yet supported.")



    def _lazy_sim_potential_tool_passes(self, tool_pass_plan: tp.ToolPassPlan, 
                                        save_name: str, save_dir: str) -> None:
        """Simulates a tool pass plan if it hasn't already been simulated. If the 
               toolpass plan was already simulated, copies the results into a new 
               directory and the names are updated accordingly. 
           
           Args:
               tool_pass_plan: The tool pass plan which may or may not need to be 
                                  re-simulated.
               save_name:      The name of the subdirectory in which the simulation 
                                   artifacts (.odb file, .sim file, etc. and the 
                                   final .cae file) will be placed. Also, the name 
                                   of the artifacts themselves.
               save_dir:       Absolute path to a directory. In this directory,
                                   a subdirectory with name save_name will be 
                                   created.
           
           Returns:
               None. 

           Raises:
               None.
        """

        match = False
        for name, already_simulated_plan, dir_ in self.commitment_phase_metadata[-1].simulated_tool_pass_plans:
            if tp.compare_tool_pass_plans(tool_pass_plan, already_simulated_plan):
                # Copy the results of the already-simulated tool pass plan
                #     into a new directory.
                shutil.copytree(dir_, save_dir)
                
                # Rename the .cae file in this new directory.
                old_path = os.path.join(dir_, name + ".cae") 
                new_path = os.path.join(save_dir, save_name + ".cae")
                shutil.move(old_path, new_path)

                # Rename the .jnl file in this new directory.
                old_path = os.path.join(dir_, name + ".jnl")
                new_path = os.path.join(save_dir, save_name + ".jnl")
                shutil.move(old_path, new_path)

                match = True

                num_commits = len(self.commitment_phase_metadata)
                dp("For commit " + str(num_commits) + ", the tool pass plan was \
                    already simulated so no additional simulation was needed!")
                break

        # Simulate the plan if necessary.
        if not match:
            self.sim_potential_tool_passes(tool_pass_plan, save_name, save_dir)



    def _start_new_commitment_phase(self, save_name: str, save_dir: str) -> None:
        """Starts a new commitment phase by propagating the results of the last
               commitment phase.
            
           Among other things, this means that the part geometry and the state
               of stress of the part, as they existed after the last tool
               pass in the last commitment phase, are propagated into an
               MDB which is the starting point for this commitment phase.

           Args:
               save_name: The name of the subdirectory and .cae file (MDB) created
                              by this function. The MDB contains the propagated
                              results of the last commitment phase. The
                              subdirectory contains the .cae file.
               save_dir:  Absolute path to a directory. The directory in which
                              this function creates a subdirectory. 
        
           Returns:
               None.

           Raises:
               RuntimeError: The subdirectory cannot be created because something 
                                already exists at the path.
        """
    
        last_commit_metadata = self.commitment_phase_metadata[-1]

        committed_tp_plan = last_commit_metadata.committed_tool_pass_plan
        path_committed_tp_plan = last_commit_metadata.committed_tool_pass_plan_path

        assert committed_tp_plan is not None
        assert path_committed_tp_plan is not None

        # Record the path of the .sim file. The .sim file may be used in the
        #     next commitment phase to set the initial stress state.
        num_jobs = len(committed_tp_plan.plan)
        sim_file_name = shim.STANDARD_JOB_PREFIX + str(num_jobs) + ".sim"
        sim_file_path = os.path.join(path_committed_tp_plan, sim_file_name) 

        # Create the new subdirectory if the target path is unoccuppied. 
        new_subdir_path = os.path.join(save_dir, save_name)
        if not pathlib.Path(new_subdir_path).exists():
            os.mkdir(new_subdir_path)
        else:
            raise RuntimeError("The subdirectory already exists!")

        # Create the MDB for the next commitment phase and the metadata that
        #     accompanies it.
        mdb = shim.create_mdb(save_name, new_subdir_path) 
        new_mdb_full_path = os.path.join(save_dir, save_name)
        mdb_metadata = abq_md.AbaqusMdbMetadata(new_mdb_full_path) 

        # Generate the names for the stuff in the MDB.
        names = naming.ModelNames(mdb_metadata, True)

        # Create a part in the pre-existing lone model which comes with a new
        #     MDB.
        odb_name = shim.STANDARD_JOB_PREFIX + str(num_jobs) + ".odb"
        odb_path = os.path.join(path_committed_tp_plan, odb_name)
        shim.create_part_from_odb(names.pre_tool_pass_part_name, names.new_model_name, odb_path, mdb_metadata, mdb)

        # Map the orphan mesh to a part geometry.
        sim.orphan_mesh_to_geometry(names.pre_tool_pass_part_name, names.new_model_name, mdb)

        # Save the MDB so it is visible in the file system.
        shim.save_mdb(mdb)

        # Create the material based on the material of the very first part in
        #     the machining process. 
        first_commit_phase_metadata = self.commitment_phase_metadata[0]
        very_first_part = first_commit_phase_metadata.init_part
        material = very_first_part.material

        # The starting point for the next commitment phase.
        abaqus_part = part.AbaqusDefinedPart(save_name, new_mdb_full_path, material)

        # The initial state of the next commitment phase depends on the result
        #     of the simulation. 
        new_commit_metadata = md.CommittedToolPassPlanMetadata(abaqus_part, abaqus_part.path_to_mdb, self.boundary_conditions)
        self.commitment_phase_metadata.append(new_commit_metadata)

        # Record the stress state after the last tool pass in the previous commit.
        self.commitment_phase_metadata[-1].path_last_commit_sim_file = sim_file_path



    def _is_new_commitment_phase(self) -> bool:
        """Checks if the last commitment phase contains committed tool passes.
               This function can be used to check if a new commitment phase has
               started."""

        last_commitment_phase_metadata = self.commitment_phase_metadata[-1]
        if last_commitment_phase_metadata.committed_tool_pass_plan is not None:
            return True
        return False






