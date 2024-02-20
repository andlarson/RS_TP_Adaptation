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
import copy

import src.core.part.part as part
import src.core.simulation.simulation as sim
import src.core.boundary_conditions.boundary_conditions as bc
import src.core.metadata.metadata as md 
import src.core.metadata.abaqus_metadata as abq_md
import src.core.metadata.naming as naming
import src.core.tool_pass.tool_pass as tp
import src.core.abaqus.abaqus_shim as shim
import src.core.material_properties.material_properties as mp
import src.core.residual_stress.residual_stress as rs

from src.util.debug import *


STRESS_RECOVERY_MDB_NAME = "recover_residual_stresses.cae"
TARGET_GEOMETRY_MDB_NAME = "target_geometry.cae"

class MachiningProcess:

    def __init__(self, init_part: part.InitialPart, boundary_conditions: list[bc.BC]) -> None:
        """Defines the machining process for a single part.
   
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

        first_tp_metadata = md.CommitmentPhaseMetadata(init_part, init_part.path_to_mdb, boundary_conditions)
        first_tp_metadata.first_commitment_phase = True
        self.commitment_phase_metadata.append(first_tp_metadata)



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
               None.
        """

        # If the user did not supply a stress profile estimate for this commitment 
        #    phase, record as such.
        if len(self.stress_profile_estimates) < len(self.commitment_phase_metadata):
            self.stress_profile_estimates.append(None)

        committment_phase_md = self.commitment_phase_metadata[-1]
        committment_phase_md.committed_tpp = (save_name, tool_pass_plan, save_dir)
        committment_phase_md.committed_tpp_mdb_metadata = abq_md.AbaqusMdbMetadata(committment_phase_md.init_part.path_to_mdb) 
        
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
               None.
        """

        committment_phase_md = self.commitment_phase_metadata[-1]
        new_dir_path = os.path.join(save_dir, save_name)
        committment_phase_md.potential_tpps.append((save_name, tool_pass_plan, new_dir_path))
        mdb_metadata = abq_md.AbaqusMdbMetadata(committment_phase_md.init_part.path_to_mdb)
        committment_phase_md.potential_tpp_mdb_metadata.append(mdb_metadata)

        self._simulate_tpp(tool_pass_plan, save_name, save_dir)



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



    def estimate_stress_via_last_tool_pass(self, path: str) -> rs.ConstantResidualStressField:
        """Estimates the residual stress tensor field which existed in the region
               of material removal due to the last committed tool pass.
           
           This function should only be called if the last committed tool pass
               plan contains exactly one committed tool pass.

           Args:
               path: Absolute path to a directory. Since some simulations are
                         necessary to estimate the stress, this directory is
                         where the simulation results will be placed.

           Returns:
               The constant residual stress field in the region of material
                   removal.

           Raises:
               None.
        """

        if len(self.commitment_phase_metadata[-1].committed_tpp[2]) != 1:
            raise RuntimeError("The committed tool pass plan doesn't contain "
                               "exactly one committed tool pass.")
        tool_pass = self.commitment_phase_metadata[-1].committed_tpp[1].plan[0]

        # TODO: Right now, the MDBs are copied along with their metadata. None of
        #     this is recorded in the commitment phase metadata, so it cannot be 
        #     recovered when this function returns.

        # The post-committed tool pass geometry always comes from the current
        #     commitment phase.
        cur_commit_phase_md = self.commitment_phase_metadata[-1]
        stress_recovery_mdb_path = os.path.join(path, STRESS_RECOVERY_MDB_NAME)
        shutil.copyfile(cur_commit_phase_md.real_world_part.path_to_mdb, stress_recovery_mdb_path)
        stress_recovery_mdb = shim.use_mdb(stress_recovery_mdb_path)
        stress_recovery_mdb_md = copy.deepcopy(cur_commit_phase_md.real_world_part_mdb_md)
        
        # In the first commitment phase, the target geometry comes from the
        #     initial part. In all other commitment phases, the target geometry
        #     comes from the real world data of the previous commitment
        #     phase.
        target_geometry_mdb_path = os.path.join(path, TARGET_GEOMETRY_MDB_NAME)
        if len(self.commitment_phase_metadata) == 1:
            shutil.copyfile(cur_commit_phase_md.init_part.path_to_mdb, target_geometry_mdb_path) 
            target_geometry_mdb_md = abq_md.AbaqusMdbMetadata(target_geometry_mdb_path)
        else:
            prev_commit_phase_md = self.commitment_phase_metadata[-2]
            shutil.copyfile(prev_commit_phase_md.real_world_part.path_to_mdb, target_geometry_mdb_path)
            target_geometry_mdb_md = copy.deepcopy(prev_commit_phase_md.real_world_part_mdb_md)
        target_geometry_mdb = shim.use_mdb(target_geometry_mdb_path)
        
        return sim.estimate_residual_stresses(stress_recovery_mdb, stress_recovery_mdb_md,
                                              target_geometry_mdb, target_geometry_mdb_md,
                                              tool_pass, cur_commit_phase_md, path)



    def add_real_world_machining_data(self, part_from_real_life: part.MinimalPart) -> None:
        """Supplies real-world machining data for the LAST commitment phase. Should
               not be called in the first commitment phase and should not be called
               twice for the same commitment phase.

           The part passed to this function should be the part which resulted
               from a scan (in real life) of the part which resulted from the
               last committed tool pass plan.
           
           Args:
               part_from_real_life: The result of performing the last committed 
                                        tool pass plan in real life. The geometric
                                        model should be constructed based on
                                        data collected from real life (such as
                                        an in-maching scan).
        
           Returns:
               None.
        
           Raises:
               None.
        """

        if len(self.commitment_phase_metadata) <= 1:
            raise RuntimeError("Real world machining data should only be passed"
                               " after the first commitment phase.")

        last_commitment_phase_md = self.commitment_phase_metadata[-2]

        if hasattr(last_commitment_phase_md, "real_world_part"):
            raise RuntimeError("The last commitment phase already had real world"
                               " data supplied for it!")

        last_commitment_phase_md.real_world_part = part_from_real_life

        # This MDB needs metadata and it is in the default initial state. 
        mdb_md = abq_md.AbaqusMdbMetadata(part_from_real_life.path_to_mdb) 
        last_commitment_phase_md.real_world_part_mdb_md = mdb_md



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
        tpps_done = self.commitment_phase_metadata[-1].potential_tpps
        for i, simulated_tpp in enumerate(tpps_done):
            name = simulated_tpp[0]
            plan = simulated_tpp[1]
            dir_ = simulated_tpp[2]
            if tp.compare_tool_pass_plans(tool_pass_plan, plan):

                new_dir_name = os.path.join(save_dir, save_name)

                # Copy the results of the already-simulated tool pass plan
                #     into a new directory.
                shutil.copytree(dir_, new_dir_name)
                
                # Rename the .cae file in this new directory.
                old_mdb_path = os.path.join(new_dir_name, name + ".cae") 
                new_mdb_path = os.path.join(new_dir_name, save_name + ".cae")
                shutil.move(old_mdb_path, new_mdb_path)

                # Rename the .jnl file in this new directory.
                old_jnl_path = os.path.join(new_dir_name, name + ".jnl")
                new_jnl_path = os.path.join(new_dir_name, save_name + ".jnl")
                shutil.move(old_jnl_path, new_jnl_path)

                match = True

                num_commits = len(self.commitment_phase_metadata)

                dump_banner("COMMITTED TOOL PASS PLAN WAS ALREADY SIMULATED")
                dp("")
                dp("For commit " + str(num_commits) + ", the tool pass plan was "
                   "already simulated in the commitment phase so no additional "
                   "simulation was needed!")
                dp("")
                dump_banner_end()

                # Even when the tool pass plan has already been simulated, it's
                #     still necessary to record the metadata about the new
                #     MDB that was just copied. 
                commitment_phase_md = self.commitment_phase_metadata[-1]
                mdb_metadata = copy.deepcopy(commitment_phase_md.potential_tpp_mdb_metadata[i])
                mdb_metadata.path_to_mdb = new_mdb_path
                commitment_phase_md.committed_tpp_mdb_metadata = mdb_metadata
                commitment_phase_md.committed_tpp = (name, copy.deepcopy(plan), dir_)
                break

        # Simulate the plan if necessary.
        if not match:
            self._simulate_tpp(tool_pass_plan, save_name, save_dir)



    def _simulate_tpp(self, tool_pass_plan: tp.ToolPassPlan, save_name: str, 
                      save_dir: str) -> None:
        """Simulates a tool pass plan. Does not make assumptions about the tool
               pass plan being a committed tool pass plan or a potential tool
               pass plan.
            
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

        commit_metadata = self.commitment_phase_metadata[-1]

        # If there is an estimated stress profile for this commitment phase, use it.
        if len(self.commitment_phase_metadata) == len(self.stress_profile_estimates):
            sim.sim_tool_pass_plan(tool_pass_plan, save_name, save_dir, 
                                   commit_metadata, self.stress_profile_estimates[-1])
        else:
            sim.sim_tool_pass_plan(tool_pass_plan, save_name, save_dir, commit_metadata)



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
    
        cur_commitment_phase_md = self.commitment_phase_metadata[-1]
        committed_tpp_mdb_md = cur_commitment_phase_md.committed_tpp_mdb_metadata

        committed_tpp = cur_commitment_phase_md.committed_tpp
        path_committed_tpp = cur_commitment_phase_md.committed_tpp_mdb_metadata.mdb_dir()

        assert committed_tpp is not None
        assert path_committed_tpp is not None

        # Record the path of the .sim file. The .sim file may be used in the
        #     next commitment phase to set the initial stress state.
        sim_file_name = naming.last_sim_file_name(committed_tpp_mdb_md)
        sim_file_path = os.path.join(path_committed_tpp, sim_file_name) 

        # Create the new subdirectory if the target path is unoccuppied. 
        new_subdir_path = os.path.join(save_dir, save_name)
        if not pathlib.Path(new_subdir_path).exists():
            os.mkdir(new_subdir_path)
        else:
            raise RuntimeError("The subdirectory already exists!")
        
        # Retrieve the path to the ODB which contains the result of the committed
        #     tool pass plan in the last committment phase.
        odb_name = naming.last_odb_file_name(committed_tpp_mdb_md)
        odb_path = os.path.join(path_committed_tpp, odb_name)

        # Use the ODB to create a new MDB with the deformed geometry in it.
        sim.create_mdb_from_odb(save_name, new_subdir_path, odb_path)

        # The material comes from the material used in the very first commitment
        #     phase.
        first_commit_phase_metadata = self.commitment_phase_metadata[0]
        very_first_part = first_commit_phase_metadata.init_part
        material = very_first_part.material

        # The starting point for the next commitment phase.
        new_mdb_path = os.path.join(new_subdir_path, save_name)
        abaqus_part = part.InitialPart(save_name, new_mdb_path, material)

        # The initial state of the next commitment phase depends on the result
        #     of the simulation. 
        new_phase_md = md.CommitmentPhaseMetadata(abaqus_part, abaqus_part.path_to_mdb, self.boundary_conditions)
        self.commitment_phase_metadata.append(new_phase_md)

        # Record the stress state which resulted from the committed tool pass 
        #     plan. 
        self.commitment_phase_metadata[-1].path_last_commit_sim_file = sim_file_path



    def _is_new_commitment_phase(self) -> bool:
        """Checks if the last commitment phase contains committed tool passes.
               This function can be used to check if a new commitment phase has
               started."""

        last_commitment_phase_md = self.commitment_phase_metadata[-1]
        if hasattr(last_commitment_phase_md, "committed_tpp"):
            return True
        return False






