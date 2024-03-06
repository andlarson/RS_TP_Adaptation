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
from typing import Any

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


DEFORMED_GEOMETRY_MDB_NAME = "recover_residual_stresses.cae"
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
        self.boundary_conditions = boundary_conditions 

        first_commit_md = md.CommitmentPhaseMetadata(init_part, init_part.path_to_mdb, boundary_conditions)
        first_commit_md.first_commitment_phase = True

        self.commitment_phase_metadata.append(first_commit_md)



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

        if self._time_for_new_commitment_phase():
            self._start_new_commitment_phase(save_name, save_dir) 
        else:
            # To make things symmetric, it's necessary to create a directory here
            #     like one is created when starting a new commitment phase.
            new_subdir_path = os.path.join(save_dir, save_name)
            os.mkdir(new_subdir_path)
        
        done, mdb_md = self._lazy_tpp_sim(tool_pass_plan, save_name, save_dir)
        if done:
            # The tool pass plan was already conducted, so its metadata
            #     can just be copied.
            commit_phase_md = self.commitment_phase_metadata[-1]
            commit_phase_md.committed_tpp = (save_name, tool_pass_plan, save_dir)
            commit_phase_md.committed_tpp_mdb_metadata = mdb_md 
        else:
            # This is inelegant but necessary. The simulation process causes the
            #     MDB to be modified, and we want to keep track of the modifications.
            #     Thus, the metadata data structure must be created in advance of
            #     the modifications.
            commit_phase_md = self.commitment_phase_metadata[-1]
            path_to_mdb = commit_phase_md.init_part.path_to_mdb
            mdb_md = abq_md.AbaqusMdbMetadata(path_to_mdb)

            mdb = shim.use_mdb(path_to_mdb)
            self._simulate_tpp(tool_pass_plan, save_name, save_dir, mdb_md, mdb)            
            shim.close_mdb(mdb)
            
            # The tool pass plan was simulated, so record the metadata. 
            commit_phase_md.committed_tpp = (save_name, tool_pass_plan, save_dir)
            commit_phase_md.committed_tpp_mdb_metadata = mdb_md 



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
        
        if self._time_for_new_commitment_phase():
            self._start_new_commitment_phase(save_name, save_dir) 
        else:
            # To make things symmetric, it's necessary to create a directory here
            #     like one is created when starting a new commitment phase.
            new_subdir_path = os.path.join(save_dir, save_name)
            os.mkdir(new_subdir_path)
        
        # This is inelegant but necessary. The simualtion process causes the
        #     MDB to be modified, and we want to keep track of the modifications.
        #     Thus, the metadata data structure must be created in advance of
        #     the modifications.
        commit_phase_md = self.commitment_phase_metadata[-1]
        path_to_mdb = commit_phase_md.init_part.path_to_mdb
        mdb_md = abq_md.AbaqusMdbMetadata(path_to_mdb)

        mdb = shim.use_mdb(path_to_mdb)
        self._simulate_tpp(tool_pass_plan, save_name, save_dir, mdb_md, mdb)            
        shim.close_mdb(mdb)
        
        # The tool pass plan was simualted, so record the metadata. 
        new_dir_path = os.path.join(save_dir, save_name)
        commit_phase_md.potential_tpps.append((save_name, tool_pass_plan, new_dir_path))
        commit_phase_md.potential_tpp_mdb_metadata.append(mdb_md)



    def use_stress_profile(self, path: str):
        """Records a stress profile to be used as the initial stress state for
               all the simulations until (and including) the next tool pass
               plan to be committed.
           
           If a stress profile is provided, it should only be provided in some
               circumstances:
               1) Before doing any simulations.
               2) Right after committing to a tool pass plan, and before
                      simulating another pass plan.
           
           If no stress profile is provided after committing to a tool pass
               plan, then all simulations until, and including, the next 
               committed tool pass plan will use the stress profile which 
               Abaqus produced from the committed tool pass plan. Since the
               user is required to specify an initial stress profile, this makes
               it possible to only ever specify the very initial stress state.

           Args:
               path: Path to the object file produced by compiling the stress
                         subroutine which defines the desired stress profile.
              
           Return:
               None.

           Raises:
               None.
        """

        cur_commit_phase_md = self.commitment_phase_metadata[-1]

        if cur_commit_phase_md.first_commitment_phase and cur_commit_phase_md.committed_tpp is None:
            cur_commit_phase_md.init_stress = path
        else:
            if cur_commit_phase_md.committed_tpp is None:
                raise RuntimeError("A stress profile should only be passed after a \
                                    tool pass plan has been committed.")

            if cur_commit_phase_md.next_phase_init_stress is not None:
                raise RuntimeError("A stress profile was already provided for the \
                                    last tool pass plan to be committed.")
        
            cur_commit_phase_md.next_phase_init_stress = path



    def estimate_stress(self, path: str) -> rs.ConstantResidualStressField:
        """Estimates the residual stress tensor field which existed in the region
               of material removal due to the last committed tool pass plan.
           
           This function should only be called after a tool pass plan has been
               committed and the real world data associated with the committed 
               tool pass plan has been collected and passed to this library. It
               should be called before simulating any additional tool pass plans.
           
           This function should only be called if the last committed tool pass
               plan contains exactly one tool pass. The underlying
               technique for the estimation is ineffective if the tool pass plan
               contains more than a single committed tool pass.

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
        
        cur_commit_phase_md = self.commitment_phase_metadata[-1]
        
        if cur_commit_phase_md.committed_tpp is None:
            raise RuntimeError("A tool pass plan needs to be committed.")

        if len(cur_commit_phase_md.committed_tpp[1]) != 1:
            raise RuntimeError("The committed tool pass plan doesn't contain \
                                exactly one committed tool pass.")

        if cur_commit_phase_md.real_world_part is None:
            raise RuntimeError("The real world data associated with the committed \
                                tool pass plan has not been passed.")

        tool_pass = self.commitment_phase_metadata[-1].committed_tpp[1].plan[0]

        # TODO: Right now, the MDBs are copied along with their metadata. None of
        #     this is recorded in the commitment phase metadata, so it is not
        #     tracked at all.
        
        # The real world data for this commitment phase is the deformed geometry.
        deformed_geometry_mdb_path = os.path.join(path, DEFORMED_GEOMETRY_MDB_NAME)
        to_copy = cur_commit_phase_md.real_world_part.path_to_mdb
        shim.copy_mdb(to_copy, deformed_geometry_mdb_path)

        # Reuse the metadata and update it. 
        deformed_geometry_mdb_md = copy.deepcopy(cur_commit_phase_md.real_world_part_mdb_md)
        deformed_geometry_mdb_md.path_to_mdb = deformed_geometry_mdb_path
        
        # In the first commitment phase, the target geometry comes from the
        #     initial part. In all other commitment phases, the target geometry
        #     comes from the real world data of the previous commitment
        #     phase.
        target_geometry_mdb_path = os.path.join(path, TARGET_GEOMETRY_MDB_NAME)
        if len(self.commitment_phase_metadata) == 1:
            # The initial workpiece geometry is the target geometry. 
            to_copy = cur_commit_phase_md.init_part.path_to_mdb
            shim.copy_mdb(to_copy, target_geometry_mdb_path)

            # The MDB containing the initial part has no metadata associated
            #     with it.
            target_geometry_mdb_md = abq_md.AbaqusMdbMetadata(target_geometry_mdb_path)
        else:
            # The real world data from the last commitment phase is the target
            #     geometry.
            prev_commit_phase_md = self.commitment_phase_metadata[-2]
            to_copy = prev_commit_phase_md.real_world_part.path_to_mdb
            shim.copy_mdb(to_copy, target_geometry_mdb_path)
            
            # Reuse the metadata and update it.
            target_geometry_mdb_md = copy.deepcopy(prev_commit_phase_md.real_world_part_mdb_md)
            target_geometry_mdb_md.path_to_mdb = target_geometry_mdb_path
        
        return sim.estimate_residual_stresses(deformed_geometry_mdb_md, target_geometry_mdb_md,
                                              tool_pass, cur_commit_phase_md, path)



    def add_real_world_machining_data(self, part_from_real_life: part.MinimalPart) -> None:
        """Supplies real-world machining data for the LAST committed tool pass
               plan. 
           
           Should not be called before the first tool pass plan is committed. Should
               not be called when the committed tool pass plan contains many
               non-contiguous or very large tool passes. Must be called between
               committing a tool pass plan and simulating another tool pass plan.

           The part passed to this function should be the part which resulted
               from a scan (in real life) of the part which resulted from the
               last committed tool pass plan.
           
           Args:
               part_from_real_life: The result of performing the last committed 
                                        tool pass plan in real life. The geometric
                                        model should be constructed based on
                                        data collected from real life (such as
                                        an in-machine scan).
        
           Returns:
               None.
        
           Raises:
               None.
        """

        cur_commitment_phase_md = self.commitment_phase_metadata[-1]

        if cur_commitment_phase_md.committed_tpp is None:
            raise RuntimeError("Real world machining data should only be passed \
                                after a tool pass plan is committed!")

        if cur_commitment_phase_md.real_world_part is not None:
            raise RuntimeError("Real world data was already supplied for the \
                                most recent committed tool pass plan!")

        cur_commitment_phase_md.real_world_part = part_from_real_life

        # This MDB needs metadata and it is in the default initial state. 
        mdb_md = abq_md.AbaqusMdbMetadata(part_from_real_life.path_to_mdb) 
        cur_commitment_phase_md.real_world_part_mdb_md = mdb_md



    def _lazy_tpp_sim(self, tpp: tp.ToolPassPlan, save_name: str, save_dir: str
                     ) -> tuple[bool, None | abq_md.AbaqusMdbMetadata]:
        """Copies the results of a simulation if it has already been done in 
               this commitment phase. 
           
           Args:
               tpp:       The tool pass plan. 
               save_name: The name of the subdirectory in which the simulation 
                              artifacts (.odb file, .sim file, etc. and the 
                              final .cae file) will be placed. Also, the name 
                              of the artifacts themselves.
               save_dir:  Absolute path to a directory. In this directory,
                              a subdirectory with name save_name will be 
                              created if one does not already exist. If a directory
                              with save_name already exists, this function
                              copies into that subdirectory, overwriting files
                              as necessary.
           
           Returns:
               Boolean which indicates if an exact copy was found and, if one
                   was found, the metadata associated with the mdb that was
                   copied.

           Raises:
               None.
        """
        
        commitment_phase_md = self.commitment_phase_metadata[-1]
        tpps_done = self.commitment_phase_metadata[-1].potential_tpps
        for i, simulated_tpp in enumerate(tpps_done):
            name = simulated_tpp[0]
            plan = simulated_tpp[1]
            dir_ = simulated_tpp[2]
            if tp.compare_tool_pass_plans(tpp, plan):

                new_dir_name = os.path.join(save_dir, save_name)
                
                # Copy the already-produced results.
                shutil.copytree(dir_, new_dir_name, dirs_exist_ok=True)
                
                # Rename the .cae file. 
                old_mdb_path = os.path.join(new_dir_name, name + ".cae") 
                new_mdb_path = os.path.join(new_dir_name, save_name + ".cae")
                shutil.move(old_mdb_path, new_mdb_path)

                # Rename the .jnl file.
                old_jnl_path = os.path.join(new_dir_name, name + ".jnl")
                new_jnl_path = os.path.join(new_dir_name, save_name + ".jnl")
                shutil.move(old_jnl_path, new_jnl_path)

                num_commits = len(self.commitment_phase_metadata)
                dump_banner("COMMITTED TOOL PASS PLAN WAS ALREADY SIMULATED")
                dp("")
                dp("For commit " + str(num_commits) + ", the tool pass plan was "
                   "already simulated in the commitment phase so no additional "
                   "simulation was needed!")
                dp("")
                dump_banner_end()

                return True, copy.deepcopy(commitment_phase_md.potential_tpp_mdb_metadata[i])

        return False, None



    def _simulate_tpp(self, tool_pass_plan: tp.ToolPassPlan, save_name: str, 
                      save_dir: str, mdb_md: abq_md.AbaqusMdbMetadata,
                      mdb: Any) -> None:
        """Simulates a tool pass plan in the current commitment phase. Does 
               not make assumptions about the tool pass plan being a committed 
               tool pass plan or a potential tool pass plan.
            
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
               mdb:            Abaqus MDB object. The MDB in which to simulate
                                   the tool pass plan. This MDB need not be empty.
               mdb_md:         The metadata associated with the MDB.
        
           Returns:
               None. 
           
           Raises:
               RuntimeError: No stress profile was supplied for the first commitment
                                 phase. If no stress profile is specified for the
                                 first commitment phase, no deformation will occur
                                 when tool paths are simulated.
        """

        # For the first commitment phase, the user must supply a stress profile.
        cur_commit_md = self.commitment_phase_metadata[-1]
        if cur_commit_md.first_commitment_phase and cur_commit_md.init_stress is None:
            raise RuntimeError("A user-specified stress profile must be supplied \
                                before any simulations can run!")

        if cur_commit_md.init_stress is not None:
            # If there is an estimated stress profile for this commitment phase, use it.
            sim.sim_tool_pass_plan(tool_pass_plan, save_name, save_dir, cur_commit_md, 
                                   mdb_md, mdb, cur_commit_md.init_stress)
        else:
            sim.sim_tool_pass_plan(tool_pass_plan, save_name, save_dir, cur_commit_md,
                                   mdb_md, mdb)



    def _start_new_commitment_phase(self, save_name: str, save_dir: str) -> None:
        """Starts a new commitment phase by propagating the results of the last
               commitment phase.
            
           Among other things, this means that the part geometry and the state
               of stress of the part, as they existed after the last tool
               pass in the last commitment phase, are propagated into an
               MDB which is the starting point for this new commitment phase.

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
    
        old_commit_phase_md = self.commitment_phase_metadata[-1]
        committed_tpp_mdb_md = old_commit_phase_md.committed_tpp_mdb_metadata

        path_committed_tpp = old_commit_phase_md.committed_tpp_mdb_metadata.mdb_dir()

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
        #     tool pass plan in the current commitment phase. 
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
        abaqus_part = part.InitialPart(new_mdb_path, material)

        # The initial state of the next commitment phase depends on the result
        #     of the simulation. 
        new_phase_md = md.CommitmentPhaseMetadata(abaqus_part, abaqus_part.path_to_mdb, self.boundary_conditions)
        self.commitment_phase_metadata.append(new_phase_md)
        
        new_commit_phase_md = self.commitment_phase_metadata[-1]

        # Record the stress state which resulted from the committed tool pass 
        #     plan in this new commitment phase.
        new_commit_phase_md.path_last_commit_sim_file = sim_file_path

        # Propagate the user-defined stress profile.
        new_commit_phase_md.init_stress = old_commit_phase_md.next_phase_init_stress



    def _time_for_new_commitment_phase(self) -> bool:
        """Determines if a new commitment phase needs to be started."""

        cur_commitment_phase_md = self.commitment_phase_metadata[-1]
        if cur_commitment_phase_md.committed_tpp is not None:
            return True
        return False
