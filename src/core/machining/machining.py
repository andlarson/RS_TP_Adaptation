"""
This module contains top-level, user-facing functionality for:
    1) Performing simulations of potential tool passes.
    2) Recovering information about the residual stress field.

The MachiningProcess object is the top-level object that a user of this library
    should create. This object's methods implement the functionality listed
    above.
"""

from __future__ import annotations
from typing import Any, TYPE_CHECKING
import shutil
import os
import pathlib
import copy

# Imports only used for static analysis / type checking.
# This prevents cyclic imports, which are sometimes necessary for type annotations,
#     from causing runtime errors. 
if TYPE_CHECKING:
    import src.core.boundary_conditions.boundary_conditions as bc
    import src.core.material_properties.material_properties as mp
    import src.core.residual_stress.residual_stress as rs

import src.core.part.part as part
import src.core.simulation.simulation as sim
import src.core.metadata.metadata as md 
import src.core.metadata.abaqus_metadata as abq_md
import src.core.metadata.naming as naming
import src.core.tool_pass.tool_pass as tp
import src.core.abaqus.abaqus_shim as shim

from src.util.debug import *



class MachiningProcess:

    def __init__(self, init_part: str, boundary_conditions: list[bc.BC],
                 material: mp.ElasticMaterial, save_path: str) -> None:
        """Defines the machining process for a single part.
   
           Args:
               init_part:           Absolute path to .stl file defining the geometry
                                        of the blank.
               boundary_conditions: The boundary conditions which should be used
                                        for all simulations in this machining
                                        process. Usually, these reflect the
                                        clamping conditions of the part.
                                        Without any boundary conditions, a
                                        finite element simulation is
                                        unconstrained. Thus, there should always
                                        be at least one boundary condition.
               material:            The material that the part is made out of.
               save_path:           Absolute path to desired save location of
                                        .cae file. This path should end with
                                        .cae.
           
           Returns:
               None.
           
           Raises:
               RuntimeError: At least one boundary condition must be specified.
        """
     
        if len(boundary_conditions) == 0:
            raise RuntimeError("At least one boundary condition must be supplied.")

        self.commitment_phase_metadata: list[md.CommitmentPhaseMetadata] = []
        self.invariants: _MachiningInvariants = _MachiningInvariants(boundary_conditions, material, init_part)
            
        mdb, mdb_md = shim.create_mdb_from_mesh(init_part, save_path)
        shim.save_mdb(mdb)
        shim.close_mdb(mdb)

        min_part = part.MinimalPart(mdb_md.path_to_mdb)

        first_commit_md = md.CommitmentPhaseMetadata(min_part)
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
        
        # The tool pass plan was simulated, so record the metadata. 
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
        
        # TODO: Technically we can do a stress estimation even if no real world
        #     data has been supplied. In this case, the stress estimation would
        #     be based on the simulation results. For now, we won't allow that.

        if cur_commit_phase_md.real_world_data is None:
            raise RuntimeError("The real world data resulting from the most recent \
                                committed tool pass plan has not been supplied.")

        if not cur_commit_phase_md.first_commitment_phase:
            if self.commitment_phase_metadata[-2].real_world_data is None:
                raise RuntimeError("The real world data resulting from the last \
                                    committed tool pass plan was not supplied. In \
                                    order to do a stress estimation, the real \
                                    world data for the last two committed tool \
                                    pass plans must be supplied.")

        tool_pass = self.commitment_phase_metadata[-1].committed_tpp[1].plan[0]
        
        if cur_commit_phase_md.first_commitment_phase:
            target_geometry_data = self.invariants.blank
        else:
            prev_commit_phase_md = self.commitment_phase_metadata[-2]
            target_geometry_data = prev_commit_phase_md.real_world_data
        
        deformed_geometry_data = cur_commit_phase_md.real_world_data

        return sim.estimate_residual_stresses(deformed_geometry_data, target_geometry_data,
                                              tool_pass, self.invariants, path)



    def add_real_world_machining_data(self, real_world_data: str) -> None:
        """Supplies real-world machining data for the LAST committed tool pass
               plan. 
           
           Should not be called before the first tool pass plan is committed. Should
               not be called when the committed tool pass plan contains many
               non-contiguous or very large tool passes. Must be called between
               committing a tool pass plan and simulating another tool pass plan.

           The data passed to this function should be the data which resulted
               from a scan of the part which resulted from the last committed 
               tool pass plan.
           
           Args:
               real_world_data: Absolute path to .stl or .odb file containing the
                                    real world data. Note that .stl is the format
                                    which will probably be generated by a true
                                    in-machine scan. The option to pass a .odb file
                                    exists because we want to provide the option
                                    to do everything in simulation.
        
           Returns:
               None.
        
           Raises:
               None.
        """

        cur_commitment_phase_md = self.commitment_phase_metadata[-1]

        if cur_commitment_phase_md.committed_tpp is None:
            raise RuntimeError("Real world machining data should only be passed \
                                after a tool pass plan is committed!")

        if cur_commitment_phase_md.real_world_data is not None:
            raise RuntimeError("Real world data was already supplied for the \
                                most recent committed tool pass plan!")

        cur_commitment_phase_md.real_world_data = real_world_data 



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
            sim.sim_tool_pass_plan(tool_pass_plan, save_name, save_dir, self.invariants, 
                                   mdb_md, mdb, cur_commit_md.init_stress)
        else:
            # Otherwise, use the stress state which resulted from the last committed
            #     tool pass plan being simulated.

            # Warn the user!
            dump_banner("NO STRESS ESTIMATE SUPPLIED")
            dp("The stress field which resulted from the simulation of the"\
               " last committed tool pass plan is being used.")
            dump_banner_end()

            prev_commit_md = self.commitment_phase_metadata[-2]
            sim.sim_tool_pass_plan(tool_pass_plan, save_name, save_dir, self.invariants, 
                                   mdb_md, mdb, prev_commit_md.sim_path)


    def _start_new_commitment_phase(self, save_name: str, save_dir: str) -> None:
        """Starts a new commitment phase by propagating the results of the last
               commitment phase.
           
           Exactly what is propagated depends on what the user supplied. For
               example, if the user supplied a stress estimate, then that stress
               estimate is used in the new commitment phase. Otherwise, the
               field produced by the simulation is used in the new commitment
               phase.
           
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
        old_commit_phase_md.sim_path = sim_file_path

        # Create the new subdirectory if the target path is unoccuppied. 
        new_subdir_path = os.path.join(save_dir, save_name)
        if not pathlib.Path(new_subdir_path).exists():
            os.mkdir(new_subdir_path)
        else:
            raise RuntimeError("The subdirectory already exists!")
        
        # The starting point for the next commitment phase.
        new_mdb_path = os.path.join(new_subdir_path, save_name)

        # If the user supplied real world data for the result of the committed
        #     tool pass plan in this commitment phase, then that real world
        #     data is the starting point for the next commitment phase. If the
        #     user didn't supply real world data, then the result (i.e. .odb)
        #     of simulating the committed tool pass plan is the starting point
        #     for the commitment phase.
        if old_commit_phase_md.real_world_data is None:
            # Retrieve the path to the ODB which contains the result of the committed
            #     tool pass plan. 
            odb_name = naming.last_odb_file_name(committed_tpp_mdb_md)
            odb_path = os.path.join(path_committed_tpp, odb_name)

            # Use the ODB to create a new MDB with the deformed geometry in it.
            mdb, mdb_md = shim.create_mdb_from_mesh(odb_path, new_mdb_path)

            # Warn the user!
            dump_banner("No real world data was supplied for the result"\
                        " of the most recent committed tool pass plan. Thus, the"\
                        " result from simulation is being used!")
            dump_banner_end()
        else: 
            mdb, mdb_md = shim.create_mdb_from_mesh(old_commit_phase_md.real_world_data, new_mdb_path)

        shim.save_mdb(mdb)

        abaqus_part = part.MinimalPart(new_mdb_path)
        new_phase_md = md.CommitmentPhaseMetadata(abaqus_part)
        self.commitment_phase_metadata.append(new_phase_md)
        
        new_commit_phase_md = self.commitment_phase_metadata[-1]
        
        if old_commit_phase_md.next_phase_init_stress is not None:
            # Propagate the user-defined stress profile.
            new_commit_phase_md.init_stress = old_commit_phase_md.next_phase_init_stress



    def _time_for_new_commitment_phase(self) -> bool:
        """Determines if a new commitment phase needs to be started."""

        cur_commitment_phase_md = self.commitment_phase_metadata[-1]
        if cur_commitment_phase_md.committed_tpp is not None:
            return True
        return False



"""
This class is basically a data structure which holds data that does not change
    during the machining process.
For example, it is assumed that the boundary conditions (aka the clamping
    conditions) do not change during the machining process.
"""
class _MachiningInvariants:
    
    def __init__(self, boundary_conditions: list[bc.BC], 
                 material: mp.ElasticMaterial, blank: str):
        """Populates the data structure.

           Args:
               boundary_conditions: The boundary conditions valid for the whole
                                        machining process.
               material:            The material that the part is made out of.
               blank:               Absolute path to .stl file defining the
                                        geometry of the blank.

           Returns:
               None.

           Raises:
               None.
        """

        self.boundary_conditions = boundary_conditions
        self.material = material
        self.blank = blank
