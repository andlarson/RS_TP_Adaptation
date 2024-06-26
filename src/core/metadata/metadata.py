"""
This module contains functionality which can be used to keep track of and retrieve 
    metadata associated with abstractions which don't natively exist in Abaqus. 
    Abaqus metadata is embedded in this metadata.
"""

import src.core.part.part as part
import src.core.boundary_conditions.boundary_conditions as bc
import src.core.metadata.abaqus_metadata as abq_md
import src.core.tool_pass.tool_pass as tp
import src.core.real_world_data.real_world_data as rwd
import src.core.material_properties.material_properties as mp


class CommitmentPhaseMetadata:

    def __init__(self, init_part: part.MinimalPart) -> None:
        """Creates a data structure associated with a single commitment phase.

           The system flow is assumed to be something like:
               1) User simulates a bunch of potential tool pass plans.
               2) User examines the results of the tool pass plans they simulated.
                      They decide to perform one of the tool pass plans in real
                      life, committing to it.
               3) User scans the workpiece after the tool pass plan is done.
               4) User recovers the residual stresses which caused the workpiece
                      to deform during/after cutting. They use this residual
                      stress information to help them decide what the next tool
                      pass plan should be.
               5) Repeat steps 1-4 until done machining.

           Thus, there are natural phases to the procedure. A phase corresponds
               to all the stuff done between two tool pass plans in real life. 
               This data structure records information about what occurred during 
               a phase.
           
           Args:
               init_part: The part at the beginning of this phase. Note that 
                              this part may have had some material removed 
                              from it already.

           Returns:
               None.

           Raises:
               None.
        """
        
        # ----- Miscellaneous / Startup -----

        # Flag indicating if this commitment phase is the very first one.
        self.first_commitment_phase: bool = False
        
        # The state of the part as it exists at the beginning of this commitment 
        #     phase.
        self.init_part: part.MinimalPart = init_part


        # ----- User-Provided Initial Stress State -----
        
        # Absolute path to object file defining a state of stress.
        # The state of stress used as initial state of stress for all simultions
        #     in this commitment phase.
        self.init_stress: str | None = None


        # ----- Potential Tool Pass Plans -----

        # Each potential tool pass plan is simulated in an MDB. This stores
        #     metadata about each of these MDBs.
        self.potential_tpp_mdb_metadata: list[abq_md.AbaqusMdbMetadata] = []

        # Optimization: To speed up committing a tool pass plan, a list of all
        #     the potential tool pass plans which have been simulated is 
        #     maintained. That way, it's possible to check if the tool pass plan 
        #     has already been simulated before commitment is done. If so, no 
        #     need to re-simulate. 
        # Name, tool pass plan, and absolute path to directory of simulation
        #     results.
        self.potential_tpps: list[tuple[str, tp.ToolPassPlan, str]] = [] 


        # ----- Committed Tool Pass Plan -----
        
        # The committed tool pass plan must also be simulated in an MDB.
        self.committed_tpp_mdb_metadata: abq_md.AbaqusMdbMetadata | None = None

        # The tool pass plan which was committed to for this commitment phase.
        # The name, tool pass plan, and absolute path to directory of simulation
        #     results.
        self.committed_tpp: tuple[str, tp.ToolPassPlan, str] | None = None

        # Absolute path to .sim file resulting from the committed tool pass plan.
        #     This .sim file may be used as the starting stress state for simulations
        #     in the next commitment phase if the user doesn't provide an estimated
        #     state of stress.
        self.sim_path: str | None = None


        # ----- Estimating Stress -----

        # The process of estimating the stress in a region of material removal
        #     requires running simulations, and therefore requires an MDB.
        self.stress_estimate_mdb_metadata: abq_md.AbaqusMdbMetadata | None = None


        # ----- Real World Data -----
        
        # Absolute path to .stl or .odb file.
        # The real world part scan data corresponding to the result of the
        #     committed tool pass plan.
        self.real_world_data: str | None = None


        # ----- User-Provided Stress State For Next Commitment Phase -----
        
        # Absolute path to object file defining a state of stress.
        # This is state of stress for the NEXT commitment phase.
        # When the next commitment phase starts, this is copied into the metadata
        #     for that new commitment phase. 
        # Yes, this design is inelegant. It results from the fact that metadata
        #     for the next commitment phase is only created when the next commitment
        #     phase actually starts, and the next commitment phase only starts
        #     when a tool pass simulation happens.
        self.next_phase_init_stress: str | None = None



        








