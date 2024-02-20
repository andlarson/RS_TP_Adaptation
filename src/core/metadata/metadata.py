"""
This module contains functionality which can be used to keep track of and retrieve 
    metadata associated with abstractions which don't natively exist in Abaqus. 
    Abaqus metadata is embedded in this metadata.
"""

import src.core.part.part as part
import src.core.boundary_conditions.boundary_conditions as bc
import src.core.metadata.abaqus_metadata as abq_md
import src.core.tool_pass.tool_pass as tp


class CommitmentPhaseMetadata:

    def __init__(self, init_part: part.InitialPart, path_to_mdb: str, 
                 BCs: list[bc.BC]) -> None:
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
           
           Also, note that there is an alternative system flow wherein everything
               is done in simulation and no machining is done in real life. However, 
               from the perspective of this data structure, nothing changes.
           
           Args:
               init_part:   The part at the beginning of this phase. Note that 
                                this part may have had some material removed 
                                from it already.
               path_to_mdb: The path to the MDB which contains the state of the 
                                part at the beginning of this commitment phase.
               BCs:         The boundary conditions active during this phase.

           Returns:
               None.

           Raises:
               None.
        """
        
        # ----- Miscellaneous / Startup -----

        # Flag indicating if this commitment phase is the very first one.
        self.first_commitment_phase: bool = False

        # The absolute path to the .sim file from the commitment phase which preceded
        #     this commitment phase.
        # This stress profile is sometimes used as the starting stress profile for
        #     this commitment phase.
        # This is really a convenience. It would not be hard to lookup the
        #     the .sim file which resulted from the committed tool pass plan in 
        #     the commitment phase which precedes this one.
        self.path_last_commit_sim_file: str
        
        # The state of the part as it exists at the beginning of this commitment 
        #     phase.
        self.init_part: part.InitialPart = init_part
        
        # The boundary conditions for this commitment phase.
        self.BCs: list[bc.BC] = BCs


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
        self.committed_tpp_mdb_metadata: abq_md.AbaqusMdbMetadata

        # The tool pass plan which was committed to for this commitment phase.
        # The name, tool pass plan, and absolute path to directory of simulation
        #     results.
        self.committed_tpp: tuple[str, tp.ToolPassPlan, str]


        # ----- Stress Estimation -----

        # The process of estimating the stress in a region of material removal
        #     requires running simulations, and therefore requires an MDB.
        self.stress_estimate_mdb_metadata: abq_md.AbaqusMdbMetadata


        # ----- Real World Data -----
       
        # A representation of the part that resulted from performing the committed
        #     tool pass plan.
        self.real_world_part: part.MinimalPart

        # The real world part exists in an MDB.
        self.real_world_part_mdb_md: abq_md.AbaqusMdbMetadata



        








