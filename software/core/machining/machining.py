import core.part.part as part
import core.abaqus.catch_all as catch_all
import storage.storage as storage

import os



class MachiningProcess:
   
    # A name is only necessary if the part is already in an Abaqus MDB.
    def __init__(self, name, part, tool_passes, in_sim):
    # type: (str, Part, ToolPasses, bool) -> None

        self.part = part
        self.tool_passes = tool_passes
        self.in_sim = in_sim

        if part.initial_part_rep.format == part.ABAQUS_MDB:
            if name != None:
                print("Name for MachiningProcess not used. The part was created\
                        as a .cae file, so it already has a name!")    
            self.mdb = catch_all.use_mdb(part.initial_part_rep.rep)
        else:
            # Save to default location with custom name. 
            self.mdb = catch_all.create_mdb(storage.MDB_save_area() + name)

    
    def sim_next_tool_pass(self):
    # type: (None) -> None
        
        # TODO: Eventually remove the roadmap below.
        # If the part has an Abaqus representation, spin up the model in Abaqus.
        #   Otherwise build the part in Abaqus.
        # 
        # Then build the next tool pass as a part in Abaqus. This must always be
        #   built.
        #
        # Then constuct the Assembly via instances of the parts. Combine the
        #   instances to simulate the post-cut geometry.
        #
        # Then simulate the stress redistribution and the associated deformation.
        #
        # Save off the resulting deformed geometry and redistributed stress
        #   profile.
   
        # Build the next tool pass path as a part.
        tool_pass_geom = self.tool_passes.pop()
        tool_pass_name = "cut_" + str(self.tool_passes.cuts_done)
        catch_all.build_part(tool_pass_name, tool_pass_geom, self.mdb)

        #  
       

    # TODO:
    # Use multiple consecutive simulations to simulate the results of multiple
    #   consecutive tool passes. We could also implement a naive-approach which
    #   could simulate multiple tool passes with a single simulation with the
    #   implicit assumption that each tool path, except possibly the last one,
    #   would not cause a large deformation.
    def sim_consecutive_tool_passes(self, cnt):
    # type: (int) -> None

        pass 
       
         

    # This does the backward direction (estimate residual stress which caused
    #   a particular deformation).
    # If real-world data is available, that data can be used.
    def estimate_stress_via_last_tool_pass(self, use_real_meas_data):
    # type: (bool) -> StressProfile

        if use_real_meas_data and self.in_sim:
            raise RuntimeError("There is no data from the real world available \ 
                                for this machining process, so the stress \
                                estimation cannot use any real world data.")


    







