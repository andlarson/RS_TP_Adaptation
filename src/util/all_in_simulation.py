"""
This file provides utilities which are useful when evaluating a full machining 
    process in simulation. In particular, the utilities make it easier to run 
    a "real life" simulation process in parallel to the potential tool path 
    simultions, the stress recovery routine, etc.

Note that many of these utilities use the Abaqus Python API directly to modify
    state, but don't record any metadata. This contrasts with the utilities 
    provided in the abaqus_shim.py file, which use the Abaqus Python API and
    therefore record and update metadata as necessary. These utilities are only 
    conveniences which make it easier to do a "real life" simulation procedure, 
    and therefore are not called within the library itself - so no metadata
    recording is necessary.
"""

from typing import Any

import src.core.abaqus.abaqus_shim as shim



def rename_model(mdb: Any) -> None:
    """Renames the single model in an MDB so that it has a standard model name.
       
       Args:
           mdb: Abaqus MDB object. The MDB to modify. Assumed to have only a
                    single model.
    
       Returns:
           None.
    
       Raises:
           None.
    """
    
    if len(mdb.models) != 1:
        raise RuntimeError("There isn't exactly one model in the MDB!")
    
    cur_name = mdb.models.keys()[0]
    mdb.models.changeKey(fromName=cur_name, toName=shim.STANDARD_MODEL_NAME)
