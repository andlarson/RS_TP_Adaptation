"""
This file provides utilities which are generally useful when evaluating a full
    machining process.
"""

import pathlib
import shutil
import os

def nuke_and_remake(path: str) -> None:
    """Deletes the directory tree rooted at path and then recreates is."""

    fs_path = pathlib.Path(path)
    if fs_path.exists():
        if fs_path.is_dir():
            shutil.rmtree(fs_path)
        else:
            raise RuntimeError("Something exists but it's not a directory!") 
    os.mkdir(fs_path)
