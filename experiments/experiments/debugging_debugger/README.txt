This directory was created to debug the Abaqus PDE debugger.

In particular, I believe that the debugger is caching files, and when the files
    change on the underlying file system, those changes are not captured by the
    debugger. This forces me to constantly restart the debugger.
