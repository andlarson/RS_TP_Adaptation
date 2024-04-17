This directory was created to test various options for creating a volumetric mesh
    on a .stl representation of a geometry.
This is currently a bottleneck in the implementation.

The mesh generation tools considered include:
    
    More Promising:
        Abaqus
        Meshlab / Pymeshlab 
        TetGen / Pyvista 
        fTetWild
        Quartet
        CGAL / Pygalmesh 
        Gmsh (uses a combination of meshing algorithms, including TetGen)
    
    Less Promising (closed source, old, not discussed in literature/online):
        NetGen
        Hypermesh
        OpenFOAM
        COMSOL Multiphysics
        Rhino 3D
        SlicerSegmentMesher

Results:
    
    For each tool, I used the .stl files included in the directories:
        00000001/ (Hard, ~40000 vertices + complex shape)
        00000007/ (Easy, ~250 vertices + simple shape)
        00000008/ (Medum, ~1500 vertices + complex shape)
        00000009/ (Easy, ~600 vertices + simple shape)
        00000011/ (Very Hard, ~250000 vertices + very complex shape)
        00000017/ (Easy, ~160 vertices + simple shape)
        00000022/ (Very Easy, ~10 vertices + cube)
        00000024/ (Hard, ~8000 vertices + complex shape)
        00000029/ (Medium, ~5000 vertices + simple shape)
    ASCII versions of these are included in the select_ascii_stl_from_ABC_dataset/ directory.
    
    THESE .STL FILES COME FROM CAD MODELS (via the ABC dataset), AND AS SUCH, MAY HAVE DIFFERENT
        CHARACTERISTICS THAN THOSE WE WILL ACTUALLY GET FROM PART SCANS!
        IN PARTICULAR, MANY OF THESE MODELS SEEM TO HAVE SELF-INTERSECTING TRIANGLES.
    
    Abaqus:
        Results in abaqus_meshing_results/ directory:
            00000001/ -> Failure. "180 not in list" failure during mesh to geometry conversion.    
            00000007/ -> Success using default free meshing technique.
            00000008/ -> Failure. "6086 not in list" failure during mesh to geometry conversion.
            00000009/ -> Success using default free meshing technique.
            00000011/ -> Failure. "98834 not in list" failure during mesh to geometry conversion. 
            00000017/ -> Failure. Invalid geometry resulting from mesh to geometry conversion.
            00000022/ -> Success using default free meshing technique. 
            00000024/ -> Failure. Invalid geometry resulting from mesh to geometry conversion. Happens when instancing into Assembly is attempted.
            00000029/ -> Failure. Invalid geometry resulting from mesh to geometry conversion.
         
    Meshlab / Pymeshlab:
        Volumetric meshing not supported.

    TetGen:
        For these tests, I give each a maximum of 5 minutes. If the maximum time
            elapses, then I write down the step at which the test got stuck.
        To give TetGen the best chance of success within the time window, I disable 
            mesh improvement (setting quality=False). I also turn on verbose output
            (setting verbose=2) to get an idea of which operations take a lot of
            time. All other settings are left as the defaults.
        To speed up the tests, I run the tests in parallel because the algorithm 
            is single threaded and CPU-bound.
        All time results are self-reported by TetGen.

        Script used to produce results in tetgen_meshing_results/ directory:
            00000001/ -> Ran out of time. Stuck on suppressing Steiner points.
            00000007/ -> Success in .006s. Success w/mesh improvement turned on in .04s.
            00000008/ -> Success in .0323s. Success w/mesh improvement turned on in .12s. 
            00000009/ -> Success in .011s. Success w/mesh improvement turned on in .3s. 
            00000011/ -> Failure due to self intersections. 
            00000017/ -> Success in .01s. Success w/mesh improvement turned on in .04s.
            00000022/ -> Success in .0004s. Success w/mesh improvement turned on in .0007s.
            00000024/ -> Failure due to self intersections.
            00000029/ -> Success in .08s. Success w/mesh improvement turned on in 2.5s.

    fTetWild:
        For these tests, I give each a maximum of 5 minutes. If the maximum time
            elapses, then I write down the step at which the test got stuck. 
        To give fTetWild the best chance of success within the time window, I let
            it use up to 8 threads. For each test, I use a stop energy of 1000 
            (less mesh optimization) and a stop energy of 10 (more mesh optimization). 
            I leave the ideal edge length and deviation epsilon as the default 
            values. I sometime require that the output mesh is manifold. It is 
            unclear if FEA solvers care about input meshes being manifold.
        Time results are wall-clock time recorded using the 'time' program.

        Results in ftetwild_meshing_results/ directory:
            00000001/ -> Stop Energy 1000 + Not Req. Manifold: Success in 1m27s.  
                         Stop Energy 10 + Not Req. Manifold:   Success in 1m35s.
                         Stop Energy 10 + Req. Manifold:       Success in 1m57s.
            00000007/ -> Stop Energy 1000 + Not Req. Manifold: Success in 1s.
                         Stop Energy 10 + Not Req. Manifold:   Success in 2s.
                         Stop Energy 10 + Req. Manifold:       Success in 2s.
            00000008/ -> Stop Energy 1000 + Not Req. Manifold: Success in 3s.
                         Stop Energy 10 + Not Req. Manifold:   Success in 4s.
                         Stop Energy 10 + Req. Manifold:       Success in 4s.
            00000009/ -> Stop Energy 1000 + Not Req. Manifold: Success in 18s.
                         Stop Energy 10 + Not Req. Manifold:   Success in 22s.
                         Stop Energy 10 + Req. Manifold:       Success in 24s.
            00000011/ -> Stop Energy 1000 + Not Req. Manifold: Success in 18m37s.
                         Did not perform additional tests. 
            00000017/ -> Skipped testing.
            00000022/ -> Skipped testing.
            00000024/ -> Stop Energy 1000 + Not Req. Manifold: Success in 11s.
                         Stop Energy 10 + Not Req. Manifold:   Success in 14s.
                         Stop Energy 10 + Req. Manifold:       Success in 14s.
            00000029/ -> Skipped testing.

        It turns out that allowing an epsilon-envelope makes our life hard. The
            trench may have moved during the meshing process! What happens to
            runtime if we set epsilon=0? The following tests include this
            modification.
            00000001/ -> Stop Energy 10 + Req. Manifold + Epsilon=0: Cancelled test after 10 minutes.
            00000007/ -> Stop Energy 10 + Req. Manifold + Epsilon=0: Cancelled test after 5 minutes.
