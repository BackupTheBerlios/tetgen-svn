tetgen:~/sync/trunk/ndelaunay si$ ../src/tetgen -pcV -C 021_145_299_pr2_cut_cut.mesh 
Opening 021_145_299_pr2_cut_cut.mesh.
Delaunizing vertices.
  Incremental inserting vertices.
Delaunay seconds:  0.05
Creating surface mesh.
  Marking acute vertices.
  594 acute vertices.
Surface meshing seconds:  0.01
Recovering boundaries.
  Delaunizing segments.
  2046 protecting points.
  Constraining facets.
  p:draw_tet(2946, 1625, 2918, 456) -- top tet.
  p:draw_tet(1625, 2946, 2935, 1437) -- bot tet.
  Found a non-Delaunay edge (2935, 1625)
    Lawson flip 15 faces.
    15 faces stacked, 0 flips.
  Delaunizing segments.
    Lawson flip 36 faces.
    36 faces stacked, 0 flips.
  2047 protecting points.
  p:draw_tet(1625, 2918, 2946, 456) -- top tet.
  p:draw_tet(2918, 1625, 2935, 1437) -- bot tet.
  Found a non-Delaunay edge (1625, 2935)
    Lawson flip 63 faces.
    63 faces stacked, 0 flips.
  Delaunizing segments.
  2047 protecting points.
  p:draw_tet(1625, 2918, 2946, 456) -- top tet.
  p:draw_tet(2918, 1625, 2935, 1437) -- bot tet.
  Found a non-Delaunay edge (1625, 2935)
    Lawson flip 15 faces.
    15 faces stacked, 0 flips.
  Delaunizing segments.
  2047 protecting points.
  Face (2229, 982, 778) - 4 is missing
  p:draw_subseg(2266, 958)
  dump 22 subfaces to facet2.lua
  dump 16 topfaces to cavity.lua
  dump 46 subfaces to facet1.lua
Writing 021_145_299_pr2_cut_cut.1.node.
constrain.cxx:2649: failed assertion `0'

tetgen:~/sync/trunk/src si$ svn ci -m "updated formedgecavity()"
sihang@svn.berlios.de's password: 
Sending        src/constrain.cxx
Sending        src/meshstat.cxx
Sending        src/tetgen.h
Transmitting file data ...
Committed revision 397.
