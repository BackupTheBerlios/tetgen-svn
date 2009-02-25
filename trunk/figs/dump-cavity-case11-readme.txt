*******************************************
-- Error Message
*******************************************

tetgen:~/sync/trunk/data si$ ../src/tetgen -pVc -C 386DXCPU.d_cut.mesh Opening 386DXCPU.d_cut.mesh.Delaunizing vertices.  Incremental inserting vertices.Delaunay seconds:  1.82Creating surface mesh.
  Marking acute vertices.
  22791 acute vertices.
Recovering boundaries.
  Delaunizing segments.
  65431 protecting points.
  Constraining facets.
  p:draw_tet(55636, 49179, 44314, 77435) -- top tet.
  p:draw_tet(49179, 55636, 55637, 33642) -- bot tet.
  Found a non-Delaunay edge (55637, 49179)
    Cut tet (30146, 56736, 33642, 21179)
  Delaunizing segments.
  65431 protecting points.
  p:draw_tet(64420, 64418, 53790, 64417) -- top tet.
  p:draw_tet(64418, 64420, 60443, 36896) -- bot tet.
  Found a non-Delaunay edge (64418, 53790)
  Delaunizing segments.
    Cut tet (60316, 32735, 35274, 19541)
  65433 protecting points.
  p:draw_tet(59732, 31080, 90681, 10522) -- top tet.
  p:draw_tet(31080, 59732, 82127, 10521) -- bot tet.
  Found a non-Delaunay edge (59732, 82127)
  Delaunizing segments.
  65433 protecting points.
  1302 subedge flips.
  4564 cavities remeshed.
Segment and facet seconds:  5.88

Writing 386DXCPU.d_cut.1.node.
Writing 386DXCPU.d_cut.1.ele.
Writing 386DXCPU.d_cut.1.face.

Output seconds:  2.53
Total running seconds:  10.24

  Checking consistency of the mesh boundary...
  Mesh boundaries connected correctly.
  Checking constrained Delaunay property of the mesh...
  !! Non-locally Delaunay (32735, 95086, 60316) - 32633, 95087
  !! !! !! !! Found 1 non-Delaunay faces.
main.cxx:263: failed assertion `0'

*******************************************
-- To reproduce this example
*******************************************

tetgen:~/sync/trunk/src si$ svn status -u
sihang@svn.berlios.de's password: 
sihang@svn.berlios.de's password: 
?                   tetgen
?                   libtet.a
M             363   delaunay.cxx
Status against revision:    363
