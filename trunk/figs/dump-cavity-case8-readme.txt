*******************************************
-- Error Message
*******************************************

tetgen:~/sync/trunk/data si$ ../src/tetgen -pV -C 386DXCPU.d_cut.mesh 
Opening 386DXCPU.d_cut.mesh.
Delaunizing vertices.
  Incremental inserting vertices.
Delaunay seconds:  1.8
Creating surface mesh.
Recovering boundaries.
  Marking acute vertices.
  22791 acute vertices.
  Delaunizing segments.
  65431 protecting points.
  Constraining facets.
  p:draw_tet(55636, 49179, 44314, 77435) -- top tet.
  p:draw_tet(49179, 55636, 55637, 33642) -- bot tet.
  Found a non-Delaunay edge (55637, 49179)
  Delaunizing segments.
  65431 protecting points.
  p:draw_tet(64420, 64418, 53790, 64417) -- top tet.
  p:draw_tet(64418, 64420, 60443, 36896) -- bot tet.
  Found a non-Delaunay edge (64418, 53790)
  Delaunizing segments.
  65433 protecting points.
  p:draw_tet(59732, 31080, 90681, 10522) -- top tet.
  p:draw_tet(31080, 59732, 82127, 10521) -- bot tet.
  Found a non-Delaunay edge (59732, 82127)
  Delaunizing segments.
  65433 protecting points.
  1302 subedge flips.
  4564 cavities remeshed.
Segment and facet seconds:  5.9

Writing 386DXCPU.d_cut.1.node.
Writing 386DXCPU.d_cut.1.ele.
Writing 386DXCPU.d_cut.1.face.

Output seconds:  2.46
Total running seconds:  10.17

  Checking consistency of the mesh boundary...
  Mesh boundaries connected correctly.
  Checking constrained Delaunay property of the mesh...
  !! Non-locally Delaunay (32735, 95086, 35274) - 32633, 60316
  !! Non-locally Delaunay (35274, 95086, 32633) - 32735, 60316
  !! Non-locally Delaunay (35274, 95086, 60316) - 32633, 32735
  !! !! !! !! Found 3 non-Delaunay faces.
main.cxx:254: failed assertion `0'

*******************************************
-- To reproduce this example
*******************************************

tetgen:~/sync/trunk/src si$ svn status -u
sihang@svn.berlios.de's password: 
sihang@svn.berlios.de's password: 
?                   tetgen
?                   libtet.a
Status against revision:    351
tetgen:~/sync/trunk/src si$ 

