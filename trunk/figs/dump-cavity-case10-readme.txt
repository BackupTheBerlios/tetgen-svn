*******************************************
-- Error Message
*******************************************

tetgen:~/sync/trunk/data si$ ../src/tetgen -pV -C shuttle9_cut.d_cut.mesh 
Opening shuttle9_cut.d_cut.mesh.
Delaunizing vertices.
  Incremental inserting vertices.
Delaunay seconds:  1.55
Creating surface mesh.
Recovering boundaries.
  Marking acute vertices.
  20558 acute vertices.
  Delaunizing segments.
  12178 protecting points.
  Constraining facets.
  p:draw_tet(36198, 36195, 29586, 36111) -- top tet.
  p:draw_tet(36195, 36198, 30993, 36197) -- bot tet.
  Found a non-Delaunay edge (30993, 36195)
  p:draw_tet(35298, 33799, 30993, 36966) -- top tet.
  p:draw_tet(33799, 35298, 29586, 35924) -- bot tet.
  Found a non-Delaunay edge (29586, 33799)
    Lawson flip 33 faces.
    33 faces stacked, 0 flips.
  p:draw_tet(35015, 31816, 25070, 34759) -- top tet.
  p:draw_tet(31816, 35015, 31817, 31874) -- bot tet.
  Found a non-Delaunay edge (25070, 35015)
  Delaunizing segments.
    Lawson flip 30 faces.
    flip 2-to-3: (32613, 37000, 33647, 31817, 35742)
    flip 3-to-2: (35742, 31817, 36999, 33647, 37000)
    42 faces stacked, 2 flips.
  12179 protecting points.
  p:draw_tet(35741, 29799, 36186, 35459) -- top tet.
  p:draw_tet(29799, 35741, 35114, 36782) -- bot tet.
  Found a non-Delaunay edge (35741, 35114)
  Delaunizing segments.
  12179 protecting points.
  p:draw_tet(35741, 29799, 36186, 35459) -- top tet.
  p:draw_tet(29799, 35741, 35114, 37001) -- bot tet.
  Found a non-Delaunay edge (35741, 35114)
  65 subedge flips.
  3493 cavities remeshed.
Segment and facet seconds:  1.67

Writing shuttle9_cut.d_cut.1.node.
Writing shuttle9_cut.d_cut.1.ele.
Writing shuttle9_cut.d_cut.1.face.

Output seconds:  0.9
Total running seconds:  4.12

  Checking consistency of the mesh boundary...
  Mesh boundaries connected correctly.
  Checking constrained Delaunay property of the mesh...
  !! Non-locally Delaunay (36998, 35668, 36444) - 36966, 35298
  !! !! !! !! Found 1 non-Delaunay faces.
main.cxx:254: failed assertion `0'

*******************************************
-- To reproduce this example
*******************************************

tetgen:~/sync/trunk/src si$ svn status -u
sihang@svn.berlios.de's password: 
sihang@svn.berlios.de's password: 
?                   tetgen
?                   libtet.a
Status against revision:    359


