################################

  -- The error message

################################

tetgen:~/benchmark/medium si$ ../tg -pcVNEFC crystal03_cut.d_cut.smesh 
Opening crystal03_cut.d_cut.smesh.
Opening crystal03_cut.d_cut.node.
Delaunizing vertices.
  Incremental inserting vertices.
Delaunay seconds:  0.68
Creating surface mesh.
  Marking acute vertices.
  9740 acute vertices.
Recovering boundaries.
  Delaunizing segments.
  11838 protecting points.
  Constraining facets.
  p:draw_tet(22481, 17720, 17366, 23219) -- top tet.
  p:draw_tet(17720, 22481, 18908, 17367) -- bot tet.
  Found a non-Delaunay edge (18908, 17720)
    Cut tet (17725, 17366, 17720, 17714)
    Cut tet (17714, 17366, 17720, 20163)
    Cut tet (17714, 17366, 20163, 17367)
    Lawson flip 18 faces.
    flip 2-to-3: (17725, 17720, 17366, 22740, 17714)
      Queue a flipped subface (17725, 17366, 17720).
    flip 2-to-3: (22740, 17366, 17720, 17714, 23545)
    flip 3-to-2: (23545, 17714, 20163, 17720, 17366)
      Queue a flipped segment (17720, 17366).
    flip 3-to-2: (17714, 23545, 17367, 17366, 20163)
    51 faces stacked, 4 flips.
Warning:  Point #23545 is duplicated with Point #23545. Ignored!
  Delaunizing segments.
  11838 protecting points.
  207 subedge flips.
  1400 cavities remeshed.
Segment and facet seconds:  1.23

NOT writing a .node file.
NOT writing an .ele file.
NOT writing an .face file.

Output seconds:  0.14
Total running seconds:  2.05

  Checking consistency of the mesh boundary...
  Mesh boundaries connected correctly.
  Checking constrained Delaunay property of the mesh...
  !! Non-locally Delaunay (17367, 23545, 19559) - 18908, 20163
  !! Non-locally Delaunay (19559, 23545, 18908) - 17367, 22481
  !! Non-locally Delaunay (19559, 23545, 20163) - 22830, 17367
  !! Non-locally Delaunay (20163, 23545, 22830) - 19559, 22481
  !! Non-locally Delaunay (22830, 23545, 19559) - 20163, 22481
  !! Non-locally Delaunay (17714, 20163, 17367) - 23545, 11757
  !! Non-locally Delaunay (22481, 23545, 22830) - 20163, 19559
  !! Non-locally Delaunay (17367, 17714, 17366) - 11757, 23545
  !! !! !! !! Found 8 non-Delaunay faces.
main.cxx:263: failed assertion `0'

############################

  -- To produce this case

############################

mac-002:~/sync/trunk/src si$ svn status -u
sihang@svn.berlios.de's password: 
sihang@svn.berlios.de's password: 
?                   tetgen
?                   libtet.a
M             364   delaunay.cxx
Status against revision:    364

