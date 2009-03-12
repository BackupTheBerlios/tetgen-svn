tetgen:~/sync/trunk/ndelaunay si$ ../src/tetgen -pV -C Y1970_BBEARING_cut.mesh 

  p:draw_subface(58327, 58326, 58478) -- 1 is missing
  p:draw_subface(58478, 58347, 58327) -- 2 is missing
  p:draw_subface(58567, 58347, 58478) -- 4 is missing
  p:draw_subface(58411, 58347, 58567) -- 5 is missing
  p:draw_subface(58326, 58478, 58047) -- 115 is missing
  p:draw_subface(58478, 58327, 58047) -- 116 is missing
  p:draw_subface(58327, 58478, 48254) -- 119 is missing
  p:draw_subface(58047, 57925, 58478) -- 125 is missing
  p:draw_subface(58047, 58478, 48254) -- 128 is missing
  p:draw_subface(34838, 58047, 57925) -- 129 is missing
  p:draw_subface(34838, 48254, 58047) -- 133 is missing
constrain.cxx:2294: failed assertion `marktested(bottet)'
Abort trap
tetgen:~/sync/trunk/ndelaunay si$ ../../../tetgen-dbg/cdt -pV Y1970_BBEARING_cut.mesh 
Opening Y1970_BBEARING_cut.mesh.
Constructing Delaunay tetrahedralization.
  Creating initial tetrahedralization.
  Incrementally inserting points.
  224208 Flips (T23 130996, T32 86842, T22 263, T44 6107)
Delaunay seconds:  1.4
Creating surface mesh.
  Constructing mapping from indices to points.
  Constructing mapping from points to tetrahedra.
  Unifying segments.
  Constructing mapping from points to subfaces.
  Merging coplanar facets.
  Marking acute vertices.
  Constructing mapping from points to segments.
  12240 acute vertices.
Perturbing vertices.
  Removing degenerate subfaces.
  295 break points.
Delaunizing segments.
  Constructing mapping from points to tetrahedra.
  Queuing missing segments.
  47854 protect points.
  R1: 9450,  R2: 33842,  R3: 4562.
Constraining facets.
  Constructing mapping from points to tetrahedra.
  The biggest cavity: 242 faces, 123 vertices
  Enlarged 0 times
Segment and facet seconds:  11.21
  Surface mesh seconds:  0.15
  Vertex perturbation seconds:  0.79
  Segment recovery seconds:  9.33
  facet recovery seconds:  0.94
Removing unwanted tetrahedra.
  Marking concavities for elimination.
  Marking neighbors of marked tetrahedra.
Warning:  Detecting an open face (25661, 1050, 1055).
  Deleting marked tetrahedra.
Hole seconds:  0.54
Repairing mesh.
Repair seconds:  0.02
Jettisoning redundants points.
  0 duplicated vertices have been removed.
  53837 unused vertices have been removed.

Writing Y1970_BBEARING_cut.1.node.
Writing Y1970_BBEARING_cut.1.ele.
Writing Y1970_BBEARING_cut.1.face.


tetgen:~/sync/trunk/src si$ svn status -u
sihang@svn.berlios.de's password: 
sihang@svn.berlios.de's password: 
?                   tetgen
?                   libtet.a
Status against revision:    406
tetgen:~/sync/trunk/src si$ 


