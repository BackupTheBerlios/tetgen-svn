tetgen:~/sync/trunk/ndelaunay si$ ../src/tetgen -pVc -C cassette_bd_pr3_cut.mesh 
Opening cassette_bd_pr3_cut.mesh.
Delaunizing vertices.
  Incremental inserting vertices.
Delaunay seconds:  0.01
Creating surface mesh.
  Marking acute vertices.
  83 acute vertices.
Surface meshing seconds:  0
Recovering boundaries.
  Delaunizing segments.
  445 protecting points.
  Constraining facets.
  Face (164, 132, 139) - 7 is missing
constrain.cxx:1720: failed assertion `0'

tetgen:~/sync/trunk/src si$ svn status -u
sihang@svn.berlios.de's password: 
sihang@svn.berlios.de's password: 
?                   tetgen
?                   libtet.a
Status against revision:    392
tetgen:~/sync/trunk/src si$ 

