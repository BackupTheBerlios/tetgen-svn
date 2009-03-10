tetgen:~/sync/trunk/ndelaunay si$ ../src/tetgen -pV -C cassette_bd_pr3_cut.mesh Opening cassette_bd_pr3_cut.mesh.
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
  p:draw_subface(164, 132, 139) -- 7 is missing
  p:draw_subseg(413, 484)
  p:draw_subface(144, 557, 353) -- 50 is missing
  p:draw_subseg(552, 161)
  p:draw_subface(557, 436, 353) -- 51 is missing
  p:draw_subseg(552, 161)
  p:draw_subface(436, 389, 353) -- 52 is missing
  p:draw_subseg(302, 161)
  p:draw_subface(240, 389, 436) -- 77 is missing
  p:draw_subseg(302, 209)
  p:draw_subface(240, 216, 389) -- 78 is missing
  p:draw_subseg(209, 475)
  dump 41 subfaces to facet1.lua
    Form edge cavity (209, 475).
  dump 16 subfaces to facet2.lua
  dump 8 topfaces to cavity.lua
    Form edge cavity (302, 209).
  dump 16 subfaces to facet2.lua
  dump 12 topfaces to cavity.lua
    Form edge cavity (302, 161).
  dump 16 subfaces to facet2.lua
  dump 7 topfaces to cavity.lua
    Form edge cavity (552, 161).
  dump 16 subfaces to facet2.lua
  dump 6 topfaces to cavity.lua
    Form edge cavity (552, 161).
  dump 16 subfaces to facet2.lua
  dump 3 topfaces to cavity.lua
    Form edge cavity (413, 484).
  dump 50 subfaces to facet2.lua
  dump 22 topfaces to cavity.lua
Writing cassette_bd_pr3_cut.1.node.
constrain.cxx:2687: failed assertion `0'
Abort trap
tetgen:~/sync/trunk/ndelaunay si$ 

tetgen:~/sync/trunk/src si$ svn ci -m "updated delaunizecavity()"
sihang@svn.berlios.de's password: 
Sending        src/constrain.cxx
Sending        src/tetgen.h
Transmitting file data ..
Committed revision 399.
tetgen:~/sync/trunk/src si$  
