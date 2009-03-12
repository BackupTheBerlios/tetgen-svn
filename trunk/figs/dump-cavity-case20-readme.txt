Opening medium/yuanli.mesh.
Delaunizing vertices.
  Incremental inserting vertices.
Delaunay seconds:  0.5
Creating surface mesh.
  Marking acute vertices.
  8436 acute vertices.
Surface meshing seconds:  0.1
Recovering boundaries.
  Delaunizing segments.
  4383 protecting points.
  Constraining facets.
  37 subedge flips.
  147 cavities remeshed.
Boundary recovery seconds:  0.36
Removing exterior tetrahedra.
  Removing flat boundary tetrahedra.
constrain.cxx:3267: failed assertion `(point) casface.tet[7] != dummypoint'
./runtest.sh: line 78:  1343 Abort trap              $1 $2 $f

  Recover facet #382: 24 subfaces, 25 vertices.
    Scout subface (13205, 12642, 15198) (30472).
    Scout subface (12642, 13205, 14664) (30472).
    Scout subface (13205, 14663, 14664) (30472).
    Scout subface (13044, 14664, 14663) (30472).
    Found a crossing tet (14664, 13044, 11491, 15211).
    Formed cavity: 11 (3) cross tets (edges).

pedge(11565,14997)
-- !!!! Should be a subsegment!!!!

(gdb) p tetarray->objects

i = 55.
    flip 2-to-2: (11565, 14997, 13615, 13750)

(gdb) p pteti(11565, 14997,13615,13750)
Tetra x2016470 with loc(0) ver(0):  ori = -1.79513e-07.
      
      [0] x175da70  loc(0) ver(4)
      [1] x2028c20  loc(1) ver(0)
      [2] x1741fc0  loc(1) ver(2)
      [3] x2042ab0  loc(1) ver(2)
      Org [0] x16dee08 (5.37405810533,0.212478607393,2.37345369556) 14997
      Dest[1] x16ceae8 (5.38280518588,0.215150463317,2.37036740283) 13615
      Apex[2] x5819c8 (5.38303453135,0.214666785205,2.37236965835) 11565
      Oppo[3] x16d0438 (5.37936,0.20362,2.37394) 13750
      [0] x17bee78 0. (*)
      [1] NULL.
      [2] NULL.
      [3] x11faf38 0.
$5 = void

tetgen:~/sync/trunk/src si$ svn ci -m "debug in constrainedfacets(), dump-cavity-case19"
sihang@svn.berlios.de's password: 
Sending        src/constrain.cxx
Sending        src/meshstat.cxx
Sending        src/tetgen.h
Transmitting file data ...
Committed revision 411.
tetgen:~/sync/trunk/src si$ 