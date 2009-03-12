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


Recover facet #515: 13 subfaces, 15 vertices.
  ....
meshstat.cxx:1009: failed assertion `0'
  tet x1725e50 (14388, 14736, 14623, 14393)Abort trap

(gdb) p pedge(14736,14388)
$2 = true
  tet x3599a0 (14388, 14736, 14389, 14392) (seg)
  tet x7b36c0 (14388, 14736, 14392, 14391) (seg)
  tet x2017460 (14388, 14736, 14391, 11996) (seg) (sub)
  tet x2017430 (14388, 14736, 11996, 14389) (seg)
(gdb) p pedge(14736,14388)
$4 = true
  tet x3599a0 (14388, 14736, 14389, 14392)
  tet x7b36c0 (14388, 14736, 14392, 14391)
  tet x2017460 (14388, 14736, 14391, 11996)
  tet x2017430 (14388, 14736, 11996, 14389)

========================================================================
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

Recover facet #382: 24 subfaces, 25 vertices.
    Scout subface (13205, 12642, 15198) (30472).
    Scout subface (12642, 13205, 14664) (30472).
    Scout subface (13205, 14663, 14664) (30472).
    Scout subface (13044, 14664, 14663) (30472).
    Found a crossing tet (14664, 13044, 11491, 15211).
    Formed cavity: 11 (3) cross tets (edges).
    Delaunizing cavity: 10 points, 11 faces.
    Create init tet (14664, 13044, 15211, 13810)



$1 = 82

   i = 86.
    flip 2-to-2: (13750, 14997, 14153, 13820)
constrain.cxx:3306: failed assertion `(point) casface.tet[7] != dummypoint'

(gdb) p pteti(13750, 14997, 13820, 14153)
Tetra x1621ea0 with loc(0) ver(0):  ori = -1.83122e-07.
      
      [0] x203bd90  loc(2) ver(4)
      [1] x2028c80  loc(2) ver(2)
      [2] x160cae0  loc(3) ver(2)
      [3] x2028c20  loc(3) ver(2)
      Org [0] x16d4fc8 (5.3689119002,0.212212475443,2.37035692585) 14153
      Dest[1] x16d1158 (5.36905897063,0.211259971421,2.37405741557) 13820
      Apex[2] x16d0438 (5.37936,0.20362,2.37394) 13750
      Oppo[3] x16dee08 (5.37405810533,0.212478607393,2.37345369556) 14997
      [0] NULL. (*)
      [1] NULL.
      [2] NULL.
      [3] NULL.
      [4] x1257b78 0.
      [5] NULL.
      [0] NULL. (*)
      [1] x1449668 0.
      [2] x16afbb8 0.
      [3] NULL.

(gdb) set searchtet.tet=0x1621ea0
(gdb) set searchtet.loc=1
(gdb) set searchtet.ver=0
(gdb) p ptet(&searchtet)
Tetra x1621ea0 with loc(1) ver(0):  ori = -1.83122e-07.
      
      [0] x203bd90  loc(2) ver(4)
      [1] x2028c80  loc(2) ver(2)
      [2] x160cae0  loc(3) ver(2)
      [3] x2028c20  loc(3) ver(2)
      Org [0] x16d4fc8 (5.3689119002,0.212212475443,2.37035692585) 14153
      Dest[3] x16dee08 (5.37405810533,0.212478607393,2.37345369556) 14997
      Apex[1] x16d1158 (5.36905897063,0.211259971421,2.37405741557) 13820
      Oppo[2] x16d0438 (5.37936,0.20362,2.37394) 13750
      [0] NULL.
      [1] NULL.
      [2] NULL.
      [3] NULL. (*)
      [4] x1257b78 0.
      [5] NULL.
      [0] NULL.
      [1] x1449668 0. (*)
      [2] x16afbb8 0.
      [3] NULL.
$7 = void
(gdb) set searchtet.loc=2
(gdb) p ptet(&searchtet)
Tetra x1621ea0 with loc(2) ver(0):  ori = -1.83122e-07.
      
      [0] x203bd90  loc(2) ver(4)
      [1] x2028c80  loc(2) ver(2)
      [2] x160cae0  loc(3) ver(2)
      [3] x2028c20  loc(3) ver(2)
      Org [1] x16d1158 (5.36905897063,0.211259971421,2.37405741557) 13820
      Dest[3] x16dee08 (5.37405810533,0.212478607393,2.37345369556) 14997
      Apex[2] x16d0438 (5.37936,0.20362,2.37394) 13750
      Oppo[0] x16d4fc8 (5.3689119002,0.212212475443,2.37035692585) 14153
      [0] NULL.
      [1] NULL.
      [2] NULL.
      [3] NULL.
      [4] x1257b78 0. (*)
      [5] NULL.
      [0] NULL.
      [1] x1449668 0.
      [2] x16afbb8 0. (*)
      [3] NULL.
$8 = void


tetgen:~/sync/trunk/src si$ svn ci -m "debug in constrainedfacets(), dump-cavity-case19"
sihang@svn.berlios.de's password: 
Sending        src/constrain.cxx
Sending        src/meshstat.cxx
Sending        src/tetgen.h
Transmitting file data ...
Committed revision 411.
tetgen:~/sync/trunk/src si$ 

