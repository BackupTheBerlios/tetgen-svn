tetgen:~/sync/trunk/ndelaunay si$ ../src/tetgen -pV -C Y565_conveyor_cut_cut.d.mesh Opening Y565_conveyor_cut_cut.d.mesh.
Delaunizing vertices.
  Incremental inserting vertices.
Delaunay seconds:  3.37
Creating surface mesh.
  Marking acute vertices.
  9893 acute vertices.
Surface meshing seconds:  0.74
Recovering boundaries.
  Delaunizing segments.
  5203 protecting points.
  Constraining facets.
  p:draw_subface(5872, 4124, 3793) -- 8 is missing
  p:draw_subseg(4125, 3965)
constrain.cxx:1623: failed assertion `dir != ACROSSVERT'


(gdb) p outnodes(0)
Writing ../ndelaunay/Y565_conveyor_cut_cut.d.1.node.

(gdb) p print_facearray(cavfaces)
p:draw_subface(6591, 420, 3952) -- 1
p:draw_subface(420, 3793, 3952) -- 2
p:draw_subface(3793, 19242, 3952) -- 3
p:draw_subface(19242, 3953, 3952) -- 4
p:draw_subface(3953, 6591, 3952) -- 5
p:draw_subface(19242, 6591, 3953) -- 6
p:draw_subface(420, 56669, 3793) -- 7
p:draw_subface(56669, 5872, 3793) -- 8
p:draw_subface(5872, 4124, 3793) -- 9
p:draw_subface(4124, 4125, 3793) -- 10
p:draw_subface(4125, 3792, 3793) -- 11
p:draw_subface(3792, 19242, 3793) -- 12
p:draw_subface(5872, 56214, 4124) -- 13
p:draw_subface(51078, 3792, 4125) -- 14
p:draw_subface(51078, 19242, 3792) -- 15
p:draw_subface(56214, 421, 4124) -- 16
p:draw_subface(421, 3965, 4124) -- 17
p:draw_subface(3965, 19238, 4124) -- 18
p:draw_subface(19238, 4125, 4124) -- 19
p:draw_subface(421, 19238, 3965) -- 20
p:draw_subface(19238, 51078, 4125) -- 21

(gdb) p psubface(5872,4125, 3965)
  !! Not exist.
$9 = void
(gdb) p psubface(4124,4125, 3793)
  sub x3371118 (4125, 4124, 3793) mark=34876
    area=0.000807931, lengths: 0.0402757, 0.04012, 0.0568485
  seg x33c58c8 (4124, 4125)
  sub x3371ef8 (4124, 3793, 3794)
  sub x3371ec8 (3793, 4125, 3792)
  tet x20a6d40 (4124, 4125, 3793, 3965) 1
  tet x2615ea0 (4125, 4124, 3793, 57947) 3
  Non-Delaunay (sign = -1e-08).

p:draw_subseg(5872,3965)  -- A crossing edge

mac-002:~/sync/trunk/src si$ svn status -u
sihang@svn.berlios.de's password: 
sihang@svn.berlios.de's password: 
?                   tetgen
?                   libtet.a
Status against revision:    398
