-- Generated by TetView. Tue Jan 20 09:25:42 2009


p=gridgen.tvCreate()
p:create_pt(775, 0.026363000000000001, -0.053999999999999999, -0.0085660000000000007)
p:create_pt(776, 0.0085660000000000007, -0.053999999999999999, -0.026363000000000001)
p:create_pt(772, -0.0085660000000000007, -0.053999999999999999, -0.026363000000000001)
p:create_pt(741, 0, -0.058000000000000003, -0.0264)
p:create_pt(759, -0.0085660000000000007, -0.053999999999999999, 0.026363000000000001)
p:create_pt(771, 0.026363000000000001, -0.053999999999999999, 0.0085660000000000007)
p:create_pt(766, -0.026363000000000001, -0.053999999999999999, -0.0085660000000000007)

p:draw_subface(772, 776, 775)
p:draw_tet(772, 776, 759, 741)
p:draw_subface(772, 776, 771)
p:draw_subface(772, 776, 766)
p:set("draw", "geom")
rnd=glvCreate(0, 0, 500, 500, "TetView")
p:plot(rnd)
glvWait()