-- Generated by TetView. Fri Jan 23 01:13:02 2009


p=gridgen.tvCreate()
p:create_pt(29, -1, 1.3090169943749475, 3.4270509831248424)
p:create_pt(31, -1, -1.3090169943749475, 3.4270509831248424)
p:create_pt(119, -2.6180339887498949, -1.3090169943749475, 2.4270509831248424)
p:create_pt(71, -3.1180339887498949, -0.5, 2.1180339887498949)
p:create_pt(79, -1.8090169943749475, -1.6180339887498949, 2.9270509831248424)
p:create_pt(5, -0.5, 0.5, 3.7360679774997898)
p:create_pt(7, -0.5, -0.5, 3.7360679774997898)

p:draw_subface(119, 31, 29)
p:draw_tet(119, 31, 79, 71)
p:draw_tet(119, 31, 5, 71)
p:draw_tet(119, 31, 7, 79)
p:set("draw", "geom")
rnd=glvCreate(0, 0, 500, 500, "TetView")
p:plot(rnd)
glvWait()