-- Generated by TetView. Tue Feb  3 00:23:55 2009


p=gridgen.tvCreate()
p:create_pt(325, 14, 9.7474699999999999, 6.5)
p:create_pt(772, 13, 11.25, 6.5)
p:create_pt(308, 13, 7, 6.5)
p:create_pt(497, 13, 15.5, 11.716279380927025)
p:create_pt(545, 11, 10.25, 15)
p:create_pt(747, 11, 15.5, 11.25)
p:create_pt(822, 13, 9.9237974281574299, 8.5)
p:create_pt(310, 13, 15.5, 8.5)
p:create_pt(311, 13, 7, 8.5)
p:create_pt(626, 11, 15.5, 7.5)
p:create_pt(911, 11, 5, 11.25)
p:create_pt(742, 11, 5, 7.5)
p:create_pt(743, 11, 5, 3.75)
p:create_pt(652, 11, 10.25, 0)
p:create_pt(785, 13, 15.5, 3.2163119219504117)
p:create_pt(641, 11, 15.5, 3.75)

p:draw_subface(308, 772, 325)
p:draw_subface(747, 545, 497)
p:draw_subface(545, 822, 497)
p:draw_subface(822, 747, 497)
p:draw_subface(310, 747, 822)
p:draw_subface(545, 311, 822)
p:draw_subface(311, 308, 822)
p:draw_subface(308, 325, 822)
p:draw_subface(325, 772, 822)
p:draw_subface(772, 310, 822)
p:draw_subface(626, 747, 310)
p:draw_subface(545, 911, 311)
p:draw_subface(772, 626, 310)
p:draw_subface(911, 742, 311)
p:draw_subface(742, 308, 311)
p:draw_subface(742, 743, 308)
p:draw_subface(743, 652, 308)
p:draw_subface(652, 325, 308)
p:draw_subface(325, 652, 772)
p:draw_subface(652, 785, 772)
p:draw_subface(785, 641, 772)
p:draw_subface(641, 626, 772)
p:draw_subface(652, 641, 785)
p:set("draw", "geom")
rnd=glvCreate(0, 0, 500, 500, "TetView")
p:plot(rnd)
glvWait()