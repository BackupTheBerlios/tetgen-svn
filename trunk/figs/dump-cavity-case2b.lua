-- Generated by TetView. Fri Jan 16 16:24:29 2009


p=gridgen.tvCreate()
p:create_pt(1005, -31, 13, 10.670999224659344)
p:create_pt(434, -30.899989999999999, 13.041425498372934, 15)
p:create_pt(1008, -30.899989999999999, 12.899990000000001, 10.67415186829335)
p:create_pt(1010, -32.662649999999999, 15.546670000000001, 11.25)
p:create_pt(283, -30.899989999999999, 16.5, 15)
p:create_pt(175, -32.662649999999999, 15.546670000000001, 15)
p:create_pt(1016, -32.965060000000001, 15.84479, 11.25)
p:create_pt(452, -30.899989999999999, 12.899990000000001, 7.5)
p:create_pt(628, -31, 13, 7.5)
p:create_pt(456, -30.899989999999999, 16.5, 7.5)
p:create_pt(723, -31, 13, 4.2737983604712007)
p:create_pt(666, -30.899989999999999, 12.899990000000001, 4.2714574123788616)
p:create_pt(635, -32.662649999999999, 15.546670000000001, 3.75)
p:create_pt(676, -32.965060000000001, 15.84479, 3.75)
p:create_pt(21, -30.899989999999999, 16.5, 0)
p:create_pt(648, -30.899989999999999, 13.041425498372934, 0)
p:create_pt(129, -32.662649999999999, 15.546670000000001, 0)

p:draw_subface(1008, 434, 1005)
p:draw_subface(434, 1010, 1005)
p:draw_subface(175, 283, 1010)
p:draw_subface(434, 175, 1010)
p:draw_subface(283, 1016, 1010)
p:draw_subface(452, 1008, 1005)
p:draw_subface(628, 452, 1005)
p:draw_subface(456, 628, 1005)
p:draw_subface(1010, 456, 1005)
p:draw_subface(1016, 456, 1010)
p:draw_subface(283, 456, 1016)
p:draw_subface(723, 452, 628)
p:draw_subface(456, 723, 628)
p:draw_subface(666, 452, 723)
p:draw_subface(456, 635, 723)
p:draw_subface(456, 676, 635)
p:draw_subface(676, 21, 635)
p:draw_subface(648, 666, 723)
p:draw_subface(635, 648, 723)
p:draw_subface(456, 21, 676)
p:draw_subface(129, 648, 635)
p:draw_subface(21, 129, 635)

p:set("draw", "geom")
rnd=glvCreate(0, 0, 500, 500, "TetView")
p:plot(rnd)
glvWait()
