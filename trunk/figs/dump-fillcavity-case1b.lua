-- Generated by TetView. Fri Jan 16 10:33:17 2009


p=gridgen.tvCreate()
p:create_pt(1020, -29.359660000000002, 11.90873, 11.25)
p:create_pt(1035, -27.410820000000001, 11.44139, 11.25)
p:create_pt(591, -29.501888444572362, 11.122324071516163, 15)
p:create_pt(1024, -29.729500000000002, 9.8638200000000005, 11.25)
p:create_pt(453, -29.729500000000002, 9.8638200000000005, 7.5)
p:create_pt(516, -29.359660000000002, 11.90873, 7.5)
p:create_pt(772, -27.410820000000001, 11.44139, 7.5)
p:create_pt(919, -27.410820000000001, 11.44139, 3.75)
p:create_pt(662, -29.359660000000002, 11.90873, 3.75)
p:create_pt(655, -29.729500000000002, 9.8638200000000005, 3.75)
p:create_pt(632, -29.501888444572362, 11.122324071516163, 0)

p:draw_tet(1024, 591, 1035, 1020)
p:draw_tet(453, 1024, 1035, 1020)
p:draw_tet(1020, 516, 1035, 453)
p:draw_tet(1035, 516, 772, 453)
p:draw_tet(516, 919, 772, 453)
p:draw_tet(516, 662, 919, 453)
p:draw_tet(919, 662, 632, 655)
p:draw_tet(655, 453, 919, 662)
p:draw_tet(1020, 1024, 591, 453)
p:draw_tet(516, 453, 1020, 591)
p:draw_tet(662, 453, 516, 632)
p:draw_tet(655, 662, 632, 453)
p:draw_tet(516, 453, 591, 632)
p:set("draw", "geom")
rnd=glvCreate(0, 0, 500, 500, "TetView")
p:plot(rnd)
glvWait()
