-- Generated by TetView. Sun Jan 18 21:02:13 2009


p=gridgen.tvCreate()
p:create_pt(806, 0.57832706901842623, -0.10674080497315996, -0.27989047473124101)
p:create_pt(1205, 0.60128772820007725, 0.036396925207835336, -0.3081050674020222)
p:create_pt(872, 0.55719560713750582, -0.11511955502964927, -0.25204615753031356)
p:create_pt(1234, 0.57742641084367374, -0.00094281009553565109, -0.27605837835613489)
p:create_pt(1228, 0.56158162309520265, -0.025737760501206894, -0.25477820432429038)
p:create_pt(679, 0.59688537159054356, 0.0025289953752336489, -0.30010332859866706)
p:create_pt(914, 0.5085000003950535, -0.10880332270349942, -0.18348749428869032)
p:create_pt(810, 0.55282493926918019, -0.090456304989228631, -0.24344219776339554)
p:create_pt(516, 0.58012510413234364, -0.066477396551184828, -0.27799713877078769)
p:create_pt(1128, 0.65678957045037922, 0.050413620503897562, -0.32338803976126074)
p:create_pt(849, 0.54079005154745263, -0.058273758683459034, -0.2268543048088576)
p:create_pt(868, 0.50965149982914426, -0.11078181273461898, -0.18474125250787732)
p:create_pt(719, 0.52304709325384768, -0.12822682582376982, -0.13595362847780601)
p:create_pt(1137, 0.64886623678652788, 0.068845827201202095, -0.30522416443299499)
p:create_pt(870, 0.81063126877551661, -0.0035163194048375285, -0.0049883722648547679)
p:create_pt(727, 0.52983564675558437, -0.13723488139038523, -0.13982381627982374)
p:create_pt(1456, 0.69941354502342235, -0.13009156009362205, 0.043006790770027503)
p:create_pt(552, 0.56741182051601047, -0.035755123198076749, -0.26112615180860033)
p:create_pt(829, 0.54478764704328442, -0.065142370688637416, -0.23120690634864188)
p:create_pt(938, 0.8207813942652632, -0.01859903610357283, -0.013130206552514623)
p:create_pt(1381, 0.52938785531419508, -0.1368771094528809, -0.13942952171308926)
p:create_pt(1457, 0.69654290315448653, -0.12782737723771431, 0.045466077950465533)
p:create_pt(1374, 0.69309518068275777, -0.12194373274422123, 0.048732514447732411)

p:draw_tet(1234, 872, 1205, 806)
p:draw_tet(872, 806, 1234, 1228)
p:draw_tet(1205, 1234, 806, 679)
p:draw_tet(914, 872, 1234, 1228)
p:draw_tet(872, 806, 1228, 810)
p:draw_tet(1234, 1228, 806, 516)
p:draw_tet(1234, 806, 679, 516)
p:draw_tet(1205, 679, 806, 1128)
p:draw_tet(914, 872, 1228, 849)
p:draw_tet(1234, 914, 1228, 849)
p:draw_tet(868, 872, 914, 810)
p:draw_tet(719, 914, 1234, 849)
p:draw_tet(1205, 1137, 679, 1128)
p:draw_tet(1137, 719, 1234, 870)
p:draw_tet(872, 1456, 806, 727)
p:draw_tet(1228, 872, 810, 849)
p:draw_tet(1228, 810, 806, 516)
p:draw_tet(1234, 1228, 516, 552)
p:draw_tet(679, 516, 806, 1128)
p:draw_tet(679, 1234, 516, 552)
p:draw_tet(872, 849, 914, 810)
p:draw_tet(1234, 849, 1228, 719)
p:draw_tet(872, 810, 868, 727)
p:draw_tet(914, 868, 810, 829)
p:draw_tet(679, 1137, 552, 1128)
p:draw_tet(1234, 870, 719, 1228)
p:draw_tet(1234, 1137, 870, 1228)
p:draw_tet(806, 727, 1456, 938)
p:draw_tet(938, 1128, 806, 516)
p:draw_tet(1228, 849, 810, 829)
p:draw_tet(1228, 810, 516, 552)
p:draw_tet(679, 516, 1128, 552)
p:draw_tet(914, 810, 849, 829)
p:draw_tet(719, 849, 1228, 870)
p:draw_tet(868, 727, 810, 829)
p:draw_tet(914, 868, 829, 849)
p:draw_tet(1381, 727, 868, 719)
p:draw_tet(1128, 1137, 552, 516)
p:draw_tet(1456, 727, 1381, 1457)
p:draw_tet(938, 1456, 1457, 1374)
p:draw_tet(1228, 849, 829, 552)
p:draw_tet(1228, 829, 810, 552)
p:draw_tet(810, 516, 552, 829)
p:draw_tet(1457, 1381, 719, 1374)
p:draw_tet(1374, 870, 938, 1456)
p:set("draw", "geom")
rnd=glvCreate(0, 0, 500, 500, "TetView")
p:plot(rnd)
glvWait()