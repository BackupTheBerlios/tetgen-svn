-- Generated by TetView. Thu Jan 22 10:03:45 2009


p=gridgen.tvCreate()
p:create_pt(6069, -0.082668000000000005, 0.18101700000000001, 0.41657894746307178)
p:create_pt(6071, -0.087165999999999993, 0.171073, 0.41873248315486722)
p:create_pt(3017, -0.070798728055994117, 0.17785276093100011, 0.39000000000000001)
p:create_pt(6073, -0.059331000000000002, 0.18260299999999999, 0.42231595730959332)
p:create_pt(5343, -0.059331000000000002, 0.18260299999999999, 0.43427292389069116)
p:create_pt(3025, -0.082668000000000005, 0.18101700000000001, 0.4431578949261436)
p:create_pt(3115, -0.087165999999999993, 0.171073, 0.44426672030792835)
p:create_pt(3001, -0.059331000000000002, 0.18260299999999999, 0.44954032886126616)
p:create_pt(5730, -0.056064999999999997, 0.190939, 0.45388955146064924)
p:create_pt(2317, -0.082668000000000005, 0.18101700000000001, 0.495)
p:create_pt(2316, -0.087165999999999993, 0.171073, 0.495)
p:create_pt(2375, -0.059331000000000002, 0.18260299999999999, 0.495)
p:create_pt(5732, -0.057755000000000001, 0.19669600000000001, 0.4553732680961276)
p:create_pt(2941, -0.075860999999999998, 0.19581899999999999, 0.49526724664035882)
p:create_pt(3471, -0.082668000000000005, 0.18101700000000001, 0.51466422186830929)
p:create_pt(3479, -0.087165999999999993, 0.171073, 0.51595628547624184)
p:create_pt(3478, -0.059331000000000002, 0.18260299999999999, 0.52967242672092041)
p:create_pt(4188, -0.08516, 0.186475, 0.53808270521141532)
p:create_pt(4410, -0.082668000000000005, 0.18101700000000001, 0.53883910946042302)
p:create_pt(5569, -0.087165999999999993, 0.171073, 0.5392387203274398)
p:create_pt(5102, -0.059331000000000002, 0.18260299999999999, 0.55627106333532783)
p:create_pt(5284, -0.082668000000000005, 0.18101700000000001, 0.55412933209531723)
p:create_pt(5267, -0.087165999999999993, 0.171073, 0.5526258831837616)
p:create_pt(4253, -0.075860999999999998, 0.19581899999999999, 0.56493566130467054)
p:create_pt(4566, -0.082668000000000005, 0.18101700000000001, 0.56941955473021144)
p:create_pt(5978, -0.087165999999999993, 0.171073, 0.56753034073325481)
p:create_pt(224, -0.059331000000000002, 0.18260299999999999, 0.59999999999999998)
p:create_pt(2969, -0.0707989073250265, 0.17785268667298165, 0.59999999999999998)
p:create_pt(3219, -0.068319000000000005, 0.17407300000000001, 0.4478901788158261)
p:create_pt(5342, -0.056064999999999997, 0.190939, 0.4425)

p:draw_subface(3017, 6071, 6069)
p:draw_subface(6073, 3017, 6069)
p:draw_subface(5343, 6073, 6069)
p:draw_subface(3025, 5343, 6069)
p:draw_subface(6071, 3025, 6069)
p:draw_subface(6071, 3115, 3025)
p:draw_subface(3001, 5343, 3025)
p:draw_subface(5730, 3001, 3025)
p:draw_subface(2317, 5730, 3025)
p:draw_subface(3115, 2317, 3025)
p:draw_subface(3115, 2316, 2317)
p:draw_subface(2375, 3001, 5730)
p:draw_subface(5732, 2375, 5730)
p:draw_subface(2317, 5732, 5730)
p:draw_subface(2941, 5732, 2317)
p:draw_subface(3471, 2941, 2317)
p:draw_subface(2316, 3471, 2317)
p:draw_subface(2941, 2375, 5732)
p:draw_subface(3471, 2375, 2941)
p:draw_subface(2316, 3479, 3471)
p:draw_subface(3478, 2375, 3471)
p:draw_subface(4188, 3478, 3471)
p:draw_subface(4410, 4188, 3471)
p:draw_subface(3479, 4410, 3471)
p:draw_subface(4410, 3478, 4188)
p:draw_subface(3479, 5569, 4410)
p:draw_subface(5102, 3478, 4410)
p:draw_subface(5284, 5102, 4410)
p:draw_subface(5267, 5284, 4410)
p:draw_subface(5569, 5267, 4410)
p:draw_subface(4253, 5102, 5284)
p:draw_subface(4566, 4253, 5284)
p:draw_subface(5978, 4566, 5284)
p:draw_subface(5267, 5978, 5284)
p:draw_subface(4566, 5102, 4253)
p:draw_subface(224, 5102, 4566)
p:draw_subface(2969, 224, 4566)
p:draw_subface(5978, 2969, 4566)
p:draw_tet(3001, 5730, 3025, 3219)
p:draw_tet(3001, 5730, 3025, 5342)
p:set("draw", "geom")
rnd=glvCreate(0, 0, 500, 500, "TetView")
p:plot(rnd)
glvWait()