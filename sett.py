import pml
import math

p_max = 126 + pml.width_p
p_max_um = 127 + pml.width_p
print("pml.width_p == %d\n", pml.width_p)
r_max = 360
t_max = 10
range = (p_max_um + 1) * (r_max + 1) * (t_max + 1)
s_max = p_max_um * (r_max + 1)
core_width = 20
# HACK:コア幅が奇数の時は偶対象（正しい）,偶数の時は奇対称となり値が正しいか不明
if math.fmod(core_width, 2) == 0:

	core_cladd_boundary = core_width // 2
	core_left = (p_max // 2) - core_cladd_boundary + 1
	core_right = (p_max // 2) + core_cladd_boundary + 1
	print("偶数モード\n")
	print("pml.width_p == ", pml.width_p)
	print("cladd_left==%d_cladd_right==%d_core_cladd_boundary==%lf\n", core_left, core_right, core_cladd_boundary)

else:

	core_cladd_boundary = core_width / 2
	core_left = 62.5 - core_cladd_boundary + 1
	core_right = 62.5 + core_cladd_boundary + 1
	print("奇数モード\n")
	print("cladd_left==%d_cladd_right==%d_core_cladd_boundary==%lf\n", core_left, core_right, core_cladd_boundary)

