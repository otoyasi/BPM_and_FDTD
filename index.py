import par
import sett
import neff
import math
import random
import numpy as np

# ///**屈折率分布を求める関数**////
# ///**光導波路を決定する関数**////
k0 = par.k0
y_connection_before = 50
y_connection_after = 300
P_MAX = sett.p_max
pml_pmax = sett.p_max_um
R_MAX = sett.r_max
cladd_width_left = sett.core_left
cladd_width_right = sett.core_right
n1 = par.n1
n2 = par.n2
n = np.zeros([R_MAX, pml_pmax])

print("導波路の種類を選んでください\n")
print("1:直線光導波路\n2:傾斜導波路\n3:Y分岐導波路\n4:テーパ構造導波路\n5:不均一導波路\n6:曲がり導波路\n")
# scanf("%d", & waveguide_type)
#waveguide_type = (input())
waveguide_type = "8"
# / 光導波路の屈折率分布を与える
#    // 直線光導波路
if "1" in waveguide_type:
    for r in range(90):
        for p in range(cladd_width_left):
            n[r, p] = n2
           #
        for p in range(cladd_width_left,cladd_width_right):
            n[r, p] = n1
        #    print("n[r,p]",r,p,n[r,p])

        for p in range(cladd_width_right,pml_pmax):
            n[r, p] = n2
            # // fprintf(waveguide_in_dat_file, "%9d   %9d   %9f\n", r, p, input_electric - 1)
    for r in range(90,R_MAX):
        for p in range(cladd_width_left):
            n[r, p] = 2
           #
        for p in range(cladd_width_left,cladd_width_right):
            n[r, p] = 2
        #    print("n[r,p]",r,p,n[r,p])

        for p in range(cladd_width_right,pml_pmax):
            n[r, p] = 2
   # // 傾斜光導波路
elif "2" in waveguide_type:
    print("角度を入力してください\n")
    angle = int(input())
    be = (n1 * k0 * math.cos(angle * math.pi / 180.0)) / k0
    for r in range(R_MAX):

        for p in range(p, cladd_width_left-i):
            n[r, p] = n2

        for p in range(cladd_width_left-i, cladd_width_right-i):
            n[r, p] = n1

        for p in range(p-i, pml_pmax):
            n[r, p] = n2

        if r < R_MAX:

            i = i + 0.02

  #  // Y分岐導波路
elif "3" in waveguide_type:
    print("角度を入力してください\n")
    angle = int(input())
    be = (n1 * k0 * math.cos(angle * math.pi / 180.0)) / k0
    for r in range(y_connection_before):
        for p in range(p, cladd_width_left):
            n[r, p] = n2

        for cladd_width_left in range(cladd_width_left, cladd_width_right):
            n[r, p] = n1

        for cladd_width_right in range(pml_pmax):
            n[r, p] = n2

    for r in range(y_connection_after):
        for p in range(cladd_width_left-i):
            n[r, p] = n2

        for p in range(cladd_width_right+i):
            n[r, p] = n1

        for p in range(p+i, pml_pmax):
            n[r, p] = n2

        if r < y_connection_after:
            i = i + 0.1

   # // Y分岐左側の屈折率設定
    for r in range(R_MAX):

        for p in range(cladd_width_left-i):
            n[r, p] = n2

            # fprint(waveguide_in_dat_file, "%9d   %9d   %9f\n", r, p, input_electric - 1)

        for p in range(cladd_width_right-i):
            n[r, p] = n1

           # // fprintf(waveguide_in_dat_file, "%9d   %9d  non\n", r, p)
        for p in range(p-i, pml_pmax):
            n[r, p] = n2

            # // fprintf(waveguide_in_dat_file, "%9d   %9d  %9f\n", r, p, input_electric - 1)

    #    // Y分岐導波路右側
        for p in range(pml_pmax/2, cladd_width_left+i):

            n[r, p] = n2

        for p in range(cladd_width_left+i, cladd_width_right+i):
            n[r, p] = n1

        for p in range(cladd_width_right+i, pml_pmax):
            n[r, p] = n2

        if r < R_MAX:
            i = i + 0.2

 #   // テーパ構造導波路
elif "4" in waveguide_type:
    for r in range(10):
        for p in range(cladd_width_left):
            n[r, p] = n2

        for p in range(cladd_width_left, cladd_width_right):
            n[r, p] = n1

        for p in range(cladd_width_right, pml_pmax):
            n[r, p] = n2

    for r in range(10, R_MAX):
        for p in range(cladd_width_left-i):

            n[r, p] = n2

        for p in range(cladd_width_left-i, cladd_width_right-i):
            n[r, p] = n1

        for p in range(cladd_width_right+i, pml_pmax):
            n[r, p] = n2

        if r < R_MAX:
            i = i + 0.02

   # // 不均一光導波路
elif "5" in waveguide_type:
    for r in range(R_MAX):
        randam_one_to_three = (random.random() % 5) - 2

        for p in range(cladd_width_left-randam_one_to_three):
            n[r, p] = n2

          #  // fprintf(waveguide_in_dat_file, "%9d   %9d   %9f\n", r, p, input_electric)

        for p in range(cladd_width_left-randam_one_to_three, cladd_width_right+randam_one_to_three):
            n[r, p] = n1

           # // fprintf(waveguide_in_dat_file, "%9d   %9d   non\n", r, p)

        for p in range(cladd_width_right+randam_one_to_three, pml_pmax):
            n[r, p] = n2


elif "6" in waveguide_type:
    R = 5e-3
    for r in range(R_MAX):

        for p in range(cladd_width_left):
            n[r, p] = n2 * (1 + (p * par.dx - 625.5 * par.dx) / par.R)

          #  // fprintf(waveguide_in_dat_file, "%9d   %9d   %9f\n", r, p, n[r,p])

        for p in range(cladd_width_left, cladd_width_right):
            n[r, p] = n1 * (1 + (p * par.dx - 625.5 * par.dx) / par.R)

           # // fprintf(waveguide_in_dat_file, "%9d   %9d   %f\n", r, p, n[r,p])

        for p in range(cladd_width_right, pml_pmax):
            n[r, p] = n2 * (1 + (p * par.dx - 625.5 * par.dx) / par.R)


elif "7" in waveguide_type:

    for r in range(R_MAX):

        for p in range(cladd_width_left):
            n[r, p] = n2 * \
                (1 + ((1251 - p) * par.dx - 625.5 * par.dx) / par.R)

        for p in range(cladd_width_left, cladd_width_right):
            n[r, p] = n1 * \
                (1 + ((1251 - p) * par.dx - 625.5 * par.dx) / par.R)

        for p in range(cladd_width_right, pml_pmax):
            n[r, p] = n2 * \
                (1 + ((1251 - p) * par.dx - 625.5 * par.dx) / par.R)

elif "8" in waveguide_type:
    for p in range(sett.p_max_um):
        for r in range(cladd_width_left):
            n[r, p] = n2
           #
        for r in range(cladd_width_left,cladd_width_right):
            n[r, p] = n1
        #    print("n[r,p]",r,p,n[r,p])

        for r in range(cladd_width_right,pml_pmax):
            n[r, p] = n2
            # // fprintf(waveguide_in_dat_file, "%9d   %9d   %9f\n", r, p, input_electric - 1)



print("参照屈折率==%lf\n", neff.value)
