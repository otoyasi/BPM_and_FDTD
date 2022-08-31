from pickle import NONE
import par
import math
import sett


mode = "TE0"
neffmode = ""
a = (par.n2 * par.n2 - par.n2 * par.n2) / (par.n1 * par.n1 - par.n2 * par.n2)
c = par.n2 * par.n2 / par.n1 / par.n1
d = c - a * (1.0 - c)
V = par.k0 * par.dx * sett.core_width * \
    math.sqrt(par.n1 * par.n1 - par.n2 * par.n2)
neffmode = mode[2:3]
modeNo = float(neffmode)
if mode[:2] == "TE":
    modeID = 1
    pcf = pfs = 1.0
else:
    modeID = 2
    pcf = d
    pfs = c


# /**参照屈折率neffを計算する関数**/
# def refrection_neff_index(par, sett, neff):
def eigenFunction_y1(a, b, c, d, b1, modeNo, V):

    if modeID == 1:
        value_1 = V * math.sqrt(1.0 - b1) - math.atan2(math.sqrt(a + b1), math.sqrt(
            1.0 - b1)) - math.atan2(math.sqrt(b1), math.sqrt(1.0 - b1)) - modeNo * math.pi
    else:
        value_1 = V * math.sqrt(1.0 - b1) - math.atan2(math.sqrt(a + b1), math.sqrt(1.0 - b1) * d) - math.atan2(math.sqrt(b1), math.sqrt(1.0 - b1) * c) - modeNo * math.pi
    return value_1

def eigenFunction_y2(a, b, c, d, b2, modeNo, V):

    if modeID == 1:
        value_2 = V * math.sqrt(1.0 - b2) - math.atan2(math.sqrt(a + b2), math.sqrt(
            1.0 - b2)) - math.atan2(math.sqrt(b2), math.sqrt(1.0 - b2)) - modeNo * math.pi
    else:
        value_2 = V * math.sqrt(1.0 - b2) - math.atan2(math.sqrt(a + b2), math.sqrt(1.0 - b2) * d) - math.atan2(math.sqrt(b2), math.sqrt(1.0 - b2) * c) - modeNo * math.pi
    return value_2

def eigenFunction(a, b, c, d, x, modeNo, V):

    if modeID == 1:
        value = V * math.sqrt(1.0 - b) - math.atan2(math.sqrt(a + b), math.sqrt(
            1.0 - b)) - math.atan2(math.sqrt(b), math.sqrt(1.0 - b)) - modeNo * math.pi
    else:
        value = V * math.sqrt(1.0 - b) - math.atan2(math.sqrt(a + b), math.sqrt(1.0 - b) * d) - math.atan2(math.sqrt(b), math.sqrt(1.0 - b) * c) - modeNo * math.pi
    return value

b = 1.0
b1 = 0.00001
b2 = 0.99999
y1 = eigenFunction_y1(a, b, c, d, b1, modeNo, V)
y2 = eigenFunction_y2(a, b, c, d, b2, modeNo, V)
print("y1==",y1)
print("y2==",y2)

# if y1 * y2 > 0.0:
#   return -1
for i in range(1000):
    b = (b1 + b2) / 2.0
    y = eigenFunction(a, b, c, d, b, modeNo, V)
    if y * y1 > 0.0:
        b1 = b
    else:
        b2 = b
    if b2 - b1 < 1.0e-14:
        break
neff_index = math.sqrt(b * (par.n1 * par.n1 - par.n2 * par.n2) + par.n2 * par.n2)
value = round(neff_index * 1e8) / 1e8
print("neff ==" ,value)
