import par
import math

width_p = 50
pml_width = 25.e-6
omega = 2 * par.pai * par.Hz
permitt = 8.8541878128e-12
R = 1e-5
sigma_max = (3 * permitt * par.c * par.n2) / (2 * pml_width) * math.log(1 / R)
p = []
pml_sigma = []
pml_s = []
for i in reversed(range(width_p, 1)):
    for k in range(width_p/2):
        p[i] = 25e-9 / i
        pml_sigma[i] = sigma_max * (pow((p[i] / pml_width), 2))
        pml_s[i] = 1 - (pml_sigma[i] /
                        (omega * permitt * pow(par.n2, 2)))*1j
