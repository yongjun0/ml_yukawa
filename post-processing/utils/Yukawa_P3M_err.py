""" 
Module for Yukawa pot P3M err calculation.
Yongjun Choi
choiyj@msu.edu
"""
import scipy.constants as const
import numpy as np
import numba as nb
import math as mt
import sys
from scipy import optimize

import Yukawa_P3M_err_sub as sub

#============================
"""
input data
"""
kappa = 2.1             # Yukawa screening parameter

tol = 1.e-6             # force error tolerance
tol = tol/np.sqrt(2)    # Because we have PM and PP force errors. tol = sqrt(PM_err**2 + PP_err**2)
#============================
alpha_min = 1e8         # small enough than kappa = 0 case.
alpha_max = 1e9         # large enough than kappa = 3 case.

rc_min = 5.3e-11        # Bohr radius. Small enough
rc_max = 1e-7           # = 0.1 um. So I think this is large enough.

alpha_ewald = (alpha_min + alpha_max)/2.0   # any number btw min and max
rc = (rc_min + rc_max)/2.0                  # any number btw min and max
alpha = alpha_ewald     # alias
#============================
"""
Initial condition. Getting from a Sarkas input file.
Refer to the manual for the meaning of each variable.
"""
num_density = 1.0e26
mass =  1.672621898e-27
Zi = 1
charge = Zi* 1.602176634e-19
MGrid = np.array([32, 32, 32])
aliases = np.array([3, 3, 3])
cao = 6.0
G_ew = alpha_ewald
N = 10000               # number of particles
#============================
eps0 = const.epsilon_0                              # permissivity
fourpie0 = 4.0*np.pi*eps0                            
twopi = 2.0*np.pi

wp = np.sqrt(charge**2*num_density/(mass*eps0))     # plasma frequency square
aws = (3.0/(4.0*np.pi*num_density))**(1./3.)        # ion sphere radius
L = aws*(4.0*np.pi*N/3.0)**(1.0/3.0)                # box length
Lv = np.array([L, L, L])                            # box length vector
box_volume = L**3
kappa /=aws                                         # de-normalize kappa. now kappa means 1/(screening distance)

# PM error
print("Now PM")
pm_err = lambda x: (sub.gf_opt(MGrid, aliases, Lv, cao, N, kappa, x, rc, fourpie0, aws, box_volume, flag="PM"))

# root finding for alpha
alpha_root = sub.root_finding(alpha_min, alpha_max, alpha, pm_err, tol)

# PP error
print("Now PP")
pp_err = lambda y: (sub.gf_opt(MGrid, aliases, Lv, cao, N, kappa, alpha_root, y, fourpie0, aws, box_volume, flag="PP"))
# root finding for rc 
rc_root = sub.root_finding(rc_max, rc_min, rc, pp_err, tol)

print("====================================")
print(f"alpha = {alpha_root:5.4e}")
print(f"rc = {rc_root:5.4e}")

PM_err = pm_err(alpha_root)
PP_err = pp_err(rc_root)
err = np.sqrt(PM_err**2 + PP_err**2)

print(f"PM err = {PM_err:5.4e}")
print(f"PP err = {PP_err:5.4e}")
print(f"total err = {err:5.4e}")
