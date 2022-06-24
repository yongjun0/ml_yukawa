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
#============================
kappa = 2.2
tol = 1.e-6
tol = tol/np.sqrt(2)
#============================
alpha_ewald = 0.4058e9   # for PM err
rc = 7.604e-9            # for PP err
#rc = 8.e-9            # for PP err
#============================
num_density = 1.0e26
mass =  1.672621898e-27
Zi = 1
charge = Zi* 1.602176634e-19
#temperature = 38.11
MGrid = np.array([32, 32, 32])
aliases = np.array([3, 3, 3])
cao = 6.0
G_ew = alpha_ewald
N = 10000
#============================
eps0 = const.epsilon_0
fourpie0 = 4.0*np.pi*eps0
twopi = 2.0*np.pi
#beta_i = 1.0/(const.Boltzmann*temperature)

wp2 = charge**2*num_density/(mass*eps0)
wp = np.sqrt(wp2)
aws = (3.0/(4.0*np.pi*num_density))**(1./3.)
L = aws*(4.0*np.pi*N/3.0)**(1.0/3.0)      # box length
Lv = np.array([L, L, L])              # box length vector
box_volume = L**3
kappa /=aws

# PM error
print("Now PM")
flag = "PM" 
pm_err = lambda x: (sub.gf_opt(MGrid, aliases, Lv, cao, N, kappa, x, rc, fourpie0, aws, box_volume, flag)-tol)
pm_root = optimize.newton(pm_err, G_ew)
print(f"{pm_root:8.4e}, {G_ew:8.4e}")
print(f"{pm_err(pm_root)+tol}, {pm_err(G_ew)+tol}")

# PP error
print("Now PP")
flag = "PP" 
pp_err = lambda x: (sub.gf_opt(MGrid, aliases, Lv, cao, N, kappa, pm_root, x, fourpie0, aws, box_volume, flag)-tol)

err_min = 1000
nloop = 100
for i in range(nloop):
  if(i%10 == 0):
      print(i, err_min)
  rc_temp = 1.0*(i)/nloop*rc + 0.7*rc
  err = pp_err(rc_temp)
  if(err < err_min):
    err_min = err
    rc_min = rc_temp

  if(err_min < 0):
      break

print(rc_min, err_min)
print(pp_err(rc_min)+tol, tol)
print(f"rc = {rc_min:8.4e}, alpha={G_ew:8.4e}")
