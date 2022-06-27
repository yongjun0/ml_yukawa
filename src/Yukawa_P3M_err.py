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
kappa = 3.0

tol = 1.e-6
tol = tol/np.sqrt(2)
#============================
alpha_ewald = 0.4875e9   # for PM err
alpha = alpha_ewald
rc = 5.604e-9            # for PP err
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
pm_err = lambda x: (sub.gf_opt(MGrid, aliases, Lv, cao, N, kappa, x, rc, fourpie0, aws, box_volume, flag="PM"))
#pm_err = lambda x: (sub.gf_opt(MGrid, aliases, Lv, cao, N, kappa, x, rc, fourpie0, aws, box_volume, flag))
#pm_root = optimize.newton(pm_err, G_ew, maxiter=150)

err_min = 1000
nloop = 1000
alpha_root = 10000
if(1):
    for i in range(nloop):
      #if(i%10 == 0):
      #    print(i, err_min)
      alpha_temp = 0.5*(i)/nloop*alpha + 1.0*alpha
      
      err = pm_err(alpha_temp)
      rel_err = abs(err-tol)/tol
       
      if(rel_err <= err_min):
        err_min = rel_err
        alpha_root = alpha_temp

      if(i%10 == 0):
          print(f"{i}, {err:5.4e}, {rel_err:5.4e}, {alpha_root:5.4e}")

      if(err_min <= 1e-2):
          print(f"!!!: {i}, {err:5.4e}, {rel_err:5.4e}, {alpha_root:5.4e}")
          break

#alpha_root = 4.9067e+08
#print("!!!pm err1 = ", pm_err(alpha_root))
# PP error
print("Now PP")
#flag = "PP" 
pp_err = lambda y: (sub.gf_opt(MGrid, aliases, Lv, cao, N, kappa, alpha_root, y, fourpie0, aws, box_volume, flag="PP"))

err_min = 1000
nloop = 100
for i in range(nloop):
  #if(i%1 == 0):
  rc_temp = 1.0*(i)/nloop*rc + 0.7*rc
  err = pp_err(rc_temp)
  rel_err = abs(err-tol)/tol
  if(rel_err <= err_min):
    err_min = rel_err
    rc_root = rc_temp

  if(i%10 == 0):
      print(f"{i}, {err:5.4e}, {rel_err:5.4e}, {rc_root:5.4e}")

  if(err_min <= 1e-2):
      print(f"!!!: {i}, {err:5.4e}, {rel_err:5.4e}, {rc_root:5.4e}")
      break

print("====================================")
print(f"alpha = {alpha_root:5.4e}")
print(f"rc = {rc_root:5.4e}")

PM_err = pm_err(alpha_root)
PP_err = pp_err(rc_root)
err = np.sqrt(PM_err**2 + PP_err**2)

print(f"PM err = {PM_err:5.4e}")
print(f"PP err = {PP_err:5.4e}")
print(f"total err = {err:5.4e}")
