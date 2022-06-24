""" 
Module for handling Coulomb interaction
"""
import scipy.constants as const
import numpy as np
import numba as nb
import math as mt
import sys
import gf_opt

num_density = 1.0e26
mass =  1.672621898e-27
Zi = 1
charge = Zi*1.6e-19
temperature = 38.11
rc = 7.604e-9
kappa = 2.1
MGrid = np.array([32,32,32])
aliases = np.array([3, 3, 3])
cao = 6.0
alpha_ewald = 0.4058e9
G_ew = alpha_ewald
N = 10000

eps0 = const.epsilon_0
fourpie0 = 4.0*np.pi*eps0
twopi = 2.0*np.pi
beta_i = 1.0/(const.Boltzmann*temperature)

wp2 = charge**2*num_density/(mass*eps0)
wp = np.sqrt(wp2)
aws = (3.0/(4.0*np.pi*num_density))**(1./3.)
L = aws*(4.0*np.pi*N/3.0)**(1.0/3.0)      # box length
Lv = np.array([L, L, L])              # box length vector
box_volume = L**3
kappa /=aws
print('--------------------')
print(kappa)
print(aliases)
print(Lv)
print(cao)
print(N)
print(G_ew)
print(rc)
print(fourpie0)
print(aws)
print('--------------------')
PM_err, PP_err = gf_opt.gf_opt(MGrid, aliases, Lv, cao, N, kappa, G_ew, rc, fourpie0)

PP_err *= np.sqrt(N)*aws**2*fourpie0
PM_err *= np.sqrt(N)*aws**2*fourpie0/box_volume**(2./3.)
F_err = np.sqrt(PM_err**2 + PP_err**2)

print(f"PP = {PP_err}")
print(f"PM = {PM_err}")
print(f"F = {F_err}")


