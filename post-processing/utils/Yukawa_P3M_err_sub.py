"""
Sub module for Yukawa pot P3M err calculation.
Yongjun Choi
choiyj@msu.edu
"""
import scipy.constants as const
import numpy as np
import numba as nb
import math as mt
import sys

#@nb.njit
def root_finding(x_min, x_max, x, func, tol):
    """
    root finding using a binary bifercation.
    Numba does not support f-statement yet.
    Most work is func call, so the performance is alsmost same without numba call.
    """
    nloop = 100
    err_low = func(x_min)
    err_up = func(x_max)

    for i in range(nloop):
        err = func(x)
        rel_err = abs(err-tol)/tol
        if(err < tol):
            x_min = x
            x = (x+x_max)/2.

        if(err >= tol):
            x_max = x
            x = (x+x_min)/2.

        print(f"{i}, {err:5.4e}, {rel_err:5.4e}, {x_min:5.4e}, {x:5.4e}, {x_max:5.4e}")

        if(rel_err <= 1e-3):
            print(f"!!!: {i}, {err:5.4e}, {rel_err:5.4e}, {x_min:5.4e}, {x:5.4e}, {x_max:5.4e}")
            break

    return x


#########################################################
@nb.njit
def gf_opt(MGrid, aliases, BoxLv, p, N, kappa, Gew, rcut, fourpie0, aws, box_volume, flag):
    """ 
    Calculate the Optimized Green Function given by eq.(22) of Ref. [2]_.

    Parameters
    ----------
    MGrid : array
        number of mesh points in x,y,z

    aliases : array
        number of aliases in each direction

    BoxLv : array
        Length of simulation's box in each direction

    p : int
        charge assignment order (CAO)

    N : int
        number of particles

    pot_matrix : array
        Potential matrix. It contains screening parameter and Ewald parameter. See potential matrix above.

    rcut : float
        Cutoff distance for the PP calculation

    fourpie0 : float
        Potential factor.

    Returns
    -------
    G_k : array_like
        optimal Green Function

    kx_v : array_like
       array of reciprocal space vectors along the x-axis

    ky_v : array_like
       array of reciprocal space vectors along the y-axis

    kz_v : array_like
       array of reciprocal space vectors along the z-axis

    PM_err : float
        Error in the force calculation due to the optimized Green's function. eq.(28) of Ref. [2]_

    PP_err : float
        Error in the force calculation due to the distance cutoff. eq.(30) of Ref. [1]_
   
    References
    ----------
    .. [2] `H.A. Stern et al. J Chem Phys 128, 214006 (2008) <https://doi.org/10.1063/1.2932253>`_
    """
    #kappa = pot_matrix[0,0,0] #params.Potential.matrix[0,0,0]
    #Gew = pot_matrix[-1,0,0] #params.Potential.matrix[3,0,0]
    rcut2 = rcut*rcut
    mx_max = aliases[0] #params.P3M.mx_max
    my_max = aliases[1] # params.P3M.my_max
    mz_max = aliases[2] #params.P3M.mz_max
    Mx = MGrid[0] #params.P3M.Mx
    My = MGrid[1] #params.P3M.My
    Mz = MGrid[2] #params.P3M.Mz
    Lx = BoxLv[0] #params.Lx
    Ly = BoxLv[1] #params.Ly
    Lz = BoxLv[2] #params.Lz
    hx = Lx/float(Mx)
    hy = Ly/float(My)
    hz = Lz/float(Mz)

    kappa_sq = kappa*kappa
    Gew_sq = Gew*Gew

    G_k = np.zeros((Mz,My,Mx))
    
    if np.mod(Mz,2) == 0:
        nz_mid = Mz/2
    else:
        nz_mid = (Mz-1)/2
    
    if np.mod(My,2) == 0:
        ny_mid = My/2
    else:
        ny_mid = (My-1)/2
    
    if np.mod(Mx,2) == 0:
        nx_mid = Mx/2
    else:
        nx_mid = (Mx-1)/2
        
    nx_v = np.arange(Mx).reshape((1,Mx))
    ny_v = np.arange(My).reshape((My,1))
    nz_v = np.arange(Mz).reshape((Mz,1,1))
    
    kx_v = 2.0*np.pi*(nx_v - nx_mid)/Lx
    ky_v = 2.0*np.pi*(ny_v - ny_mid)/Ly
    kz_v = 2.0*np.pi*(nz_v - nz_mid)/Lz
    
    PM_err = 0.0
    
    for nz in range(Mz):
        nz_sh = nz-nz_mid
        kz = 2.0*np.pi*nz_sh/Lz
        
        for ny in range(My):
            ny_sh = ny-ny_mid
            ky = 2.0*np.pi*ny_sh/Ly
            
            for nx in range(Mx):
                nx_sh = nx-nx_mid
                kx = 2.0*np.pi*nx_sh/Lx
                           
                k_sq = kx*kx + ky*ky + kz*kz
                
                if k_sq != 0.0:
                
                    U_k_sq = 0.0
                    U_G_k = 0.0

                    # Sum over the aliases
                    for mz in range(-mz_max,mz_max+1):
                        for my in range(-my_max,my_max+1):
                            for mx in range(-mx_max,mx_max+1):
                                  
                                kx_M = 2.0*np.pi*(nx_sh + mx*Mx)/Lx
                                ky_M = 2.0*np.pi*(ny_sh + my*My)/Ly
                                kz_M = 2.0*np.pi*(nz_sh + mz*Mz)/Lz
                            
                                k_M_sq = kx_M**2 + ky_M**2 + kz_M**2
                                
                                if kx_M != 0.0:
                                    U_kx_M = np.sin(0.5*kx_M*hx)/(0.5*kx_M*hx)
                                else:
                                    U_kx_M = 1.0

                                if ky_M != 0.0:
                                    U_ky_M = np.sin(0.5*ky_M*hy)/(0.5*ky_M*hy)
                                else: 
                                    U_ky_M = 1.0
                                    
                                if kz_M != 0.0:
                                    U_kz_M = np.sin(0.5*kz_M*hz)/(0.5*kz_M*hz)
                                else:
                                    U_kz_M = 1.0
                                
                                U_k_M = (U_kx_M*U_ky_M*U_kz_M)**p
                                U_k_M_sq = U_k_M*U_k_M
                                
                                G_k_M = np.exp(-0.25*(kappa_sq + k_M_sq)/Gew_sq)/(kappa_sq + k_M_sq)/fourpie0
                                
                                k_dot_k_M = kx*kx_M + ky*ky_M + kz*kz_M

                                U_G_k += (U_k_M_sq * G_k_M * k_dot_k_M)
                                U_k_sq += U_k_M_sq
                                
                    # eq.(22) of Ref.[2]_
                    G_k[nz,ny,nx] = U_G_k/((U_k_sq**2)*k_sq)
                    Gk_hat = np.exp(-0.25*(kappa_sq + k_sq)/Gew_sq)/(kappa_sq + k_sq)/fourpie0      

                    # eq.(28) of Ref.[2]_
                    PM_err = PM_err + Gk_hat*Gk_hat*k_sq - U_G_k**2/((U_k_sq**2)*k_sq)

    PP_err = 2.0/np.sqrt(Lx*Ly*Lz)*np.exp(-0.25*kappa_sq/Gew_sq)*np.exp(-Gew_sq*rcut2)/np.sqrt(rcut)/fourpie0
    PM_err = np.sqrt(PM_err)/Lx
    PP_err *= np.sqrt(N)*aws**2*fourpie0
    PM_err *= np.sqrt(N)*aws**2*fourpie0/box_volume**(2./3.)
    if(flag == "PM"):
        err = PM_err
    if(flag == "PP"):
        err = PP_err

    return err

