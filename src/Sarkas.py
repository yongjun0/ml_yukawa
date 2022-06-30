"""
                       Sarkas

A code that executes molecular dynamics simulation for a Yukawa plasma using
the efficient Particle-Particle-Particle-Mesh algorithm for force computation.
The code constitutes a number of functions that are in separate files.

Developed by the research group of:
Professor Michael S. Murillo
murillom@msu.edu
Dept. of Computational Mathematics, Science, and Engineering,
Michigan State University
"""

# Python modules
import numpy as np
import time
import sys
import os

# Importing MD modules, non class
import S_thermostat as thermostat

# import MD modules, class
from S_integrator import Integrator, calc_pot_acc
from S_particles import Particles
from S_verbose import Verbose
from S_params import Params
from S_checkpoint import Checkpoint

time_stamp= np.zeros(10)
its = 0
time_stamp[its] = time.time(); its += 1
input_file = sys.argv[1]
params = Params()
params.setup(input_file)                # Read initial conditions and setup parameters


###############################
#params.P3M.G_ew, params.Potential.rc = sub.Yukawa_err(params)
#print("=================")
#print(f"alpha1 = {params.P3M.G_ew:5.4e}")
#print(f"rc = {params.Potential.rc:5.4e}")
#sys.exit()
###############################



integrator = Integrator(params)
checkpoint = Checkpoint(params)         # For restart and pva backups.

verbose = Verbose(params)
#######################
Nt = params.Control.Nstep    # number of time steps
N = params.total_num_ptcls
Ndump = int(Nt/params.Control.dump_step)
print("Ndump = ", Ndump)
pos_all = np.zeros((Ndump+1, N, 3))
vel_all = np.zeros((Ndump+1, N, 3))


#######################
# Un-comment the following if you want to calculate n(q,t)
# n_q_t = np.zeros((params.Control.Nt, params.Nq, 3), dtype="complex128")
# initializing the wave vector vector qv
#qv = np.zeros(params.Nq)
#for iqv in range(0, params.Nq, 3):
#    iq = iqv/3.
#    qv[iqv] = (iq+1.)*params.dq
#    qv[iqv+1] = (iq+1.)*np.sqrt(2.)*params.dq
#    qv[iqv+2] = (iq+1.)*np.sqrt(3.)*params.dq

########################
# Print simulation's parameters to screen
if(params.Control.verbose):
    verbose.sim_setting_summary()   # simulation setting summary

# array for temperature, total energy, kinetic energy, potential energy
t_Tp_E_K_U = np.zeros((1, 5))

# Initializing particle positions and velocities
time_stamp[its] = time.time(); its += 1
ptcls = Particles(params)
ptcls.load()
ptcls.assign_attributes(params)

###############################
if(0):
    ptcls.pos[0, 0] = 0e-3
    ptcls.pos[0, 1] = 2e-9
    ptcls.pos[0, 2] = 0

    ptcls.vel[0, 0] = 1000
    ptcls.vel[0, 1] = 0
    ptcls.vel[0, 2] = 0

    ptcls.pos[1, 0] = 0
    ptcls.pos[1, 1] = 0
    ptcls.pos[1, 2] = 0

    ptcls.vel[1, 0] = 0
    ptcls.vel[1, 1] = 0
    ptcls.vel[1, 2] = 0

    m1 = ptcls.mass[0]
    m2 = ptcls.mass[1]
    q1 = ptcls.charge[0]
    q2 = ptcls.charge[1]
    eps0 = 8.8541878178176e-12
#b90 = q1*q2/(4*np.pi*eps0)*1/(mr*ptcls.vel[0, 0]**2)
#print("b90 = ", b90)
#ptcls.pos[0, 1] = 3*b90

    ptcls.species_name[0] = "P1"
    ptcls.species_name[1] = "P2"
################
# Calculate initial forces and potential energy
U = calc_pot_acc(ptcls,params)
print("U = ", U)



# Calculate initial kinetic energy and temperature
Ks, Tps = thermostat.calc_kin_temp(ptcls.vel,ptcls.species_num,ptcls.species_mass,params.kB)
K = np.ndarray.sum(Ks)
Tp = np.ndarray.sum(Tps)/params.num_species
            
E = K + U
thermostat.remove_drift(ptcls.vel,ptcls.species_num, ptcls.species_mass)
if ( params.Control.screen_output):
    #print("\nInitial: T = {:2.6e}, E = {:2.6e}, K = {:2.6e}, U = {:2.6e}".format(Tp, E, K, U) )
    print("\nInitial: T, E, K, U", Tp, E, K, U)

time_stamp[its] = time.time(); its += 1

# Un-comment the two following lines if you want to save Thermalization data for debugging
#f_output_E = open(params.Control.checkpoint_dir+"/"+"ThermEnergy_"+params.Control.fname_app + ".out", "w")
#f_xyz = open(params.Control.checkpoint_dir+"/"+"ThermPVA_"+params.Control.fname_app + ".xyz", "w")
if not (params.load_method == "restart"):
    if ( params.Control.screen_output):
        print("\n------------- Equilibration -------------")
    for it in range(params.Control.Neq):
        # Calculate the Potential energy and update particles' data
        U = integrator.update(ptcls,params)
        # Thermostate
        thermostat.Berendsen(ptcls,params,it)

        # Print Energies and Temperature to screen
        if (it % params.Control.dump_step == 0 and params.Control.screen_output):
            Ks, Tps = thermostat.calc_kin_temp(ptcls.vel,ptcls.species_num,ptcls.species_mass,params.kB)
            K = np.ndarray.sum(Ks)
            Tp = np.ndarray.sum(Tps)/params.num_species
            
            E = K + U
            if ( (it% params.Control.dump_step) == 0 and params.Control.screen_output):
                print("Equilibration: timestep {:6}, T = {:2.6e}, E = {:2.6e}, K = {:2.6e}, U = {:2.6e}".format(it, Tp, E, K, U) )

        # Un-comment the following if-statement if you want to save Thermalization data for debugging    
        #if params.Control.writexyz == 1:
        #    f_xyz.writelines("{0:d}\n".format(N))
        #    f_xyz.writelines("name x y z vx vy vz ax ay az\n")
        #    np.savetxt(f_xyz, np.c_[ptcls.species_name, ptcls.pos/params.Lx, ptcls.vel, ptcls.acc], 
        #        fmt="%s %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e")
    thermostat.remove_drift(ptcls.vel,ptcls.species_num, ptcls.species_mass)
pos_all[0] = ptcls.pos
vel_all[0] = ptcls.vel
checkpoint.dump(ptcls, 0)

time_stamp[its] = time.time(); its += 1
# Turn on magnetic field, if not on already, and thermalize
if (params.Magnetic.on == 1 and params.Magnetic.elec_therm == 1):
    params.Integrator.type == "Magnetic_Verlet"
    integrator = Integrator(params)

    if ( params.Control.screen_output):
        print("\n------------- Magnetic Equilibration -------------")
    
    for it in range(params.Magnetic.Neq_mag):
        # Calculate the Potential energy and update particles' data
        U = integrator.update(ptcls,params)
        # Thermostate
        thermostat.Berendsen(ptcls,params,it)

        # Print Energies and Temperature to screen
        if ( (it% params.Control.dump_step) == 0 and params.Control.screen_output):
            Ks, Tps = thermostat.calc_kin_temp(ptcls.vel,ptcls.species_num,ptcls.species_mass,params.kB)

            K = np.ndarray.sum(Ks)
            Tp = np.ndarray.sum(Tps)/params.num_species
            
            E = K + U
            print("Magnetic Equilibration: timestep {:6}, T = {:2.6e}, E = {:2.6e}, K = {:2.6e}, U = {:2.6e}".format(it, Tp, E, K, U) )

    thermostat.remove_drift(ptcls.vel,ptcls.species_num, ptcls.species_mass)
    # saving the 0th step
    checkpoint.dump(ptcls, 0)
    time_stamp[its] = time.time(); its += 1

# Close thermalization files. Un-comment if you want to save Thermalization data for debugging
#f_output_E.close()
#f_xyz.close()

# Open output files

if (params.load_method == "restart"):
    it_start = params.load_restart_step + 0
    f_output_E = open(params.Control.checkpoint_dir+"/"+"Energy_" + params.Control.fname_app + ".out", "a+")
    
    # Particles' positions scaled by Wigner-Seitz radius for Radial Distribution Function calculation
    #f_pos = open(params.Control.checkpoint_dir+"/"+"pos_" + params.Control.fname_app + ".out","a+")

    if (params.Control.writexyz == 1):
        # Particles' positions, velocities, accelerations for OVITO
        f_xyz = open(params.Control.checkpoint_dir+"/"+"pva_" + params.Control.fname_app + ".xyz", "a+")

else:
    it_start = 0
    f_output_E = open(params.Control.checkpoint_dir+"/"+"Energy_" + params.Control.fname_app + ".out", "w+")
    # Particles' positions scaled by Wigner-Seitz radius for Radial Distribution Function calculation
    #f_pos = open(params.Control.checkpoint_dir+"/"+"pos_" + params.Control.fname_app + ".out","w+")

    if (params.Control.writexyz == 1):
        # Particles' positions, velocities, accelerations for OVITO
        f_xyz = open(params.Control.checkpoint_dir+"/"+"pva_" + params.Control.fname_app + ".xyz", "w+")


ic_dump = 0
if ( params.Control.screen_output):
    print("\n------------- Production -------------")

for it in range(it_start, Nt):
    # Move the particles and calculate the potential
    U = integrator.update(ptcls,params)

    # Calculate Kinetic Energy and Temperature
    Ks, Tps = thermostat.calc_kin_temp(ptcls.vel,ptcls.species_num,ptcls.species_mass,params.kB)

    K = np.ndarray.sum(Ks)
    Tp = np.ndarray.sum(Tps)/params.num_species
            
    # Calculate the total Energy
    E = K + U

    # Print progress to screen
    if (it % params.Control.dump_step == 0 and params.Control.screen_output):
        print("Production: timestep, T, E, K, U = ", it, Tp, E, K, U)

    # Write Energies and Temperature to file
    t_Tp_E_K_U = np.array([params.Control.dt*it, Tp, E, K, U]).reshape(1,5)
    np.savetxt(f_output_E, t_Tp_E_K_U)
    f_output_E.flush()


    # Save particles' data for restart
    if ((it+1) % params.Control.dump_step == 0):
        ic_dump += 1
        print("it = ", it+1, ic_dump)
        pos_all[ic_dump] = ptcls.pos
        vel_all[ic_dump] = ptcls.vel
        #checkpoint.dump(ptcls, it+1)

    # Un-comment for n(q,t) calculation
    # Spatial Fourier transform for (Dynamic) Structure Factor
    #for iqv in range(params.Nq):
    #   q_p = qv[iqv]
    #   n_q_t[it, iqv, 0] = np.sum(np.exp(-1j*q_p*ptcls.pos[:, 0]))
    #   n_q_t[it, iqv, 1] = np.sum(np.exp(-1j*q_p*ptcls.pos[:, 1]))
    #   n_q_t[it, iqv, 2] = np.sum(np.exp(-1j*q_p*ptcls.pos[:, 2]))
    
    # Write particles' data to XYZ file for OVITO Visualization 
    if params.Control.writexyz == 1:
        #np.savetxt(f_pos,ptcls.pos/params.aws, fmt="%.6e %.6e %.6e")

        f_xyz.writelines("{0:d}\n".format(N))
        f_xyz.writelines("name x y z vx vy vz ax ay az\n")
        if(0):
            np.savetxt(f_xyz, np.c_[ptcls.species_name, \
                                    ptcls.pos/params.aws, \
                                    ptcls.vel/(params.wp*params.aws), \
                                    ptcls.acc/(params.aws*params.wp**2)], 
                                    fmt="%s %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e")
        np.savetxt(f_xyz, np.c_[ptcls.species_name, \
                                ptcls.pos, \
                                ptcls.vel, \
                                ptcls.acc], 
                                fmt="%s %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e")

# np.save("n_qt", n_q_t)
time_stamp[its] = time.time(); its += 1
# close output files
f_output_E.close()
#f_pos.close()

if (params.Control.writexyz):
    f_xyz.close()

# Print elapsed times to screen
if (params.Control.verbose):
    verbose.time_stamp(time_stamp)

file_name = params.Control.checkpoint_dir+"/"+"test"
np.savez(file_name, pos=ptcls.pos, vel=ptcls.vel)

# end of the code
