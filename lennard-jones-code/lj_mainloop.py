# NVT 20

import numpy as np
import random as rd

#............................................................................#

loc_bs, vel = initiate_loc_vel_acc(ndim, npart, boxsize)

loc_bs, center_of_mass = ref_to_com(ndim, npart, loc_bs)

#............................................................................#

# Main loop 
for step in range(1, nscans + 1):
    
    # Application of periodic boundary conditions (PBC)
    loc_bs = pbc(loc_bs)
    
    # Calculate acceleration, virial (for pressure), and potential energy
    acc, virial, pot_en = compute_force(ndim, npart, loc_bs, 
                                        boxsize, eps, sig, r_cutoff, 
                                        u_at_cutoff, du_at_cutoff)
    
    # Update the location using velocity Verlet algorithm
    loc_bs = update_loc(loc_bs, vel, acc, dt)

    # Calculate instant tempeature and scaling factor (for temperature)
    average_kin_en, instant_temp, scaling_factor = calculate_temp(ndim, 
                                                                  npart, vel, 
                                                                  boxsize, 
                                                                  reference_temp)
    # Is step is 20% of the scans, scale the velocities
    if step <= (nscans * 0.2):
        vel = update_vel(vel, scaling_factor, dt, acc)
    
        # Again, calculate instant temperature
        average_kin_en, instant_temp, scaling_factor = calculate_temp(ndim, 
                                                                      npart, vel, 
                                                                      boxsize, 
                                                                      reference_temp)
    
    # If step is beyond 20% of the scans, don't scale the velocities
    else:
        vel = update_vel_unscaled(vel, dt, acc)
        
        average_kin_en, instant_temp, scaling_factor = calculate_temp(ndim, 
                                                                      npart, vel, 
                                                                      boxsize, 
                                                                      reference_temp)
    
    virial = - (virial / ndim)
    average_pot_en = sum(pot_en) / npart
    average_tot_en = average_kin_en + average_pot_en
    pressure = (rho * instant_temp) + virial / volume
    comp_factor = pressure * volume / (npart * instant_temp)
    
    # Monitor the change of various parameters
    step1 = step - 1
    
    check_temp20[step1] = instant_temp
    check_pressure20[step1] = pressure
    check_kin_en20[step1] = average_kin_en
    check_pot_en20[step1] = average_pot_en
    check_total_en20[step1] = average_kin_en + average_pot_en
    
    # Write the data in files
    f_tpzsc.write("%d %e %e %e %e\n" % (step, step * dt, instant_temp, 
                                        pressure, comp_factor))
    f_kpesc.write("%d %e %e %e %e\n" % (step, step * dt, average_kin_en, 
                                        average_pot_en, average_tot_en))
    f_xyzsc.write("%d\n\n" % (step))
    
    for iatom in range(npart):
        f_xyzsc.write("%e %e %e\n" % (loc_bs[0, iatom] * boxsize, 
                                      loc_bs[1, iatom] * boxsize, 
                                      loc_bs[2, iatom] * boxsize))
    
    # Use a condition to eliminate initial extraordinary values
    if step > ignore_steps:
        total_pressure += pressure
        total_kin_en += average_kin_en
        total_pot_en += average_pot_en
        total_temp += instant_temp
    
    # Allow monitoring of progress of iteration
    if step % (nscans * 0.1) == 0:
        print(">>>>>> Iteration # %d out of %d" %(step, nscans))

effective_scans = nscans - ignore_steps
print("\n>>>>>> Statistical Averages <<<<<<\n") #real-time progress of the compilation
print("Temperature = %f K " %(total_temp / effective_scans))
print("Kinetic energy = %f eps" %(total_kin_en / effective_scans))
print("Potential energy = %f eps" %(total_pot_en / effective_scans))
print("Total energy = %f eps" %((total_kin_en + total_pot_en) / effective_scans))
print("Pressure = %f eps/(sig^3) \n" %(total_pressure / effective_scans))
