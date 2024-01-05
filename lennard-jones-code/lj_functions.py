def initiate_loc_vel_acc(ndim, npart, boxsize):
    '''
    Providing locations and velocities in all three axes from random
    distribution. Variables that have been used should be self-explanatory.

	For each particle, a 3D coordinate and 3D velocity vector have been 
 	picked up. Once done, the location has been copied to another variable
	and that has been divided by the box size to reduce the 3D coordinates in 
 	box length. 
    
    That copying is necessary since changing the numpy array in the future would
	not impact the arrays that were created in the past. To avoid any overlapping
 	or conflicts, taking a .copy() of the array is prudent.  

  	In this entire script, all the variables contain words that are related to
   	the processes that have been already carried out on the particle coordinates.
    '''
    loc = np.zeros((ndim, npart))
    vel = np.zeros((ndim, npart))
    
    for ncolumn in range(npart):
        for nrow in range(ndim):
            loc[nrow, ncolumn] = rd.uniform(0, boxsize)
            vel[nrow, ncolumn] = rd.uniform(0, 1)
            
    loc_bs = loc.copy() 
    
    loc_bs /= boxsize
    
    return loc_bs, vel

def ref_to_com(ndim, npart, loc_bs):
    '''
    Atom after atom, the center of mass(com) is updated by adding
    the location numbers along a single axis.
    
    nrow = 0, 1, 2 means X, Y, Z axes, respectively.
    
    loc_bs_com has been created to indicate that the referencing
    to com is done. As already mentioned, just assigning a new variable doesn't
    work, a .copy() would be necessary so change of one won't change the
    other.
    '''
    
    center_of_mass = np.zeros(ndim)

    for nrow in range(ndim):
        for ncolumn in range(npart):
            center_of_mass[nrow] += loc_bs[nrow, ncolumn]
    
    center_of_mass /= npart
   
    loc_bs_com = loc_bs.copy()
        
    for nrow in range(ndim):
        loc_bs_com[nrow, :] -= center_of_mass[nrow]
                
    return loc_bs_com, center_of_mass

def pbc(loc_bs_com):
    
    '''
    This function allows the calculation of locations with
    respect to the periodic cell boundaries (acronym PBC). That is 
	one reason why transforming the locations in box size was important.
    
    If the locations are beyond half of the box, then, its
    mirror image is taken up (considering mirror image convention, MIC).
    
    1.0 is supposed to be the length of the box. Locations are 
    fractions since they are represented in terms of the box
    length.
    '''
    loc_bs_com_pbc = loc_bs_com.copy() # create a copy to avoid changing both
    
    loc_bs_com_pbc[loc_bs_com_pbc > 0.5] -= 1.0
    loc_bs_com_pbc[loc_bs_com_pbc < -0.5] += 1.0

    return loc_bs_com_pbc

def within_box(x, y):
    '''
    This function would allow us to check if the distance is zero or
    positive, and then, 1.0 would be subtracted from the distance,
    following minimum image convention. 
	
 	If the distance is negative, 1.0 would be added.
    '''
    if y >= 0:
        return x
    else:
        return -x
    
def compute_force(ndim, npart, loc_bs_com_pbc, boxsize, eps, sig, r_cutoff,
                  u_at_cutoff, du_at_cutoff):
    '''
    Acceleration, potential energy, and the virial are considered
    to be zero. Every time, a calculation is done for a new step,
    all of these parameters are turned to zero.

 	All the parameters, that have been returned, have been calculated
  	relying on standard equations of Lennard-Jones potential. 
    
    When the absolute distance is greater than half of the box, 
    then, minimum image convention (MIC) is applied in order to
    calculate the minimum distance between a pair of particles.
    
    When calculating r2, boxsize has been multiplied to make sure
    distance is no longer calculated in box length unit, rather in
    real units. 
    '''
    
    acc = np.zeros((ndim, npart))
    pot_en = np.zeros(npart)
    virial = 0.0

    for iatom in range(npart - 1):
        for jatom in range(iatom + 1, npart):
            r_box_unit = [loc_bs_com_pbc[0, iatom] - loc_bs_com_pbc[0, jatom],
                          loc_bs_com_pbc[1, iatom] - loc_bs_com_pbc[1, jatom],
                          loc_bs_com_pbc[2, iatom] - loc_bs_com_pbc[2, jatom]]

            if abs(r_box_unit[0]) > 0.5:
                r_box_unit[0] -= within_box(1.0, r_box_unit[0])
            if abs(r_box_unit[1]) > 0.5:
                r_box_unit[1] -= within_box(1.0, r_box_unit[1])
            if abs(r_box_unit[0]) > 0.5:
                r_box_unit[2] -= within_box(1.0, r_box_unit[2])

            r2 = np.dot(r_box_unit, r_box_unit) * (boxsize ** 2)

            if r2 < r_cutoff ** 2:
                r1 = r2 ** 0.5
                ri2 = 1.0 / r2
                ri6 = ri2 ** 3
                ri12 = ri2 ** 6
                sig6 = sig ** 6
                sig12 = sig ** 12

                u = 4.0 * eps * ((sig12 * ri12) - (sig6 * ri6)) - u_at_cutoff \
                    - (r1 * du_at_cutoff)
                du = 24.0 * eps * ri2 * ((2.0 * sig12 * ri12) - (sig6 * ri6)) \
                    + (du_at_cutoff * (ri2 ** 0.5))

                pot_en[jatom] += u
                virial -= du * r2

                for naxis in range(ndim):
                    acc[naxis, iatom] += du * r_box_unit[naxis]
                    acc[naxis, jatom] -= du * r_box_unit[naxis]

    return acc, virial, pot_en

def update_loc(loc_bs_com_pbc, vel, acc, dt):
    '''
	The location of the particles needs to be updated at each step. And 
 	the same has been done with the Velocity Verlet algorithm. 
	'''
    loc_bs_upd = loc_bs_com_pbc.copy()
    
    loc_bs_upd += (dt * vel) + (0.5 * acc * dt ** 2) 
    
    return loc_bs_upd

def calculate_temp(ndim, npart, vel, boxsize, reference_temp):
    '''
	The temperature has been calculated with the implementation of
 	average kinetic energy calculation. 

  	Scaling factor is calculated later on to make sure that the
   	temperature doesn't soar exponentially at the initial point.

 	At the beginning, starting with randomly distributed particles,
  	some may overlap and get immediately repelled at the initial time step. 
   	The repulsion force then heats the system up, sometimes, up to 
	the infinite temperature. To keep that in check, velocity rescaling
 	is done. 
 	'''
    kin_en = np.zeros((npart))
    v2 = np.zeros((npart))
    
    for iatom in range(npart):
        for naxis in range(ndim):
            v2[iatom] += (vel[naxis, iatom] * boxsize) ** 2
        kin_en[iatom] = 0.5 * v2[iatom]
        
    average_kin_en = sum(kin_en) / npart
    instant_temp = 2.0 * average_kin_en / ndim
    scaling_factor = (reference_temp / instant_temp) ** 0.5
    
    return average_kin_en, instant_temp, scaling_factor 

def update_vel(vel, scaling_factor, dt, acc):
    '''
    The velocity of the atoms is updated with the Velocity
    Verlet algorithm.
    '''
    vel_upd = vel.copy()
    
    vel_upd = (scaling_factor * vel_upd) + (0.5 * dt * acc) * 2
    
    return vel_upd

def update_vel_unscaled(vel, dt, acc):
    '''
 	This is a duplicated function, created to make sure that running NVE simulation
  	becomes easier. With X% steps running as NVT ensemble and (100-X)% running as 
   	NVE ensemble, changing the function from update_vel_scaled to update_vel_unscaled
	makes the reading of the code easier. 
    '''
    vel_upd = vel.copy()
    
    vel_upd = vel_upd + (0.5 * dt * acc) * 2
    
    return vel_upd

#............................................................................
#............................................................................
#............................................................................

# NVT 20

import numpy as np
import random as rd


#............................................................................
#............................................................................
#............................................................................

loc_bs, vel = initiate_loc_vel_acc(ndim, npart, boxsize)

loc_bs, center_of_mass = ref_to_com(ndim, npart, loc_bs)

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
print("\n>>>>>> Statistical Averages <<<<<<\n")
print("Temperature = %f K " %(total_temp / effective_scans))
print("Kinetic energy = %f eps" %(total_kin_en / effective_scans))
print("Potential energy = %f eps" %(total_pot_en / effective_scans))
print("Total energy = %f eps" %((total_kin_en + total_pot_en) / effective_scans))
print("Pressure = %f eps/(sig^3) \n" %(total_pressure / effective_scans))
