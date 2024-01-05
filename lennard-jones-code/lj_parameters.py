
'''
Parameters should be monitored and modulated according to the needs. Using
indecent parameters can make the system crash or provide abysmal results. 

To avoid that, check the parameter details. Use this set of parameters as the 
starting point and then, change them accordingly. 
'''
# parameters

ndim = 3
nscans = 1000
npart = 100
dt = 1.0e-4 #ps
reference_temp = 1.0 #K
sig = 1.0
eps = 1.0
r_cutoff = 3.5 * sig #nm

# calculate LJ potential
u_at_cutoff = 4.0 * eps * ((sig / r_cutoff) ** 12 - (sig / r_cutoff) ** 6)

# calculate the derivative of potential energy at cutoff
term1 = (sig ** 12) / (r_cutoff ** 7)
term2 = (sig ** 6) / (r_cutoff ** 7)
du_at_cutoff = 48.0 * eps * (term1 - 0.5 * term2)

boxsize = 10.0 * sig #nm
volume = boxsize ** 3
rho = npart / volume

[total_pressure, total_kin_en, total_pot_en, total_temp] = np.zeros([4])

# Create the files to write the data
f_tpzsc = open("tpzsc.out", "w")
f_kpesc = open("kpesc.out", "w")
f_xyzsc = open("loc_sc.xyz", "w")

# ignore steps to avoid initial temperature spike
ignore_steps = 30

# functions to check 20% NVT and 80% NVE
check_temp20 = np.zeros(nscans)
check_pressure20 = np.zeros(nscans)
check_kin_en20 = np.zeros(nscans)
check_pot_en20 = np.zeros(nscans)
check_total_en20 = np.zeros(nscans)
