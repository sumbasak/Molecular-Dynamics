# Created on 13-Oct-2023

file_path = '/home/sumon/MEGA/mount/zeus/q6/ice1h_202012.gro'
f_q6atomwise = open('/home/sumon/MEGA/mount/zeus/q6/q6atomwise.out', "w")
cutoff = 0.35 #nm, remember to change it accordingly
total_atoms = len(ow_xyz)
row_dist = total_atoms - 1
imaginary = 0.0 + 1.0j
pi = 22.0 / 7.0
q6total, cn_total, cn_ice = 0, 0, 0
