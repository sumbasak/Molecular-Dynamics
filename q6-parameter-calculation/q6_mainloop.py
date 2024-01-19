# Created on 13-Oct-2023

# main code for q6total calculation for one gro file
ow_xyz, all_lines = read_file(file_path) #coordinates in nm scale
box_len = check_box_dimension(all_lines) #box length in nm scale
q6total, cn_total, cn_ice = q6atomwise(total_atoms, ow_xyz, box_len, 
                                       cn_total, q6total, cn_ice, f_q6atomwise, cutoff) 

# data visualization with data frame
ox_display = pd.DataFrame(ow_xyz)
ox_display
