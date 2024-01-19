# Created on 13-Oct-2023

def read_file(file_path):
    
    # open the file without the need to close the file
    with open(file_path,'r') as file_content:
        line_count = sum(1 for num_lines in file_content) #this line allows checking of all the lines in the file
        file_content.seek(0) #brings the cursor back to the start after complete iteration
        all_lines = [next(file_content) for _ in range(line_count)] #reads and stores all the lines from the file

    start_counting = 2 #start the counting of the lines, first two lines are redundant
    stop_counting = line_count - 1 #the last line is also redundant 

    '''
    In this particular file, only the coordinates and velocities are required. Nothing else is necessary.
    In this code, all the lines have been split into their component elements and have been stored as a list of elements.

    The list is comprised of all the elements, like Residue Number, Residue Name, Atom Number, X, Y, Z, Vx, Vy, Vz.
    Since there is a bit of anomaly after the 10000th line, only the last six columns are being called to store. To do that,
    an empty numpy array has been created with 6 columns and a requisite number of rows (line_count - 3). The first two and the last
    lines are of no use. 

    .split() allows splitting all the elements from all the lines.
    '''
    ow_only = []

    for lines1 in all_lines[start_counting:stop_counting:5]:
        ow_only += [lines1.split()]

    row_ow = len(ow_only)
    store_data = np.zeros((row_ow, 3))

    for i in range(len(ow_only)):
        store_data[i][0] = ow_only[i][-3] #x coordinate
        store_data[i][1] = ow_only[i][-2] #y coordinate
        store_data[i][2] = ow_only[i][-1] #z coordinate

    ow_xyz = store_data.copy() #store a copy to avoid making changes in all the data
    
    return ow_xyz, all_lines

def check_box_dimension(all_lines):
    
    '''
    This function just allows tracing the box dimensions in .gro file that would be used in the final 
    q6total calculation. 
    '''

    box_len = np.zeros(3)
    last_line = all_lines[-1].split()

    for i in range(len(box_len)):
        box_len[i] = last_line[i]

    return box_len

def q6atomwise(total_atoms, ow_xyz, box_len, cn_total, q6total, cn_ice, f_q6atomwise, cutoff):

    '''
    For iatom, all jatom distance would be calculated at first. Then, those jatom distance values would be
    taken into account whose o_range values are smaller than 0.35, the cut-off box_lennce. 

    For those jatoms, who make the cut, let's call them jcatom. Now, for every jatom inside the if loop,
    there would be one extra coordination that would be counted. Each jatom would require a calculation
    of all the y values. Once the y value calculation for one jcatom is done, then, other jatom would be
    iterated. If this jatom also makes the cut, then, this jcatom would add another coordination, and
    y values would be updated with the new calculations.

    Now, these y values contain values for 2 jatom and the coordination number is 2. After the calculation
    of all jatom, we would have total y values. If 10 jatom are there, then, the coordination number would be 10
    and all y values of all those jcatom would be summed up.

    Once calculation for that iatom is done with all other jatom, all those y values would be divided by 
    the number of all jatom. These are <y> values for that precise iatom. Subsequently, all the <y> values 
    would be used to calculate their conjugates.

    Using those conjugate values, the final q6 value of that iatom would be calculated. 

    Now, if the q6 value satisfies a certain condition, say, it lies between 0.0 and 1.0, that q6 value would
    be added to the q6total and that iatom would be added to the number of atoms that has
    q6 value fulfilling this condition, 0.0 < q6 < 1.0.

    Finally, q6total would mean summing all those iatom whose q6 had fulfilled the criterion and their
    total count, say coord_total. By dividing q6total with the coord_total, the q6average value
    can be calculated. 
    
    Note: All the calculations are being performed on a nanometer scale. Remember, all other analyses using MDAnalysis
    would use Angstrom.
    '''

    for iatom in range(total_atoms - 1):

        y1, y2, y3, y4, y5, y6, y7 = [0] * 7
        y8, y9, y10, y11, y12, y13 = [0] * 6
        coord_num = 0

        for jatom in range(iatom + 1, total_atoms):

            o_distance = ow_xyz[iatom] - ow_xyz[jatom]
            o_distance_box = o_distance - box_len * np.round(o_distance / box_len)
            o_range = np.sqrt(np.dot(o_distance_box, o_distance_box))

            if o_range < cutoff:

                coord_num += 1
                x, y, z = o_distance_box[0], o_distance_box[1], o_distance_box[2]

                r = np.sqrt(x**2 + y**2 + z**2)
                y0 = (x - y * imaginary) / r
                y00 = (x + y * imaginary) / r
                y1 += ((np.sqrt(3003.0 / pi)) / 64.0) * (y0**6)
                y2 += (((np.sqrt(1001.0 / pi)) * 3.0) / 32.0) * (y0**5) * (z / r)
                y3 += (((np.sqrt(91.0 * 0.5 / pi)) * 3.0) / 32.0) * (y0**4) * (((11.0 * (z**2)) - (r**2)) / (r**2))
                y4 += ((np.sqrt(1365.0 / pi)) / 32.0) * (y0**3) * z * (((11.0 * (z**2)) - (3.0 * (r**2)))/(r**3))
                y5 += ((np.sqrt(1365.0 / pi)) / 64.0) * (y0**2) * (((33.0 * (z**4)) - (18.0 * (z**2) * (r**2)) + (r**4)) / r**4)
                y6 += ((np.sqrt(273.0 * 0.5 / pi)) / 16.0) * y0 * z * (((33.0 * (z**4)) - (30.0 * (z**2) * (r**2))+5.0*(r**4))/(r**5))
                y7 += ((np.sqrt(13.0 / pi)) / 32.0) * (((231.0 * (z**6) - 315.0 * (z**4) * (r**2) + 105.0 * (z**2) * (r**4) - 5.0 * (r**6))) / (r**6))
                y8 += (-(((np.sqrt(273.0 * 0.5 / pi)) / 16.0) * y00 * z * (((33.0 * (z**4)) - (30.0 * (z**2) * (r**2)) + (5.0 * (r**4)))/(r**5))))
                y9 += ((np.sqrt(1365.0 / pi)) / 64.0) * (y00**2) * ((33.0 * (z**4) - 18.0 * (z**2) * (r**2) + (r**4)) / (r**4))
                y10 += (-(((np.sqrt(1365.0 / pi)) / 32.0) * (y00**3) * z * ((11.0 * ( z**2) - (3.0 * (r**2)))/(r**3))))
                y11 += (((np.sqrt(91.0 * 0.5 / pi)) * 3.0) / 32.0) * (y00**4) * ((11.0 * (z**2)-(r**2))/(r**2))
                y12 += (-((((np.sqrt(1001.0 / pi)) * 3.0) / 32.0) * (y00**5) * (z / r)))
                y13 += ((np.sqrt(3003.0 / pi)) / 64.0) * (y00**6)

        if coord_num != 0:

            #calculation of y-values
            y1 /= coord_num
            y2 /= coord_num
            y3 /= coord_num
            y4 /= coord_num
            y5 /= coord_num
            y6 /= coord_num
            y7 /= coord_num
            y8 /= coord_num
            y9 /= coord_num
            y10 /= coord_num
            y11 /= coord_num
            y12 /= coord_num
            y13 /= coord_num

            # calculation of conjugates
            y1sq = y1 * np.conj(y1)
            y2sq = y2 * np.conj(y2)
            y3sq = y3 * np.conj(y3)
            y4sq = y4 * np.conj(y4)
            y5sq = y5 * np.conj(y5)
            y6sq = y6 * np.conj(y6)
            y7sq = y7 * np.conj(y7)
            y8sq = y8 * np.conj(y8)
            y9sq = y9 * np.conj(y9)
            y10sq = y10 * np.conj(y10)
            y11sq = y11 * np.conj(y11)
            y12sq = y12 * np.conj(y12)
            y13sq = y13 * np.conj(y13)

            #calculation of q6 value
            q6 = np.sqrt(((4.0 * pi) / 13.0) * (y1sq + y2sq + y3sq + y4sq + y5sq + y6sq + y7sq + y8sq + y9sq + y10sq + y11sq + y12sq + y13sq))
            
            #check the value of q6 and make sure it is in the right range
            if 0.0 < q6 < 1.0:
                cn_total += 1
                q6total += q6 
                
            #check the value of q6 and measure the proportion of ice
            if 0.55 < q6 < 1.0:
                cn_ice += 1
                 
            f_q6atomwise.write("%d %d %e %e %d %d\n" % (iatom, coord_num, q6, q6total, cn_total, cn_ice))
            
    return q6total, cn_total, cn_ice
