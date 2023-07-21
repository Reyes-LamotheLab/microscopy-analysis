import numpy as np
import matplotlib.pyplot as plt
from cellpose import utils
import statistics as s
import os
import pandas as pd
import math
import copy
from matplotlib.ticker import PercentFormatter
import seaborn as sns
import tkinter as tk
from tkinter import filedialog, messagebox


def point_in_poly(x, y, poly):
    # Check if a point (x,y) is inside a polygon defined by a set of vertices

    n = len(poly)
    inside = False

    p1x, p1y = poly[0]

    for i in range(1, n + 1):
        p2x, p2y = poly[i % n]
        if y > min(p1y, p2y):
            if y <= max(p1y, p2y):
                if x <= max(p1x, p2x):
                    if p1y != p2y:
                        xints = (y - p1y) * (p2x - p1x) / (p2y - p1y) + p1x
                    if p1x == p2x or x <= xints:
                        inside = not inside
        p1x, p1y = p2x, p2y

    return inside


def cell_lenght(outline):
    # ---initialize the min, max values of the coordinates of the outline of one cell
    maxx = float('-inf')
    maxy = float('-inf')
    minx = float('inf')
    miny = float('inf')

# --This for loop finds the max and min value for x and y coordinates: maxx= max (x) pixel coord, minx= min (x) pixel coord, maxy= max (y) pixel coord, miny= "" ""
    for i in range(len(outline)):

        x = outline[i][0]
        y = outline[i][1]

        if x > maxx:
            maxx = x
        if y > maxy:
            maxy = y
        if y < miny:
            miny = y
        if x < minx:
            minx = x

    # Find the lenght using the Pythagoras formula
    lenght = math.sqrt((maxx-minx)**2+(maxy-miny)**2)
    # Find the approximate center of the cell
    center_coordinate = [minx+(maxx-minx)/2, miny+(maxy-miny)/2]

    return lenght, center_coordinate


def table_creation(data, header, outlines, size, min):

    # -----------Initialization of variables -----------
    counts = np.empty(0, dtype=int)  # The number of particles per cell
    # the final table with 10 columns and 1000 rows (1000 is arbitrary to not cause error)
    table = np.zeros((1000, 10))
    count = 0
    one = []  # one ,two , ... are lists that will contain the counts for each different classes of cell lenght chosen
    two = []
    three = []
    four = []
    rel_dist = []
    z = 0  # this var
    coloc = {}
    cellnum = 1

    if header[2] == "X_(px)":
        col_x = 2
        col_y = 3
        col_int = 12
    else:
        col_x = 1
        col_y = 2
        col_int = 11

    # -----------create tables
    for i in range(len(outlines)):
        poly = outlines[i]
        lenght, center = cell_lenght(poly)

        if lenght > min:
            count = 0

            for j in range(len(data)):

                # CHANGE THIS SO THE COLUMN IS ALIGNED WITH X_(px)
                x = data[j][col_x]
                # CHANGE THIS SO THE COLUMN IS ALIGNED WITH Y_(px)
                y = data[j][col_y]
                coord = [x, y]

                if point_in_poly(x, y, poly):  # if particle inside the cell segmentation

                    distance = math.sqrt((y-center[1])**2+(x-center[0])**2)
                    rel_dist.append(distance/(lenght/2))

                    count += 1

                    table[z][0] = i+1
                    table[z][1] = count
                    table[z][3] = lenght
                    table[z][4] = x
                    table[z][5] = y
                    # CHANGE THIS SO THE COLUMN IS ALIGNED IntegratedInt
                    table[z][6] = data[j][col_int]
                    table[z][7] = center[0]
                    table[z][8] = center[1]
                    table[z][9] = distance

                    if cellnum not in coloc:

                        coloc[cellnum] = []
                        coloc[cellnum].append(coord)

                    else:

                        coloc[cellnum].append(coord)

                    z += 1
            counts = np.append(counts, count)
   

            if (count == 0):
                table[z][0] = i+1
                table[z][3] = lenght
                table[z][7] = center[0]
                table[z][8] = center[1]  # the rest of the columns will be 0
                coloc[cellnum] = []

                z += 1

            else:
                y = z
                for w in range(0, count):
                    table[y-1][2] = count
                    y -= 1

            cellnum += 1

            # -----------change variable "function" to True or False depending if you want histograms for each diff cell size------

            function = True

            if (function):

                nanometer = size*lenght

                if (nanometer >= 2000 and nanometer <= 2800):
                    one.append(count)
                elif (nanometer > 2800 and nanometer <= 3600):
                    two.append(count)
                elif (nanometer > 3600 and nanometer <= 4400):
                    three.append(count)
                elif (nanometer > 4400):
                    four.append(count)
                else:
                    continue
        else:
            continue


    return table, counts, average, one, two, three, four, coloc, rel_dist, cellnum


def counts_histogram(total_1, total_2):
    global number_of_cells

    global results_folder_path
    if total_2.size == 0:
        fig = plt.figure()
        bins1 = np.arange(0, int(np.max(total_1)), 1)
        ax = fig.add_subplot(111)
        ax.hist(total_1, weights=np.ones(len(total_1)) / len(total_1),  alpha=0.6, label="channel 1,  n=" +
                str(number_of_cells), ec="black", rwidth=0.5,  bins=bins1-0.5, facecolor='yellow')
        # ax.hist(total_2, weights=np.ones(len(total_2)) / len(total_2), alpha=0.75, ec="black", label='"channel' 2,  n="+str(sum(total_2)), rwidth = 0.5,  bins=bins1-0.25, facecolor='red')
        plt.xticks(bins1)

        ax.set_title("Particles per cell for each channel")
        ax.set_xlabel('Particles per cell')
        ax.set_ylabel('Percentage for each channel')

        plt.gca().yaxis.set_major_formatter(PercentFormatter(1))
        plt.legend(loc="upper right")
        plt.savefig(str(results_folder_path)+'/Counts.png')

    else:

        fig = plt.figure()
        bins1 = np.arange(0, int(max(np.max(total_1), np.max(total_2))), 1)
        ax = fig.add_subplot(111)
        ax.hist(total_1, weights=np.ones(len(total_1)) / len(total_1),  alpha=0.6, label="channel 1,  n=" +
                str(sum(total_1)), ec="black", rwidth=0.5,  bins=bins1-0.5, facecolor='yellow')
        ax.hist(total_2, weights=np.ones(len(total_2)) / len(total_2), alpha=0.75, ec="black",
                label="channel 2,  n="+str(sum(total_2)), rwidth=0.5,  bins=bins1-0.25, facecolor='red')
        plt.xticks(bins1)

        ax.set_title("Particles per cell for each channel")
        ax.set_xlabel('Particles per cell')
        ax.set_ylabel('Percentage for each channel')

        plt.gca().yaxis.set_major_formatter(PercentFormatter(1))
        plt.legend(loc="upper right")
        plt.savefig(str(results_folder_path)+'/Counts.png')
        plt.show()


def size_histogram(one, two, three, four, chan, total):

    global results_folder_path
    fig = plt.figure()
    bins1 = np.arange(6)
    ax = fig.add_subplot(111)
    ax.hist(one, weights=np.ones(len(one)) / len(total),  alpha=0.6, label='2 ≤size< 2.8' +
            ', n=' + str(len(one)), ec="black", rwidth=0.2,  bins=bins1-0.8, facecolor='yellow')
    ax.hist(two, weights=np.ones(len(two)) / len(total), alpha=0.75, label='2.8 ≤size< 3.6' +
            ', n=' + str(len(two)), ec="black", rwidth=0.2,  bins=bins1-0.60, facecolor='red')
    ax.hist(three, weights=np.ones(len(three)) / len(total),  alpha=0.6, label='3.6 ≤size< 4.4' +
            ', n=' + str(len(three)), ec="black", rwidth=0.2,  bins=bins1-0.40, facecolor='blue')
    ax.hist(four, weights=np.ones(len(four)) / len(total), alpha=0.75, label='4.4 ≤size' +
            ', n='+str(len(four)), ec="black", rwidth=0.2,  bins=bins1-0.20, facecolor='green')

    ax.set_title(
        "Particle per cell for different cell sizes for channel:"+str(chan))
    ax.set_xlabel('Particles per cell')
    ax.set_ylabel('Proportion of particles')
    plt.legend(fontsize='small')
    plt.xticks(bins1)
    plt.gca().yaxis.set_major_formatter(PercentFormatter(1))
    plt.legend(loc="upper right")
    plt.savefig(str(results_folder_path)+"/"+str(chan)+'Cell_size.png')
    plt.show()


def coloc_histogram(data, tot_channel, title, pix, channel):

    global results_folder_path
    total = sum(tot_channel)
    fig, ax = plt.subplots()
    print("There is "+str(total)+" particles for the channel: " + str(channel))
    x_axis = []
    for size in data.keys():

        x_axis.append(size*pix)

    bars = ax.bar(x_axis, [j/total for j in data.values()],
                  label='n=' + str(total)+' channel: ' + str(channel),  ec="black", width=10)

    # set the y-axis limit
    ax.set_ylim(0, 1)
    # set the x-axis label
    ax.set_xlabel('Max distance in nanometers')
    # set the y-axis label
    ax.set_ylabel('Colocalization rate')
    # set the title
    ax.set_title(str(title))

    for bar in bars:
        height = bar.get_height()
        ax.text(bar.get_x() + bar.get_width() / 2, height,
                f'{height:.2f}', ha='center', va='bottom')

    plt.xticks(x_axis)

    plt.legend(loc="upper right")
    plt.savefig(str(results_folder_path)+'/Colocalization.png')
    plt.show()


def distance_histogram(data, total, channel):

    global results_folder_path
    fig, ax = plt.subplots()
    bin_width = 1/40
    bin_edges = np.arange(0, 0.51, bin_width)
    hist, edges = np.histogram(data/2, bins=bin_edges)

# plot the data as a bar chart
    plt.bar(edges[:-1], hist, width=bin_width, align='edge',
            label='n=' + str(sum(total)), ec="black")

    # set the y-axis limit
    ax.set_xlim(0, 0.5)
    # set the x-axis label
    ax.set_xlabel('Absolute relative distance to the center')
    # set the y-axis label
    ax.set_ylabel('Quantity of particles')
    # set the title
    ax.set_title(
        "Histogram of the relative distance of particles\n from the cell center for Channel: "+str(channel))

    plt.legend(loc="upper right")
    # show the plot
    plt.savefig(str(results_folder_path)+"/"+'chan'+str(channel)+'Distance.png')
    plt.show()


def coloc_size(dict_size, tot1, tot2, tot3, tot4, pixel, channel): 
    global results_folder_path
    colc = []
    for val in dict_size.values():
        sum_by_size = np.sum(val, axis=0)
        colc.append(sum_by_size)

    colc2 = []
    tots = [sum(tot1), sum(tot2), sum(tot3), sum(tot4)]
    for a in colc:
        quotient = []
        for i in range(0, 4):
            quotient.append(a[i]/tots[i])
        colc2.append(quotient)

# -------------HISTOGRAM FOR THE CELL SIZE COLOCS-------------------------
    x_labels = [x*pixel for x in [0.5, 1, 1.5, 2, 2.5, 3]]
    data = [list(x) for x in zip(*colc2)]
    sizes = ["2 to 2.8 microns", "2.8 to 3.6 microns",
            "3.6 to 4.4 microns", "4.4 and above microns"]
    lili=-1
    for sublist in data:
      
        lili+=1
        # Multiply weights by 100
        heights = [weight * 100 for weight in sublist]

        plt.bar(x_labels, heights, width=10, alpha=0.5,  ec="black")

        # Display the exact frequency above each bar
        for x, y in zip(x_labels, heights):
            plt.text(x, y, str(round(float(y), 2)), ha='center', va='bottom')

        plt.ylim(0, 100)
        plt.xticks(x_labels)
        plt.xlabel('Max distance in nanometers')
        plt.ylabel('Colocalization rate')
        plt.title(sizes[data.index(sublist)])
        plt.savefig(str(results_folder_path)+'/'+'coloc_per_size'+str(sizes[lili])+'chan'+str(channel)+'.png')
        plt.show()


def integrated_intensity(folder_integrated):
    global results_folder_path
 
    # This iterate over files in directory_res and put it in a list
    files = []
    for filename in os.listdir(folder_integrated):
        f = os.path.join(folder_integrated, filename)
        if not filename.endswith(".DS_Store"):
            files.append(f)

    files.sort()
    IntDen = []
    for m in range(len(files)):

        file_int = np.genfromtxt(files[m], delimiter=',', skip_header=1)
        for n in range(len(file_int)):
            IntDen.append(file_int[n][4])

    sns.swarmplot(data=IntDen, size=3)
    plt.xlabel('Data')
    plt.ylabel('Values')
    plt.title('Strip Plot')
    plt.savefig(str(results_folder_path)+"/"+'IntDen_violinplot.png')


def allowed_dist(cell1, cell2, maxcol_dist):
    dist = math.sqrt((cell1[0]-cell2[0])**2+(cell1[1]-cell2[1])**2)
    if (dist > maxcol_dist):
        return False, dist
    else:
        return True, dist


def colocalization(Input_chan, other_chan, pixel_dist):
    result = 0
    boo = False
    candidate = float('inf')

    for key1 in Input_chan.keys():  # iterates through all the cells of the image

        # if the cell chosen chanel has no particle go to next cell
        if len(Input_chan[key1]) == 0:
            continue

        # if the cell other chanel has no particle return result=0
        elif len(other_chan[key1]) == 0:
            continue

        else:

            for j in range(0, len(Input_chan[key1])):
                smallest_dist = float('inf')

                for y in range(0, len(other_chan[key1])):
                    allowed, dist = allowed_dist(
                        Input_chan[key1][j], other_chan[key1][y], pixel_dist)

                    if allowed and dist <= smallest_dist:
                        smallest_dist = dist
                        candidate = y
                        boo = True

                    else:
                        continue
                if boo:
                    del other_chan[key1][candidate]
                    result += 1
                boo = False

    return result


def size_colocalization(Input_chan, other_chan, pixel_dist, table, size):
    df = pd.DataFrame(table)
    result = [0, 0, 0, 0]
    boo = False
    candidate = float('inf')
    x = True
    j = 0
    for key1 in Input_chan.keys():  # iterates through all the cells of the image

        # if the cell chosen chanel has no particle go to next cell
        if len(Input_chan[key1]) == 0:
            continue

        if len(other_chan[key1]) == 0:
            # if the cell other chanel has no particle return result=0
            continue

        else:

            for j in range(0, len(Input_chan[key1])):
                smallest_dist = float('inf')

                for y in range(0, len(other_chan[key1])):
                    allowed, dist = allowed_dist(
                        Input_chan[key1][j], other_chan[key1][y], pixel_dist)

                    if allowed and dist <= smallest_dist:
                        smallest_dist = dist
                        candidate = y
                        boo = True

                    else:
                        continue
                if boo:
                
                    if candidate < len(other_chan[key1]):
                        del other_chan[key1][candidate]
                    else:
                        print(len(other_chan[key1]))
                        print("Invalid candidate index:", candidate)
                    # del other_chan[key1][candidate]

                    value = df.set_index('cell_id').loc[key1, 'cell_lenght']

                    if value.size >= 2:
                        value = value.iloc[0]

                    lenght = value*size

                    if (lenght >= 2000 and lenght <= 2800):
                        result[0] = result[0]+1
                    elif (lenght > 2800 and lenght <= 3600):
                        result[1] = result[1]+1
                    elif (lenght > 3600 and lenght <= 4400):
                        result[2] = result[2]+1
                    elif (lenght > 4400):
                        result[3] = result[3]+1
                    else:
                        continue

                boo = False

    return result



def main1(directory_seg, directory_res, Coloc_bysize , folder_integrated, parameter, channel, selected_folder_path):

    global results_folder_path
    global number_of_cells

    results_folder_path=selected_folder_path
  


    files=[]
    for filename in os.listdir(directory_res):
        f = os.path.join(directory_res, filename)
        if not filename.endswith(".DS_Store"):
            files.append(f)
    files.sort()

    if folder_integrated != "":
        integrated_intensity(folder_integrated)


# SOME INPUT TO ASK THE USER
    pixel = parameter['pixel_size']
    too_small = parameter['cell_length']
    names = ['cell_id', 'part_id', 'total_count', 'cell_lenght',
             'Part_X_(px)', 'Part_Y_(px)', 'IntegatedIntensity', 'Cellcenter_X(px)', 'Cellcenter_Y(px)', 'Distance']

# ------ #total particles for each channel --------------------------------------------------------------------------------------------------------------------------------------------
    total_1 = np.empty(0)
    total_2 = np.empty(0)

# ------ Variables to make the histogram for particles per cell size------------------------------------------------------------------------------------------------------------------
    tot_1 = np.empty(0) # for channel 1 divided in cell size
    tot_2 = np.empty(0)
    tot_3 = np.empty(0)
    tot_4 = np.empty(0)
    tot_1_2 = np.empty(0)  # for 'channel' 2 divided in cell
    tot_2_2 = np.empty(0)
    tot_3_2 = np.empty(0)
    tot_4_2 = np.empty(0)

# ------Variable to make the histogram of colocalization rate per distance size in pixels----------------------------------------------------------------------------------------------
    # Here you can change the number of distances you want by deleting or adding the desired pixel distance,
    new_dict = {0.5: 0, 1.0: 0, 1.5: 0, 2: 0, 2.5: 0, 3: 0}
    dict_size = {0.5: [], 1.0: [], 1.5: [], 2: [], 2.5: [], 3: []}
    dist1 = np.empty(0)
    dist2 = np.empty(0)

    i = 0
    number_of_cells=0
    for filename in sorted(os.listdir(directory_seg)):

        if not filename.endswith(".DS_Store"):
            print(filename)

            seg_file = os.path.join(directory_seg, filename)
            seg_masks = np.load(seg_file, allow_pickle=True).item()
            #masks = seg_masks['masks']
            outlines = utils.outlines_list(seg_masks['masks'])
            number_of_cells += len(outlines)

            # input channel 1 CHANGE COMMA OR TAB IF NEEDED
            fiji_1 = np.genfromtxt(files[2*i], delimiter='\t', skip_header=1)
            header = np.genfromtxt(files[2 * i], delimiter='\t', dtype=str, max_rows=1)
            # input channel 2 CHANGE COMMA OR TAB IF NEEDED
            fiji_2 = np.genfromtxt(files[2*i+1], delimiter='\t', skip_header=1)
            i += 1

            table_1, counts_1, average_1, one, two, three, four, coloc, distance, cellnum1 = table_creation(
                fiji_1,header, outlines, pixel, too_small)

            table_2, counts_2, average_2, one2, two2, three2, four2, coloc2, distance2, cellnum2 = table_creation(
                fiji_2,header, outlines, pixel, too_small)  # channel 2


# ------------------# of particle per cell histogram for both channel -------------

            total_1 = np.concatenate((total_1, counts_1), axis=None)
            total_2 = np.concatenate((total_2, counts_2), axis=None)

# --------------------PARTICLE PER CELL SIZE HISTOGRAM-----------------------------
            tot_1 = np.concatenate((tot_1, one), axis=None)
            tot_2 = np.concatenate((tot_2, two), axis=None)
            tot_3 = np.concatenate((tot_3, three), axis=None)
            tot_4 = np.concatenate((tot_4, four), axis=None)

            tot_1_2 = np.concatenate((tot_1_2, one2), axis=None)
            tot_2_2 = np.concatenate((tot_2_2, two2), axis=None)
            tot_3_2 = np.concatenate((tot_3_2, three2), axis=None)
            tot_4_2 = np.concatenate((tot_4_2, four2), axis=None)

            dist1 = np.concatenate((dist1, distance),  axis=None)
            dist2 = np.concatenate((dist2, distance2),  axis=None)

# -----------------Saving results of tables created here in CSV file ---------------

            df1 = pd.DataFrame(table_1, columns=names)
            file_path = os.path.join(
                selected_folder_path, 'First_channel'+str(i)+'.csv')
            df1.to_csv(file_path, index=True, header=True, sep=',')
            df2 = pd.DataFrame(table_2, columns=names)
            file_path = os.path.join(
                selected_folder_path, 'Second_channel'+str(i)+'.csv')
            df2.to_csv(file_path, index=True, header=True, sep=',')

    
# ----------------- COLOCALIZATION-------------------------

            for key in new_dict.keys():

                inco = copy.deepcopy(coloc)
                othco = copy.deepcopy(coloc2)
                rate = colocalization(inco, othco, key)
                new_dict[key] += rate

            for key in new_dict.keys():

                inco = copy.deepcopy(coloc)
                othco = copy.deepcopy(coloc2)
                res = size_colocalization(inco, othco, key, df1, size=pixel)
                res = np.array(res)
                dict_size[key].append(res)

 # --------------Plot the histrogram of counts per cell
    counts_histogram(total_1, total_2)

# --------------Plot histogram for each size of cell
    size_histogram(tot_1, tot_2, tot_3, tot_4, 1, total_1)
    size_histogram(tot_1_2, tot_2_2, tot_3_2, tot_4_2, 2, total_2)

# -------------Plot histogram for colocalization
    if Coloc_bysize:

        if channel == 1:
     
            coloc_size(dict_size, tot_1, tot_2, tot_3, tot_4, pixel, channel)
        else:
            coloc_size(dict_size, tot_1_2, tot_2_2, tot_3_2, tot_4_2, pixel, channel)

    if channel == 1:
        coloc_histogram(new_dict, total_1,
                        "Channel 1 on channel 2 ", pixel, channel)
    else:
        coloc_histogram(new_dict, total_2,
                        "Channel 2 on channel 1", pixel, channel)
        
# distance histogram from the center for both channels 
    distance_histogram(dist1, total_1, channel=1)
    distance_histogram(dist2, total_2, channel=2)



def main2(directory_seg, directory_res, folder_integrated, parameters, selected_folder_path):
    global results_folder_path
    results_folder_path=selected_folder_path
    global number_of_cells

    pixel = parameters['pixel_size']
    too_small = parameters['cell_length']
    names = ['cell_id', 'part_id', 'total_count', 'cell_lenght',
             'Part_X_(px)', 'Part_Y_(px)', 'IntegatedIntensity', 'Cellcenter_X(px)', 'Cellcenter_Y(px)', 'Distance']

    files=[]
    for filename in os.listdir(directory_res):
        f = os.path.join(directory_res, filename)
        if not filename.endswith(".DS_Store"):
            files.append(f)
    files.sort()

    if folder_integrated != "":
        integrated_intensity(folder_integrated)

    #------- #total particles for each channel --------------------------------------------------------------------------------------------------------------------------------------------
    total_1=np.empty(0)
    total_2=np.empty(0)
    #------ Variables to make the histogram for particles per cell size------------------------------------------------------------------------------------------------------------------
    tot_1=np.empty(0) #for channel 1 divided in cell size
    tot_2=np.empty(0)
    tot_3=np.empty(0) 
    tot_4=np.empty(0)
    #------ Variable to make the histogram of colocalization rate per distance size in pixels----------------------------------------------------------------------------------------------

    dist1=np.empty(0)

    i=0
    number_of_cells=0
    for filename in sorted(os.listdir(directory_seg)):
    
        if not filename.endswith(".DS_Store"):
            print(filename)
        
            seg_file = os.path.join(directory_seg, filename)
            seg_masks=np.load(seg_file, allow_pickle=True).item()
            masks=seg_masks['masks']
            outlines= utils.outlines_list(seg_masks['masks'])
            number_of_cells+=len(outlines)
            with open(files[i], 'r', newline='') as file:
                file_content = file.read().splitlines()

            with open(files[i], 'r', newline='') as file1:
                file_content1 = file1.read()
            
                
            fiji_1= np.genfromtxt(file_content[1:], delimiter='\t')   # input channel 1
            if fiji_1.ndim==1: 
                fiji_1 = fiji_1.reshape(1, -1)
                print("changed it to an arrray of dim", fiji_1.shape)


            header = np.genfromtxt([file_content[0]], delimiter='\t', dtype=str)

            i+=1

            table_1, counts_1, average_1, one, two, three, four, coloc, distance, cellnum1=table_creation(fiji_1,header, outlines,pixel, too_small)
        #------------------num of particle per cell histogram dataset-------------

            total_1 =np.concatenate((total_1, counts_1), axis=None) 
        #--------------------particle per cell size dataset-----------------------------

            tot_1=np.concatenate((tot_1, one), axis=None) 
            tot_2=np.concatenate((tot_2, two), axis=None) 
            tot_3=np.concatenate((tot_3, three), axis=None) 
            tot_4=np.concatenate((tot_4, four), axis=None) 

        #--------------distance from the center dataset--------------
            dist1=np.concatenate((dist1,distance),  axis=None)

        #-----------------Saving results of tables created here in CSV file ---------------
            df1 = pd.DataFrame(table_1, columns=names)
            file_path = os.path.join(selected_folder_path, 'table'+str(i)+'.csv')
            df1.to_csv(file_path, index=True, header=True, sep=',')
        
    #--------------Plot the histrogram of counts per cell--------------
    counts_histogram (total_1, total_2)
    #--------------Plot histogram for each size of cell--------------
    size_histogram (tot_1, tot_2, tot_3, tot_4, 1, total_1)
    #------------Plot histogram for relative distance from the cell center--------------
    distance_histogram(dist1, total_1, 1)



class App:
    def __init__(self):

        #OPEN START WINDOW
        self.root = tk.Tk()
        self.root.title("Analysis")

        # Initialize attributes
        self.selected_folder_path = None
        self.directory_seg=None
        self.directory_res=None
        self.folder_integrated=None
        self.checkbox_values = {}
        self.parameters = {}

        # Welcome message
        self.message_label = tk.Label(self.root, text="Welcome! Before starting make sure your segmentations are valid and that your particle detection is good! Make sure that you have a folder of _seg.npy files only and one for the output of particle_detections.ijm script!")
        self.message_label.pack()

        # Folder selection options
        self.folder_selection = tk.StringVar()
        self.folder_selection.set("new")  # Default to create new folder

        self.new_folder_radiobutton = tk.Radiobutton(self.root, text="Create New Folder",
                                                     variable=self.folder_selection, value="new")
        self.new_folder_radiobutton.pack()

        self.existing_folder_radiobutton = tk.Radiobutton(self.root, text="Choose Existing Folder",
                                                          variable=self.folder_selection, value="existing")
        self.existing_folder_radiobutton.pack()

        self.proceed_button = tk.Button(self.root, text="Proceed", command=self.proceed)
        self.proceed_button.pack()
        
    def proceed(self):
        selected_option = self.folder_selection.get()

        if selected_option == "new":
            self.create_new_folder()
        elif selected_option == "existing":
            self.choose_existing_folder()
        
            
    def choose_existing_folder(self):
        self.root.destroy()
        if self.show_messages("Choose the folder where you want results saved"):
            self.selected_folder_path = filedialog.askdirectory()
            self.selected_folder_path = self.selected_folder_path.replace('\\', '/')
            self.select_analysis()


    def create_new_folder(self):

        self.root.destroy()
        # Open folder selection window
        self.folder_window = tk.Tk()
        self.folder_window.title("Folder Selection")
        self.folder_window_label=tk.Label(self.folder_window, text="Enter the name of the new folder which will contain all the results (Tables and Figures):" )
        self.folder_window_label.pack(pady=5)
        self.entry_newfolder=tk.Entry(self.folder_window)        
        self.entry_newfolder.pack(pady=5)

        self.status_newf = tk.Label(self.folder_window, text="")
        self.status_newf.pack(pady=5)
        self.btn_create=tk.Button(self.folder_window, text="Create Folder", command= self.create_folder)
        self.btn_create.pack(pady=10)


    def create_folder(self):

        folder_name = self.entry_newfolder.get()

        if folder_name:
            # Implement your logic here to create the new folder
            target_path= filedialog.askdirectory(title="Select a directory to place the folder")
            if target_path:
                new_path=os.path.join(target_path, folder_name)
                try:
                    os.mkdir(new_path)
                    self.selected_folder_path=new_path
                    self.status_newf.config(text=f"Folder '{folder_name}' created successfully.")
                    self.select_analysis()
                    self.folder_window.destroy()
                    

                except OSError as e:
                    self.status_newf.config(text=f"Error: {e}")
            else:

                self.status_newf.config(text="No directory selected.")
        else:
            self.status_newf.config(text="Please enter a folder name.")
        


    def update_question2_var(self):
        self.question2_var.set(True)
    def update_question1_var(self):
        self.question1_var.set(True)

    def select_analysis(self):
  
        # New window for questions
        self.questions_window = tk.Tk()
        self.questions_window.title("Select the Analysis you need")

        # Questions with Yes/No options
        self.question1_label = tk.Label(self.questions_window, text="Do you want to analyse the Integrated Intensity of the cells assuming you have results from the script_intensity.ijm code?")
        self.question1_label.pack()
        self.question1_var = tk.BooleanVar()
        self.question1_var.set(False)
        self.question1_checkbox_yes = tk.Radiobutton(self.questions_window, text="Yes", variable=self.question1_var, value=True, command=self.update_question1_var)
        self.question1_checkbox_yes.pack()
        self.question1_checkbox_no= tk.Radiobutton(self.questions_window, text="No", variable=self.question1_var, value=False)
        self.question1_checkbox_no.pack()

        self.question2_label = tk.Label(self.questions_window, text="Do you want to have the colocalization of your particles?")
        self.question2_label.pack()
        self.question2_var = tk.BooleanVar()
        self.question2_var.set(False)
        self.question2_checkbox_yes= tk.Radiobutton(self.questions_window, text="Yes", variable=self.question2_var, value=True, command=self.update_question2_var)
        self.question2_checkbox_yes.pack()
        self.question2_checkbox_no= tk.Radiobutton(self.questions_window, text="No", variable=self.question2_var, value=False)
        self.question2_checkbox_no.pack()

        # Submit button
        self.submit_button = tk.Button(self.questions_window, text="Submit", command=self.save_checkbox_values)
        self.submit_button.pack()


    def save_checkbox_values(self):
        # Save checkbox values

        self.checkbox_values['Integrated'] = self.question1_var.get()
       
        self.checkbox_values["Coloc"] = self.question2_var.get()
        

        if self.checkbox_values["Coloc"]:

        # New Window for questions about colocalization
            self.questions_window.destroy()
            self.colc_window=tk.Tk()
            self.colc_window.title("Colocalization customization")

        # Questions with Yes/No options
            self.q1_label = tk.Label(self.colc_window, text="Additionally, would you like the colocalization rate based on cell size category?")
            self.q1_label.pack()
            self.q1_var = tk.StringVar()
            
            self.q1_yes= tk.Radiobutton(self.colc_window, text="Yes", variable=self.q1_var, value="Yes", command=lambda: self.coloc_param)
            self.q1_yes.pack()
            self.q1_no= tk.Radiobutton(self.colc_window, text="No", variable=self.q1_var, value="No",  command=lambda: self.coloc_param)
            self.q1_no.pack()
            
            self.q2_label = tk.Label(self.colc_window, text="Which Channel do you want to colocalization be based on?")
            self.q2_label.pack()
            self.q2_var = tk.IntVar()

            self.q2_chan1 = tk.Radiobutton(self.colc_window, text="Channel 1", variable=self.q2_var,  command=lambda: self.coloc_param, value=1)
            self.q2_chan1.pack()
            self.q2_chan2 = tk.Radiobutton(self.colc_window, text="Channel 2", variable=self.q2_var, command=lambda: self.coloc_param, value=2)
            self.q2_chan2.pack()

            # Submit button
            self.submit_button = tk.Button(self.colc_window, text="Submit", command=self.coloc_param)
            self.submit_button.pack()
     
        else:
            self.questions_window.destroy()
            self.parameters_entry()
            
       

    def coloc_param(self):
        
        self.checkbox_values['coloc_size']=self.q1_var.get()

        self.checkbox_values['channel']=self.q2_var.get()
  
        self.parameters_entry()
        self.colc_window.destroy()

    def parameters_entry(self):

        # New window for parameter input
        self.parameter_window = tk.Tk()
        self.parameter_window.title("Parameter Input")

        # Integer parameter inputs
        self.param1_label = tk.Label(self.parameter_window, text="Pixel size in nanometer (integer):")
        self.param1_label.pack()
        self.param1_entry = tk.Entry(self.parameter_window)
        self.param1_entry.pack()

        self.param2_label = tk.Label(self.parameter_window, text="Minimum accepted cell length for the cell to be used in the analysis (integer):")
        self.param2_label.pack()
        self.param2_entry = tk.Entry(self.parameter_window)
        self.param2_entry.pack()

        # Submit button
        self.submit_button = tk.Button(self.parameter_window, text="Submit", command=self.save_parameters)
        self.submit_button.pack()


        self.param1_entry.insert(tk.END, "65")
        self.param2_entry.insert(tk.END, "0")

    def save_parameters(self):
        # Save parameter values
    
        self.parameters['pixel_size'] = int(self.param1_entry.get())
        self.parameters['cell_length'] = int(self.param2_entry.get())
        self.call_other_functions()
        
    

    def show_messages(self, files):

        response=messagebox.showinfo("Folder/directory to input", "Please Enter the directory for the "+files+" files")
        return response=="ok"

    def call_other_functions(self):
        self.parameter_window.destroy()
        # Example: Displaying the saved information
        print("Selected Directory:", self.selected_folder_path)
        print("Checkbox Values:", self.checkbox_values)
        print("Parameters:", self.parameters)
        
        if self.show_messages("Segmentations (_seg.npy)"):
            self.directory_seg = filedialog.askdirectory()
            self.directory_seg = self.directory_seg.replace('\\', '/')
        
        if self.show_messages("Particle detection results (.csv or .txt)"):
            self.directory_res = filedialog.askdirectory()
            self.directory_res = self.directory_res.replace('\\', '/')

        if self.checkbox_values['Integrated']:
            self.show_messages("Integrated Intensity results (.csv or .txt)") 
            self.folder_integrated = filedialog.askdirectory()
            self.folder_integrated = self.folder_integrated.replace('\\', '/')
        else:
            self.folder_integrated =""

    def run(self):
        self.root.mainloop()
        if self.checkbox_values["Coloc"]:
            if self.checkbox_values['coloc_size']=="Yes": 
                messagebox.showinfo("Analysis", "You have two color channels. Output: Cells Length Distribution, Particle's Relative Distance from the Center, Counts of Particles per Cell, Colocalization, Colocalization per cell length size ")
                main1(self.directory_seg, self.directory_res, True, self.folder_integrated, self.parameters, self.checkbox_values['channel'], self.selected_folder_path)
            else:
                messagebox.showinfo("Analysis", "You have two color channels. Output: Cells Length Distribution, Particle's Relative Distance from the Center, Counts of Particles per Cell, Colocalization")
                main1(self.directory_seg, self.directory_res, False, self.folder_integrated, self.parameters,self.checkbox_values['channel'], self.selected_folder_path)
        else:
            messagebox.showinfo("Analysis", "You have one color channel, Output: Cells Length Distribution, Particle's Relative Distance from the Center and Counts of Particles per Cell")
            main2(self.directory_seg, self.directory_res, self.folder_integrated, self.parameters, self.selected_folder_path)


if __name__ == '__main__':
    app = App()
    app.run()
