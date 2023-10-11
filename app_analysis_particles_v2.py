import glob
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


bin1=2000  #less than 
bin2=2500 
bin3=3000
bin4=3500
bin5=4000 #more than

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
    five=[]
    six=[]
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

                if (nanometer < bin1 ):
                    one.append(count)
                elif (nanometer > bin1 and nanometer <= bin2):
                    two.append(count)
                elif (nanometer > bin2 and nanometer <= bin3):
                    three.append(count)
                elif (nanometer > bin3 and nanometer <= bin4):
                    four.append(count)
                elif(nanometer > bin4 and nanometer <= bin5):
                    five.append(count)
                elif(nanometer > bin5):
                    six.append(count)
                
                else:
                    continue
        else:
            continue


    return table, counts, one, two, three, four, five, six, coloc, rel_dist, cellnum


def counts_histogram(total_1, total_2):
    global number_of_cells
    global results_folder_path

    if total_2.size == 0:
        fig = plt.figure()
        bins1 = np.arange(0, int(np.max(total_1)), 1)
        ax = fig.add_subplot(111)
        ax.hist(total_1, weights=np.ones(len(total_1)) / len(total_1),  alpha=0.6, label="channel 1,  number of particles=" +
                str(sum(total_1)), ec="black", rwidth=0.5,  bins=bins1-0.5, facecolor='yellow')
        # ax.hist(total_2, weights=np.ones(len(total_2)) / len(total_2), alpha=0.75, ec="black", label='"channel' 2,  n="+str(len(total_2)), rwidth = 0.5,  bins=bins1-0.25, facecolor='red')
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


def size_histogram(one, two, three, four, five, six, chan, total):

    global results_folder_path
    fig = plt.figure()
    bins1 = np.arange(6)
    ax = fig.add_subplot(111)
    ax.hist(one, weights=np.ones(len(one)) / len(total),  alpha=0.6, label=f'{bin1} <size' +
            ', n=' + str(sum(one)), ec="black", rwidth=0.2,  bins=bins1-0.8, facecolor='yellow')
    ax.hist(two, weights=np.ones(len(two)) / len(total), alpha=0.75, label=f'{bin1} <size≤ {bin2}' +
            ', n=' + str(sum(two)), ec="black", rwidth=0.2,  bins=bins1-0.60, facecolor='red')
    ax.hist(three, weights=np.ones(len(three)) / len(total),  alpha=0.6, label=f'{bin2} <size≤ {bin3}' +
            ', n=' + str(sum(three)), ec="black", rwidth=0.2,  bins=bins1-0.40, facecolor='blue')
    ax.hist(four, weights=np.ones(len(four)) / len(total), alpha=0.75, label=f'{bin3} <size≤ {bin4} ' +
            ', n='+str(sum(four)), ec="black", rwidth=0.2,  bins=bins1-0.20, facecolor='green')

    ax.set_title(
        "Particle per cell for different cell sizes for channel:"+str(chan))
    ax.set_xlabel('Particles per cell')
    ax.set_ylabel('Proportion of particles')
    plt.legend(fontsize=7)
    plt.xticks(bins1)
    plt.gca().yaxis.set_major_formatter(PercentFormatter(1))
    plt.legend(loc="upper right")
    plt.savefig(str(results_folder_path)+"/"+str(chan)+'Cell_size.png')


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


 # Function to filter bins with less than 5 counts

def filter_bins(data, range_data):
    counts, bin_edges = np.histogram(data, bins='auto', range=range_data)
    filtered_bins = [bin_edges[i] for i, count in enumerate(counts) if count >= 5]
    return filtered_bins

def integrated_int_histogram():
    global results_folder_path
    csv_files = glob.glob(os.path.join(results_folder_path, '*.csv'))

    # Initialize DataFrames for each channel
    first_channel_df = pd.DataFrame()
    second_channel_df = pd.DataFrame()
    color_channel1='cyan'
    color_channel2='yellow'


    # Iterate through the CSV files
    for csv_file in csv_files:
        channel_name = None
        if os.path.basename(csv_file).startswith('First_channel'):
            channel_name = 'First_channel'
        elif os.path.basename(csv_file).startswith('Second_channel'):
            channel_name = 'Second_channel'
    
        if channel_name:
            # Read the CSV file
            df = pd.read_csv(csv_file)
            
            # Filter out rows with 'IntegratedIntensity' >= 200
            df = df[df['IntegratedIntensity'] >= 200]
            
            # Append the data to the corresponding DataFrame
            if channel_name == 'First_channel':
                first_channel_df = pd.concat([first_channel_df, df['IntegratedIntensity']], axis=0)
            else:
                second_channel_df = pd.concat([second_channel_df, df['IntegratedIntensity']], axis=0)

    channel1_counts=first_channel_df.shape[0]
    channel2_counts=second_channel_df.shape[0]

    mean_first_channel = first_channel_df[0].mean(axis=0)
    std_dev_first_channel = first_channel_df[0].std(axis=0)
    range_first_channel = (0, mean_first_channel + 3 * std_dev_first_channel)

    mean_second_channel = first_channel_df[0].mean(axis=0)
    std_dev_second_channel = first_channel_df[0].std(axis=0)
    range_second_channel = (0, mean_second_channel + 3 * std_dev_second_channel)

    filtered_bins_first_channel = filter_bins(first_channel_df, range_first_channel)
    filtered_bins_second_channel = filter_bins(second_channel_df, range_second_channel)

    plt.figure(figsize=(10, 5))
    plt.hist(first_channel_df.dropna(), bins=filtered_bins_first_channel, alpha=0.5,ec="black", color=color_channel1, label=f'First_channel (Total Counts: {channel1_counts})')
    plt.xlabel('Integrated Intensity')
    plt.ylabel('Counts')
    plt.title('Histogram of Integrated Intensity - First_channel')
    plt.legend()
    output_file1 = os.path.join(results_folder_path, f'_First_channel_particleInt.png')
    plt.savefig(output_file1)

    plt.figure(figsize=(10, 5))
    plt.hist(second_channel_df.dropna(),bins=filtered_bins_second_channel, alpha=0.5,ec="black", color=color_channel2, label=f'Second_channel (Total Counts: {channel2_counts})')
    plt.xlabel('Integrated Intensity')
    plt.ylabel('Counts')
    plt.legend()
    plt.title('Histogram of Integrated Intensity - Second_channel')
    output_file = os.path.join(results_folder_path, f'_Second_channel_particleInt.png')

    plt.savefig(output_file)
    plt.show()


def coloc_size_histo(dict_size, tot1, tot2, tot3, tot4,tot5,tot6, pixel, channel): 
    global results_folder_path
    colc = []   #list of the number of colocalizations per cell size
    for i, val in enumerate(dict_size.values()):
        sum_by_size = np.sum(val, axis=0) # sum all columns, since each column corresponds to a size
        colc.append(sum_by_size)

    colc2 = []  #list of np.arrays() each np  array is for a certain coloc distance
    tots = np.array([sum(tot1), sum(tot2), sum(tot3), sum(tot4),sum(tot5), sum(tot6)])
    print("TOTALS ARRAY", tots)
    for a in colc:
        a=np.array(a)
        quotient = np.array(a/tots)

       # for i in range(0, 4):
        #    quotient.append(a[i]/tots[i])
        colc2.append(quotient)

# -------------HISTOGRAM FOR THE CELL SIZE COLOCS for different colocalization distances-------------------------

    # x_labels = [x*pixel for x in [0.5, 1, 1.5, 2, 2.5, 3]]
    # data = [list(x) for x in zip(*colc2)]
    # sizes_fig = ["2 to 2.8 microns", "2.8 to 3.6 microns",
    #         "3.6 to 4.4 microns", "4.4 and above microns"]
    # lili=-1
    # for sublist in data:
    #     plt.figure()
    #     lili+=1
    #     # Multiply weights by 100
    #     heights = [weight * 100 for weight in sublist]

    #     plt.bar(x_labels, heights, width=10, alpha=0.5,  ec="black")

    #     # Display the exact frequency above each bar
    #     for x, y in zip(x_labels, heights):
    #         plt.text(x, y, str(round(float(y), 2)), ha='center', va='bottom')

    #     plt.ylim(0, 100)
    #     plt.xticks(x_labels)
    #     plt.xlabel('Max distance in nanometers')
    #     plt.ylabel('Colocalization rate')
    #     plt.title(sizes[data.index(sublist)])
    #     plt.savefig(str(results_folder_path)+'/'+'coloc_per_size'+str(sizes[lili])+'chan'+str(channel)+'.png')

# -------------HISTOGRAM FOR THE CELL SIZE COLOCS for 130 nanometers only -------------------------
    #x_labels = [x*pixel for x in [0.5, 1, 1.5, 2, 2.5, 3]]

    data = colc2[3]  # Select the list in the third index

    print("DATA COLOC FOR 130", data)
    x_labels = [f"{bin1}<size", f"{bin1}< size <{bin2}", f"{bin2}< size <{bin3}", f"{bin3}< size <{bin4}", f"{bin4}< size <{bin5}", f"{bin5}< size"]
    x_values = np.arange(len(x_labels))  # Numerical x-axis values

    heights = [weight * 100 for weight in data]
    plt.figure()

    # Create the bar chart with custom x-axis labels
    plt.bar(x_values, heights, alpha=0.5, edgecolor="black", label=f'bin1 n={sum(tot1)},\nbin2 n={sum(tot2)},\nbin3 n={sum(tot3)},\nbin4 n={sum(tot4)},\nbin5 n={sum(tot5)},\nbin6 n={sum(tot6)},\n')
    plt.xticks(x_values, x_labels, fontsize=5)

    # Display the exact frequency above each bar
    for x, y in zip(x_values, heights):
        plt.text(x, y, str(round(float(y), 2)), ha='center', va='bottom')
    plt.legend(fontsize=8)

    plt.ylim(0, 100)
    plt.xlabel('Cell_length')
    plt.ylabel('Colocalization rate')
    plt.title('Histogram of the colocalization rate for \n a distance of maximum 130 nanometers (2 pixels)')
    plt.savefig(str(results_folder_path) + '/' + 'coloc_per_size for 130 nanometers' + str(channel) + '.png')


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

        if (len(Input_chan[key1]) != 0 and len(other_chan[key1]) != 0):

            for j,can1 in enumerate(Input_chan[key1]):

                smallest_dist = float('inf')
                if len(other_chan[key1]) != 0:

                    for y, can2 in enumerate(other_chan[key1]):

                        allowed, dist = allowed_dist(
                            can1, can2, pixel_dist)

                        if allowed and dist <= smallest_dist:
                            smallest_dist = dist
                            boo = True

                        else:
                            continue
                    if boo:
                        del other_chan[key1][y]
                        result += 1
                    boo = False
                else:
                    continue

    return result



def size_colocalization(Input_chan, other_chan, pixel_dist, table, size):
    df = pd.DataFrame(table)
    result = [0, 0, 0, 0, 0, 0]
    boo = False
    candidate = float('inf')
    x = True
    j = 0
    for key1 in Input_chan.keys():  # iterates through all the cells of the image


        if (len(Input_chan[key1]) != 0 and len(other_chan[key1]) != 0):

            for j,can1 in enumerate(Input_chan[key1]):
                smallest_dist = float('inf')
                if len(other_chan[key1]) != 0:

                    for y, can2 in enumerate(other_chan[key1]):
                        allowed, dist = allowed_dist(
                            can1, can2, pixel_dist)

                        if allowed and dist <= smallest_dist:
                            smallest_dist = dist
                            boo = True

                        else:
                            continue
                    if boo:
                    
                        if len(other_chan[key1])==0:
                            print("Invalid candidate index:", candidate)
                 
                        else:
                            del other_chan[key1][y]
           
                        value = df.set_index('cell_id').loc[key1, 'cell_lenght']

                        if value.size >= 2:
                            value = value.iloc[0]

                        lenght = value*size

                        if (lenght < bin1):
                            result[0] = result[0]+1  
                        elif (lenght > bin1 and lenght <= bin2):
                            result[1] = result[1]+1
                        elif (lenght > bin2 and lenght <= bin3):
                            result[2] = result[2]+1
                        elif (lenght > bin3 and lenght <=bin4):
                            result[3] = result[3]+1
                        elif (lenght > bin4 and lenght <=bin5):
                            result[4] = result[4]+1
                        elif (lenght > bin5 ):
                            result[5] = result[5]+1
                            
                        else:
                            continue

                    boo = False
                else:
                    continue

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
    too_small = parameter['cell_length']  # this is in pixel length 
    names = ['cell_id', 'part_id', 'total_count', 'cell_lenght',
             'Part_X_(px)', 'Part_Y_(px)', 'IntegratedIntensity', 'Cellcenter_X(px)', 'Cellcenter_Y(px)', 'Distance']

# ------ #total particles for each channel --------------------------------------------------------------------------------------------------------------------------------------------
    total_1 = np.empty(0)
    total_2 = np.empty(0)

# ------ Variables to make the histogram for particles per cell size------------------------------------------------------------------------------------------------------------------
    tot_1 = np.empty(0) # for channel 1 divided in cell size
    tot_2 = np.empty(0)
    tot_3 = np.empty(0)
    tot_4 = np.empty(0)
    tot_5 = np.empty(0)
    tot_6 = np.empty(0)
    tot_1_2 = np.empty(0)  # for 'channel' 2 divided in cell
    tot_2_2 = np.empty(0)
    tot_3_2 = np.empty(0)
    tot_4_2 = np.empty(0)
    tot_5_2 = np.empty(0)
    tot_6_2 = np.empty(0)

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
            #fiji_1 = np.genfromtxt(files[2*i], delimiter='\t', skip_header=1)
            #header = np.genfromtxt(files[2 * i], delimiter='\t', dtype=str, max_rows=1)
            # input channel 2 CHANGE COMMA OR TAB IF NEEDED
            #fiji_2 = np.genfromtxt(files[2*i+1], delimiter='\t', skip_header=1)

            with open(files[2*i], 'r', newline='') as file:
                file_content = file.read().splitlines()
            
                
            fiji_1= np.genfromtxt(file_content[1:], delimiter='\t')   # input channel 1
            if fiji_1.ndim==1: 
                fiji_1 = fiji_1.reshape(1, -1)
                #print("changed it to an arrray of dim", fiji_1.shape)

            header = np.genfromtxt([file_content[0]], delimiter='\t', dtype=str)

            with open(files[2*i+1], 'r', newline='') as file2:
                file_content2 = file2.read().splitlines()
            
                
            fiji_2= np.genfromtxt(file_content2[1:], delimiter='\t')   # input channel 1
            if fiji_2.ndim==1: 
                fiji_2 = fiji_2.reshape(1, -1)
                #print("changed it to an arrray of dim", fiji_2.shape)


            i += 1

            table_1, counts_1, one, two, three, four,five, six, coloc, distance, cellnum1 = table_creation(
                fiji_1,header, outlines, pixel, too_small)

            table_2, counts_2, one2, two2, three2, four2,five2, six2, coloc2, distance2, cellnum2 = table_creation(
                fiji_2,header, outlines, pixel, too_small)  # channel 2


# ------------------# of particle per cell histogram for both channel -------------

            total_1 = np.concatenate((total_1, counts_1), axis=None)
            total_2 = np.concatenate((total_2, counts_2), axis=None)

# --------------------PARTICLE PER CELL SIZE HISTOGRAM-----------------------------
            tot_1 = np.concatenate((tot_1, one), axis=None)
            tot_2 = np.concatenate((tot_2, two), axis=None)
            tot_3 = np.concatenate((tot_3, three), axis=None)
            tot_4 = np.concatenate((tot_4, four), axis=None)
            tot_5 = np.concatenate((tot_5, five), axis=None)
            tot_6 = np.concatenate((tot_6, six), axis=None)

            tot_1_2 = np.concatenate((tot_1_2, one2), axis=None)
            tot_2_2 = np.concatenate((tot_2_2, two2), axis=None)
            tot_3_2 = np.concatenate((tot_3_2, three2), axis=None)
            tot_4_2 = np.concatenate((tot_4_2, four2), axis=None)
            tot_5_2= np.concatenate((tot_5, five2), axis=None)
            tot_6_2= np.concatenate((tot_6, six2), axis=None)


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
    size_histogram(tot_1, tot_2, tot_3, tot_4,tot_5,tot_6, 1, total_1)
    size_histogram(tot_1_2, tot_2_2, tot_3_2, tot_4_2, tot_5_2, tot_6_2, 2, total_2)

# -------------Plot histogram for colocalization
    if Coloc_bysize:

        if channel == 1:
     
            coloc_size_histo(dict_size, tot_1, tot_2, tot_3, tot_4,tot_5,tot_6, pixel, channel)
        else:
            coloc_size_histo(dict_size, tot_1_2, tot_2_2, tot_3_2, tot_4_2, tot_5_2, tot_6_2, pixel, channel)

    if channel == 1:
        coloc_histogram(new_dict, total_1,
                        "Channel 1 on channel 2 ", pixel, channel)
    else:
        coloc_histogram(new_dict, total_2,
                        "Channel 2 on channel 1", pixel, channel)
        
# distance histogram from the center for both channels 
    distance_histogram(dist1, total_1, channel=1)
    distance_histogram(dist2, total_2, channel=2)
    integrated_int_histogram()
"""
    res_file=[]
    for file in os.listdir(results_folder_path):
        f = os.path.join(results_folder_path, file)
        if not file.endswith(".DS_Store"):
            res_file.append(f)
    res_file.sort()

    for i in res_file:
        res_table= np.genfromtxt(i[1:], delimiter='\t')   # input channel 1
        if res_table.ndim==1: 
            res_table = res_table.reshape(1, -1)
        header = np.genfromtxt([i[0]], delimiter='\t', dtype=str)
        print(header, len(header), len(res_table[1]))"""

def main2(directory_seg, directory_res, folder_integrated, parameters, selected_folder_path):
    global results_folder_path
    results_folder_path=selected_folder_path
    global number_of_cells

    pixel = parameters['pixel_size']
    too_small = parameters['cell_length'] 
    names = ['cell_id', 'part_id', 'total_count', 'cell_lenght',
             'Part_X_(px)', 'Part_Y_(px)', 'IntegratedIntensity', 'Cellcenter_X(px)', 'Cellcenter_Y(px)', 'Distance']

    files=[]
    for filename in os.listdir(directory_res):
        f = os.path.join(directory_res, filename)
        if not filename.endswith(".DS_Store"):
            files.append(f)
    files.sort()

    if folder_integrated != "":
        integrated_intensity(folder_integrated)

    #------- #total particles for each channel (total 2 = empty in this case--------------------------------------------------------------------------------------------------------------
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

            table_1, counts_1, one, two, three, four, coloc, distance, cellnum1, integrated_int1=table_creation(fiji_1,header, outlines,pixel, too_small)
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
    integrated_int_histogram()
"""
    res_file=[]
    for file in os.listdir(results_folder_path):
        f = os.path.join(results_folder_path, file)
        if not (file.endswith(".DS_Store") and file.endswith(".png")):
            res_file.append(f)
    res_file.sort()

    for i in res_file:

        res_table= np.genfromtxt(i[1:], delimiter='\t')   # input channel 1
        if res_table.ndim==1: 
            res_table = res_table.reshape(1, -1)
        header = np.genfromtxt([i[0]], delimiter='\t', dtype=str)
        print(header, len(header), len(res_table[1]))

"""


class App:
    def __init__(self):

        #OPEN START WINDOW
        self.root = tk.Tk()
        self.root.title("Analysis")

        # Initialize attributes
        self.selected_folder_path = None
        self.directory_seg=None  #directory that contains the segmentations
        self.directory_res=None  #where we store the results
        self.folder_integrated=None #folder that contains integrated intensity values 
        self.channel=0
        self.parameters={}
        self.col_var=""
        self.int_var=""
        self.size_var=""
        self.chan_var=""

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
        


    def select_analysis(self):
  
        # New window for questions
        self.questions_window = tk.Tk()
        self.questions_window.title("Select the Analysis you need")

        # Questions with Yes/No options
        self.int_label = tk.Label(self.questions_window, text='Do you want to analyse the Integrated Intensity of the cells assuming you have results from the script_intensity.ijm code? Write "y" or "n" ' )
        self.int_label.pack(pady=2)
        self.int_entry=tk.Entry(self.questions_window )
        self.int_entry.pack(pady=2)



        self.col_label = tk.Label(self.questions_window, text='Do you want to have the colocalization of your particles? Write "y" or "n" ')
        self.col_label.pack(pady=2)
        self.col_entry=tk.Entry(self.questions_window )
        self.col_entry.pack(pady=2)



        # Submit button
        self.submit_button = tk.Button(self.questions_window, text="Submit", command=self.coloc_analysis_customizations)
        self.submit_button.pack(pady=5)
        self.status_analysis = tk.Label(self.questions_window, text="")
        self.status_analysis.pack(pady=5)


    def coloc_analysis_customizations(self):
        
        self.int_var=self.int_entry.get()
        self.col_var=self.col_entry.get()
        print(self.col_var, self.int_var)

        if (self.col_var=="y" or self.col_var=="n") and (self.int_var=="y" or self.int_var=="n"):
            print("here")
            if self.col_var=="y":
                print("here2")
                self.questions_window.destroy()
            # New Window for questions about colocalization
                self.colc_window=tk.Tk()
                self.colc_window.title("Colocalization customization")

            # Questions with Yes/No options
                #ask for colocalization by size
                self.size_label = tk.Label(self.colc_window, text="Additionally, would you like the colocalization rate based on cell size category? Write 'y' or 'n'")
                self.size_label.pack(pady=2)
                self.size_entry=tk.Entry(self.colc_window)
                self.size_entry.pack(pady=2) 


                #ask for the channel
                self.chan_label = tk.Label(self.colc_window, text="Which Channel do you want to colocalization be based on? Enter 1 or 2")
                self.chan_label.pack(pady=2)
                self.chan_entry=tk.Entry(self.colc_window)
                self.chan_entry.pack(pady=2)
           
                # Submit button
                self.submit_button = tk.Button(self.colc_window, text="Submit", command=self.parameters_entry)
                self.submit_button.pack(pady=2)
            else:
                print('here')
                self.parameters_entry()

        else:
            self.status_analysis.config(text='Please enter only the letters "y" and "n" corresponding to yes or no.')
      

    def parameters_entry(self):

        # New window for parameter input
        if self.col_var=="y":
            self.size_var=self.size_entry.get()
            self.chan_var=self.chan_entry.get()
            self.colc_window.destroy()

        self.parameter_window = tk.Tk()
        self.parameter_window.title("Parameter Input")

        # Integer parameter inputs
        self.param1_label = tk.Label(self.parameter_window, text="Pixel size of the microscope in nanometer (integer):")
        self.param1_label.pack()
        self.param1_entry = tk.Entry(self.parameter_window)
        self.param1_entry.pack()

        self.param2_label = tk.Label(self.parameter_window, text="Minimum accepted cell length (in pixel unit) for the cell to be used in the analysis. Enter the length as an integer only.")
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
       # print("Selected Directory:", self.selected_folder_path)
       # print("Checkbox Values:", self.checkbox_values)
      #  print("Parameters:", self.parameters)
        
        if self.show_messages("Segmentations (_seg.npy)"):
            self.directory_seg = filedialog.askdirectory()
            self.directory_seg = self.directory_seg.replace('\\', '/')
        
        if self.show_messages("Particle detection results (.csv or .txt)"):
            self.directory_res = filedialog.askdirectory()
            self.directory_res = self.directory_res.replace('\\', '/')

        if self.int_var=="y":
            self.show_messages("Integrated Intensity results (.csv or .txt)") 
            self.folder_integrated = filedialog.askdirectory()
            self.folder_integrated = self.folder_integrated.replace('\\', '/')
        else:
            self.folder_integrated =""

    def run(self):
        self.root.mainloop()
        if self.col_var=="y":
            if self.size_var=="y": 
                print("ENTERED IN THE FIRST COLOC SIZE CONDITION")
                if self.chan_var=="1":
                    self.channel=1
                else:
                    self.channel=2
                messagebox.showinfo("Analysis", "You have two color channels. Output: Cells Length Distribution, Particle's Relative Distance from the Center, Counts of Particles per Cell, Colocalization, Colocalization per cell length size ")
                main1(self.directory_seg, self.directory_res, True, self.folder_integrated, self.parameters, self.channel, self.selected_folder_path)
            else:
                messagebox.showinfo("Analysis", "You have two color channels. Output: Cells Length Distribution, Particle's Relative Distance from the Center, Counts of Particles per Cell, Colocalization")
                main1(self.directory_seg, self.directory_res, False, self.folder_integrated, self.parameters,self.channel, self.selected_folder_path)
        else:
            messagebox.showinfo("Analysis", "You have one color channel, Output: Cells Length Distribution, Particle's Relative Distance from the Center and Counts of Particles per Cell")
            main2(self.directory_seg, self.directory_res, self.folder_integrated, self.parameters, self.selected_folder_path)


if __name__ == '__main__':
    app = App()
    app.run()
