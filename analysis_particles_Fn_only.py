import glob
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from cellpose import utils
import statistics as s
import os
import pandas as pd
import math
import copy          
from matplotlib.ticker import PercentFormatter
import seaborn as sns
from collections import defaultdict


matplotlib.rcParams['pdf.fonttype'] = 42

# BINNING FOR THE SIZE OF THE CELLS, THE USER CAN CHANGE THE VALUES BUT THE NUMBER OF BINS SHOULD BE THE SAME
bin1 ,bin2, bin3, bin4, bin5=[2000, 2500, 3000, 3500, 4000]



color_channel1 = 'cyan'
color_channel2 = 'yellow'

# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
####################################################################################################################################################################################################################
# ----------------------------------------------------------------------# FUNCTIONS TO CREATE THE PLOTS USING MATPLOTLIB -------------------------------------------------------------------------------------------#
####################################################################################################################################################################################################################
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#

def counts_histogram(total_1, total_2): 
    """
    Output:
    This function creates a histogram of the number of particles per cell for each channel and saves it in the results folder.
    Input:
    - total_1: The number of particles per cell for channel 1.
    - total_2: The number of particles per cell for channel 2.
    """

    global number_of_cells
    global results_folder_path

    if total_2.size == 0:
        fig = plt.figure()
        bins1 = np.arange(0, int(np.max(total_1))+1, 1) 
        print(int(np.max(total_1))+1, 'int(np.max(total_1))+1')
        ax = fig.add_subplot(111)
        ax.hist(total_1, weights=np.ones(len(total_1)) / len(total_1),  alpha=0.6, label="channel 1,  # of particles=" +
                str(sum(total_1)), ec="black", rwidth=0.5,  bins=bins1-0.5, facecolor=color_channel1)
        plt.xticks(bins1)

        ax.set_title("Particles per cell for each channel")
        ax.set_xlabel('Particles per cell')
        ax.set_ylabel('Percentage for each channel')

        plt.gca().yaxis.set_major_formatter(PercentFormatter(1))
        plt.legend(loc="upper right")
        plt.savefig(str(results_folder_path)+'/Counts.pdf', format='pdf')

    else:

        fig = plt.figure()
        bins1 = np.arange(0, int(max(np.max(total_1)+1, np.max(total_2)))+1, 1)

        ax = fig.add_subplot(111)
        ax.hist(total_1, weights=np.ones(len(total_1)) / len(total_1),  alpha=0.6, label="channel 1,  # of particles=" +
                str(sum(total_1)), ec="black", rwidth=0.5,  bins=bins1-0.5, facecolor=color_channel1)
        ax.hist(total_2, weights=np.ones(len(total_2)) / len(total_2), alpha=0.75, ec="black",
                label="channel 2,  n="+str(sum(total_2)), rwidth=0.5,  bins=bins1-0.25, facecolor=color_channel2)
        plt.xticks(bins1)

        ax.set_title("Particles per cell for each channel")
        ax.set_xlabel('Particles per cell')
        ax.set_ylabel('Percentage for each channel')

        plt.gca().yaxis.set_major_formatter(PercentFormatter(1))
        plt.legend(loc="upper right")
        plt.savefig(str(results_folder_path)+'/Counts.pdf', format='pdf')


def size_histogram(one, two, three, four, five, six, chan, total):
    """
    Output:
    This function creates a histogram of the number of particles per cell for A channel and saves it in the results folder.
    However it will divide the particles per cell in 4 different bins depending on the size of the cell.
    It ignores the cells that are too small. 
    Input:
    - The number of particles per cell for each cell size (one, two, three, four, five, six)
    - The channel number.
    - The total number of particles for that channel.
    """


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
    plt.savefig(str(results_folder_path)+"/"+str(chan)+'Cell_size.pdf', format='pdf')


def coloc_histogram(data, tot_channel, title, pix, channel):
    """
    Output:
    This function creates a histogram of the colocalization rate per distance for a channel and saves it in the results folder.
    Input:
    - The data of the colocalization rate per maximal distance of colocalization (a dictionnary).
    - The total number of particles for that channel.
    - The title of the plot.
    - The pixel size and the channel number.
    """

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
    plt.savefig(str(results_folder_path)+'/Colocalization.pdf', format='pdf')


def distance_histogram(data, total, channel,axis):
    """
    Output:
    This function creates a histogram of the relative distance of particles from the cell center for a channel 
    and saves it in the results folder.

    Input:
    - The data of the relative distance of particles from the cell center, 
    - The total number of particles for that channel and the channel number.
    """
    print(data)
    global results_folder_path
    fig, ax = plt.subplots()
    bin_width = 1/40
    bin_edges = np.arange(0, 0.51, bin_width)
    hist, edges = np.histogram(data/2, bins=bin_edges)  # divide by 2 because the whole cell size was 2 and now we want to have the relative distance from the center of a cell size normalized to a length of 1. 

# plot the data as a bar chart
    plt.bar(edges[:-1], hist, width=bin_width, align='edge',
            label='# of particles =' + str(sum(total)), ec="black")

    # set the y-axis limit
    ax.set_xlim(0, 0.5)
    # set the x-axis label
    if axis=="center":
        ax.set_xlabel('Lengthwise relative distance from the middle of the center of the cell')
    else:
        ax.set_xlabel('Lengthwise relative distance from the middle of the center of the cell')
    # set the y-axis label
    ax.set_ylabel('Quantity of particles')
    # set the title
    ax.set_title(
        "Histogram of the Lengthwise Distance From the Middle\n of the Center of the Cell for Channel: "+str(channel))

    plt.legend(loc="upper right")
    # show the plot
    plt.savefig(str(results_folder_path)+"/"+'chan'+str(channel)+'Distance.pdf', format='pdf')


def filter_bins(data, range_data):
    """
    Output:
    This function filters the bins of a histogram by removing the bins that have less than 5 counts. 
    Input: 
    -The data of the histogram and the range of the data.
    """

    counts, bin_edges = np.histogram(data, bins='auto', range=range_data)
    filtered_bins = [bin_edges[i] for i, count in enumerate(counts) if count >= 5]
    return filtered_bins


def integrated_int_histogram(filename, color):
    """
    Output:
    This function creates a histogram of the integrated intensity of the particles for a channel and saves it in the results folder.
    Input:
    - The filename: only the suffix of the filename which are the CSV files that contain the data of the particles.
    - The results of this script from the function 'table_creation' and the color of the channel.
    """

    global results_folder_path
    csv_files = glob.glob(os.path.join(results_folder_path, '*.csv'))
    # Initialize DataFrames for each channel
    channel_df = pd.DataFrame()
    color_channel=color

    # Iterate through the CSV files
    for csv_file in csv_files:
        file = None
        if os.path.basename(csv_file).startswith(filename):
            file = filename
        if file:
            # Read the CSV file
            df = pd.read_csv(csv_file)
            
            # Filter out rows with 'IntegratedIntensity' >= 200
            df = df[df['IntegratedIntensity'] >= 200]

            channel_df = pd.concat([channel_df, df['IntegratedIntensity']], axis=0)

    channel_counts=channel_df.shape[0]
    mean = channel_df[0].mean(axis=0)
    std_dev = channel_df[0].std(axis=0)
    range = (0, mean + 3 * std_dev)

    filtered_bins = filter_bins(channel_df, range)

    plt.figure(figsize=(10, 5))
    plt.hist(channel_df.dropna(),bins=filtered_bins, alpha=0.5,ec="black", color=color_channel, label=f'Second_channel (Total Counts: {channel_counts})')
    plt.xlabel('Integrated Intensity')
    plt.ylabel('Counts')
    plt.legend()
    plt.title(f'Histogram of Integrated Intensity - {filename}')
    output_file = os.path.join(results_folder_path, f'{filename}_particleInt.pdf')
    plt.savefig(output_file)
    plt.show()


def coloc_size_histo(dict_size, tot1, tot2, tot3, tot4,tot5,tot6, pixel, channel): 
    """
    Output:
    This function creates a histogram of the colocalization rate for a distance of 2 pixels per cell size for a channel and saves it in the results folder.
    Input:
    - dict_size: Dictionnary of the colocalization rate per cell size
    - tot1, tot2, tot3, tot4, tot5, tot6: The total number of particles for each cell size
    - pixel: The pixel size
    - channel: The channel number
    """

    global results_folder_path
    colc = []   #list of the number of colocalizations per cell zize
    for i, val in enumerate(dict_size.values()):
        sum_by_size = np.sum(val, axis=0) # sum all columns, since each column corresponds to a size
        colc.append(sum_by_size)

    colc2 = []  #list of np.arrays() each np  array is for a certain coloc distance
    tots = np.array([sum(tot1), sum(tot2), sum(tot3), sum(tot4),sum(tot5), sum(tot6)])
    print("TOTALS ARRAY", tots)
    for a in colc:
        a=np.array(a)
        quotient = np.array(a/tots)
        colc2.append(quotient)

    #x_labels = [x*pixel for x in [0.5, 1, 1.5, 2, 2.5, 3]]

    colocalization_130nm = colc2[3]  # Select the list in the third index since we are interested in the colocalization rate for a distance of 2 pixels only

    # print("DATA COLOC FOR 130", data)
    x_labels = [f"{bin1}<size", f"{bin1}< size <{bin2}", f"{bin2}< size <{bin3}", f"{bin3}< size <{bin4}", f"{bin4}< size <{bin5}", f"{bin5}< size"]
    x_values = np.arange(len(x_labels))  # Numerical x-axis values

    heights = [weight * 100 for weight in colocalization_130nm]
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
    plt.savefig(str(results_folder_path) + '/' + 'coloc_per_size for 130 nanometers' + str(channel) + '.pdf', format='pdf')




def coloc_integrated_plot(filename):
    """
    Output:
    Histogram of the ratio of colocalized particles/total particles per integrated intensity bin. 
    """
    plt.style.use('ggplot')
    global results_folder_path
    csv_files = glob.glob(os.path.join(results_folder_path, '*.csv'))
    big_df=pd.DataFrame()

    for csv_file in csv_files:
        if os.path.basename(csv_file).startswith(filename):
            df = pd.read_csv(csv_file)
            big_df = pd.concat([big_df, df], ignore_index=True) 

    # for each integrated intensity value that falls in the bin, calculate the ratio of 1 to total number of elements in the bin 
    min_value = big_df['IntegratedIntensity'].min()
    max_value = big_df['IntegratedIntensity'].max()
    ratios = []

    # Ensure that the min and max values are valid numbers and min is not greater than max
    if np.isnan(min_value) or np.isnan(max_value) or min_value > max_value:
        print("Invalid min or max values. Cannot create bin edges.")
    else:
        print("Creating bin edges for the integrated intensity values.")
    # Define the number of quantiles
        num_quantiles = 50  # This will create 10 quantile bins

        # Use np.linspace to create evenly spaced quantiles between 0 and 1
        quantiles = np.linspace(0, 1, num_quantiles + 1)

        # Use np.quantile to find the bin edges based on the quantiles
        bin_edges = np.quantile(big_df['IntegratedIntensity'], quantiles)
            # Remove non-unique edges
        bin_edges = np.unique(bin_edges)

        if len(bin_edges) < 2:
            print("Not enough unique bin edges available after removing duplicates.")
            return

    # Use pd.cut to bin the data
    big_df['bins'] = pd.cut(big_df['IntegratedIntensity'], bins=bin_edges)

    
    ratios=big_df.groupby('bins')['coloc'].agg(lambda x: x.sum() / len(x) if len(x) > 0 else np.nan)
    total = big_df.groupby('bins')['total_count'].sum()
   
    

    # Create the primary axis for the ratios
    fig, ax1 = plt.subplots(figsize=(10, 6))

    # Plot the ratios as a bar plot on the primary axis
    ratios.plot(kind='bar', ax=ax1, color='blue', alpha=0.6)

    # Set the labels for the primary axis
    ax1.set_xlabel('Integrated Intensity Bins')
    
    ax1.set_ylabel('Ratio of Colocalization', color='blue')
    ax1.tick_params(axis='y', labelcolor='blue')

    # Define the bin labels for the x-axis
    bin_labels = [f"{int(left)}-{int(right)}" for left, right in zip(bin_edges[:-1], bin_edges[1:])]
    ax1.set_xticklabels(bin_labels, rotation=45, ha='right', fontsize=7)

    # Create a second y-axis to plot the total count
    ax2 = ax1.twinx()

    # Plot the total count as a line plot on the secondary axis
    total.plot(kind='line', ax=ax2, color='red', marker='o', linewidth=2, markersize=6)

    # Set the labels for the secondary axis
    ax2.set_ylabel('Total Number of Particles', color='red')
    ax2.tick_params(axis='y', labelcolor='red')

    # Set the title for the plot
    ax1.set_title('Ratio of Colocalization and Total Number of Particles vs Integrated Intensity')

    plt.tight_layout()  # Adjust layout to fit all labels


    plt.savefig(str(results_folder_path) + '/' + 'coloc_per_integratedIntensity.pdf', format='pdf')


def coloc_integrated_number(filename):
    """
    Output:
    Histogram of the ratio of colocalized particles/total particles per integrated intensity bin. 
    """
    plt.style.use('ggplot')
    global results_folder_path
    csv_files = glob.glob(os.path.join(results_folder_path, '*.csv'))
    big_df=pd.DataFrame()

    for csv_file in csv_files:
        file = None
        print(csv_file)
        if os.path.basename(csv_file).startswith(filename):
            file = filename
        if file:
            df = pd.read_csv(csv_file)
            big_df = pd.concat([big_df, df], axis=0)

    # for each integrated intensity value that falls in the bin, calculate the ratio of 1 to total number of elements in the bin 
    min_value = big_df['IntegratedIntensity'].min()
    max_value = big_df['IntegratedIntensity'].max()
    ratios = []

    # Ensure that the min and max values are valid numbers and min is not greater than max
    if np.isnan(min_value) or np.isnan(max_value) or min_value > max_value:
        print("Invalid min or max values. Cannot create bin edges.")
    else:
    # Define the number of quantiles
        num_quantiles = 50  # This will create 10 quantile bins

        # Use np.linspace to create evenly spaced quantiles between 0 and 1
        quantiles = np.linspace(0, 1, num_quantiles + 1)

        # Use np.quantile to find the bin edges based on the quantiles
        bin_edges = np.quantile(big_df['IntegratedIntensity'], quantiles)
            # Remove non-unique edges
        bin_edges = np.unique(bin_edges)

        if len(bin_edges) < 2:
            print("Not enough unique bin edges available after removing duplicates.")
            return

    # Use pd.cut to bin the data
    big_df['bins'] = pd.cut(big_df['IntegratedIntensity'], bins=bin_edges)

    
    ratios=big_df.groupby('bins')['coloc'].agg(lambda x: x.sum() / len(x) if len(x) > 0 else np.nan)
    total = big_df.groupby('bins')['total_count'].sum()
   
    # Create the primary axis for the ratios
    fig, ax1 = plt.subplots(figsize=(10, 6))

    # Plot the ratios as a bar plot on the primary axis
    ratios.plot(kind='bar', ax=ax1, color='blue', alpha=0.6)

    # Set the labels for the primary axis
    ax1.set_xlabel('Integrated Intensity Bins')
    
    ax1.set_ylabel('Ratio of Colocalization', color='blue')
    ax1.tick_params(axis='y', labelcolor='blue')

    # Define the bin labels for the x-axis
    bin_labels = [f"{int(left)}-{int(right)}" for left, right in zip(bin_edges[:-1], bin_edges[1:])]
    ax1.set_xticklabels(bin_labels, rotation=45, ha='right', fontsize=7)

    # Create a second y-axis to plot the total count
    ax2 = ax1.twinx()

    # Plot the total count as a line plot on the secondary axis
    total.plot(kind='line', ax=ax2, color='red', marker='o', linewidth=2, markersize=6)

    # Set the labels for the secondary axis
    ax2.set_ylabel('Total Number of Particles', color='red')
    ax2.tick_params(axis='y', labelcolor='red')

    # Set the title for the plot
    ax1.set_title('Ratio of Colocalization and Total Number of Particles vs Integrated Intensity')

    plt.tight_layout()  # Adjust layout to fit all labels


    plt.savefig(str(results_folder_path) + '/' + 'coloc_per_integratedIntensity.pdf', format='pdf')




def integrated_intensity(folder_integrated):
    """
    Output:
    This function creates a histogram of the integrated intensity of the ENTIRE CELL for a channel and saves it in the results folder.
    Input:
    - folder_integrated: The folder that contains the integrated intensity files created from the Integrated_intensity_analysis.py script.
    """

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
    plt.savefig(str(results_folder_path)+"/"+'IntDen_violinplot.pdf', format='pdf')


def particle_distance_against_cell_size(filename, channel, pixel_size):
    global results_folder_path

    csv_files = glob.glob(os.path.join(results_folder_path, '*.csv'))
    big_df=pd.DataFrame()

    for csv_file in csv_files:
        if os.path.basename(csv_file).startswith(filename):
            df = pd.read_csv(csv_file)
            big_df = pd.concat([big_df, df], ignore_index=True) 

    # Create a scatter plot of the relative distance of the particles from the cell center against the cell size, 
    # separate the color of points by the number of particles per cell. If the particle is is a cell that has x number of particles, the color of the point will be different.
    plt.figure(figsize=(10, 6))

    #the number of particles per cell serves as the color of the points except when the number of particles is 0, then ignore the point
    big_df = big_df[big_df['total_count'] != 0]

    # For the same cell id, if there is one particle only in total_count column, put the distance in a negative value 

    #big_df['Distance_from_center'] = big_df['Distance_from_center'].where(big_df['total_count'] == 1, -big_df['Distance_from_center'])
    print(big_df[['total_count', 'Distance_from_center']].head())


    # Apply the condition to negate Distance_from_center where total_count is 1 and it's not already negative
    big_df['Distance_from_center'] = np.where(
        (big_df['total_count'] == 1) & (big_df['Distance_from_center'] > 0),
        -big_df['Distance_from_center'],
        big_df['Distance_from_center']
    )

    print(big_df[['total_count', 'Distance_from_center']].head())

    # Determine the maximum 'Distance_from_center' for each 'cell_id'
    max_distance = big_df.groupby('cell_id')['Distance_from_center'].transform('max')

    # Apply the conditions using np.where
    big_df['Distance_from_center'] = np.where(
        big_df['total_count'] > 1,
        big_df['Distance_from_center'],  # Keep the value unchanged
        np.where(
            (big_df['Distance_from_center'] == max_distance) & (big_df['Distance_from_center'] > 0),
            -big_df['Distance_from_center'],  # Negate the value if it's the maximum and positive
            big_df['Distance_from_center']    # Keep the value unchanged
        )
    )
    plt.scatter(big_df['cell_length'] *  pixel_size, big_df['Distance_from_center'] * pixel_size, c=big_df['total_count'], cmap='viridis', s=10)

    #plot the color legend for the number of particles per cell
        
    #plot two lines to show the cell size considering we want to have the cell length dividec by two and plot it from the center of the cell
    cell_size = big_df['cell_length']
    half_cell_size = cell_size/2

    cell_size_microns= cell_size*pixel_size 
    half_cell_size_microns= half_cell_size*pixel_size

   
    plt.plot(cell_size_microns , half_cell_size_microns, color='red', linestyle='-', linewidth=0.5)
    plt.plot(cell_size_microns , -half_cell_size_microns, color='red', linestyle='-', linewidth=0.5)

    plt.colorbar(label='Number of particles per cell')
    plt.xlabel('Cell Size')
    plt.ylabel('Relative Distance from the Cell Center')
    plt.title('Relative Distance of Particles from the Cell Center vs Cell Size')
    plt.savefig(str(results_folder_path) + '/' + 'distance_against_cell_size.pdf', format='pdf')
    plt.show()



# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
####################################################################################################################################################################################################################
# ----------------------------------------------------------------------# HELPER FUNCTIONS TO CREATE THE LOGIC OF COLOCALIZATIONS---------------------------------------------------------------------------------------------------#
####################################################################################################################################################################################################################
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#

def point_in_poly(x, y, poly):
    """
    Output:
    This function checks if a point (x,y) is inside a polygon defined by a set of vertices.
    Input:
    - x: The x coordinate of the point.
    - y: The y coordinate of the point.
    - poly: The polygon defined by a set of vertices, which is a list of coordinates given by the function 'outlines_list' from the cellpose package.   
    """

    # Check if a point (x,y) is inside a polygon defined by a set of vertices
    n = len(poly)
    inside = False

    if len(poly) < 3:
        # Polygon must have at least 3 points
        return False

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


def cell_length(outline):
    """
    Output:
    - The lenght of the cell (in pixels)
    - The approximate center of the cell (in pixels)
    Input:
    - outline: The polygon defined by a set of vertices, which is a list of coordinates given by the function 'outlines_list' from the cellpose package.
    """

    # initialize the min, max values of the coordinates of the outline of one cell
    maxx = float('-inf')
    maxy = float('-inf')
    minx = float('inf')
    miny = float('inf')

    #This for loop finds the max and min value for x and y coordinates: maxx= max (x) pixel coord, minx= min (x) pixel coord, maxy= max (y) pixel coord, miny= "" ""
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
    cell_len = math.sqrt((maxx-minx)**2+(maxy-miny)**2)
    # Find the approximate center of the cell
    center_coordinate = [minx+(maxx-minx)/2, miny+(maxy-miny)/2]

    return cell_len, center_coordinate, maxx, maxy, minx, miny

import numpy as np
import math

def table_creation(data, header, outlines, size, min, bin_edges):
    """
    Output:
    - Creates a table with the data of the particles per cell for a channel and saves it in the results folder.
    - Multiple lists of number of particles per cell for each cell size are also created. (one, two, three, four, five, six)
    - A dictionary of the coordinates of the particles per unique cell ID is also created. It will be used to calculate the colocalization rate.
    - A list of the relative distance of particles from the cell center is also created.
    - The total number of cells is also calculated.

    Input:
    - data: The data of the particles for a channel.
    - header: The header of the data.
    - outlines: The outlines of the cells.
    - size: The pixel size.
    - min: The minimum cell size accepted to be part of the analysis in pixels.
    """

    bins = [[] for _ in range(len(bin_edges) + 1)]
    counts = np.empty(0, dtype=int)  # The number of particles per cell
    distance_from_center = 0
    distance_from_yaxis = 0
    table = []  # Initialize an empty list for the table
    one, two, three, four, five, six = [[] for _ in range(6)]  # Lists for the counts for each different class of cell length chosen
    rel_dist_center = []  # List of the relative distance of particles from the cell center
    rel_dist_axis = []  # List of the relative distance of particles from the cell center
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
        cell_len, cell_center, maxx, maxy, minx, miny = cell_length(poly)

        if cell_len > min:  # if the cell size is greater than the minimum cell size accepted
            count = 0

            for j in range(len(data)):

                if data is None or data.size == 0:
                    print("No data")
                    continue
                # if data row is empty, skip it

                #if np.isnan(data[j]).any():
                #    print("Nan values")
                #    continue
                # if the particle is outside the cell, skip it
                if data[j][col_x] == 0 and data[j][col_y] == 0:
                    print("Particle outside the cell")
                    continue

                x = data[j][col_x]
                y = data[j][col_y]

                # Find if the particle is inside the cell segmentation
                if point_in_poly(x, y, poly):
                    coord = [x, y]  # if particle inside the cell segmentation
                    if x > cell_center[0]:  # if particle is on the right side of the cell center
                        cell_middle_vector = np.array([maxx - cell_center[0], maxy - cell_center[1]])
                        vector_particle_cell = np.array([x - cell_center[0], y - cell_center[1]])
                    else:  # if particle is on the left side of the cell center
                        cell_middle_vector = np.array([cell_center[0] - minx, cell_center[1] - miny])
                        vector_particle_cell = np.array([cell_center[0] - x, cell_center[1] - y])

                    # do a scalar projection of the vector of the particle and the center of the cell on the vector of the cell
                    scalar_projection = np.dot(vector_particle_cell, cell_middle_vector.T) / np.linalg.norm(cell_middle_vector)

                    # calculate the relative distance of the particle from the cell center
                    distance_from_yaxis = scalar_projection  # can be negative depending if its left or right of the axis of the cell
                    distance_from_center = math.sqrt((y - cell_center[1]) ** 2 + (x - cell_center[0]) ** 2)  # always positive

                    #if the particle is on the right side of the cell center, the distance from the center is positive
                    if x > cell_center[0]:
                        distance_from_center = distance_from_center
                    else:  # if the particle is on the left side of the cell center, the distance from the center is negative
                        distance_from_center = -distance_from_center

                    # append the relative distance of the particle from the cell center to the lists
                    rel_dist_center.append(distance_from_center / (cell_len / 2))  # normalized distance by half the cell length, this gives us a normalized value between 0 and 1.
                    rel_dist_axis.append(abs(distance_from_yaxis) / (cell_len / 2))

                    count += 1  # increment the number of particles per cell

                    row = [
                        i + 1, count, 0, cell_len, x, y,
                        data[j][col_int], cell_center[0], cell_center[1],
                        distance_from_center, distance_from_yaxis
                    ]
                    table.append(row)  # append the row to the table

                    if cellnum not in coloc:
                        coloc[cellnum] = []  # create a list for each cell number in the dictionary
                        coloc[cellnum].append(coord)  # add the coordinates of the particle to the list
                    else:
                        coloc[cellnum].append(coord)

            counts = np.append(counts, count)  # add the number of particles per cell to the list counts

            if count == 0:  # if there is no particle in the cell
                row = [i + 1, 0, 0, cell_len, 0, 0, 0, cell_center[0], cell_center[1], 0, 0]
                table.append(row)  # append the row to the table
                coloc[cellnum] = []  # create a list for each cell number in the dictionary
            else:  # if there is at least one particle in the cell find the cell number and add the number of particles per cell for each cell
                for k in range(count):
                    table[-(k + 1)][2] = count

            cellnum += 1

            nanometer = size * cell_len  # cell_length in pixels, size= pixel size in nanometers


            if nanometer < bin1:
                one.append(count)
            elif bin1 <= nanometer <= bin2:
                two.append(count)
            elif bin2 < nanometer <= bin3:
                three.append(count)
            elif bin3 < nanometer <= bin4:
                four.append(count)
            elif bin4 < nanometer <= bin5:
                five.append(count)
            elif nanometer > bin5:
                six.append(count)

        else:  # cell size is too small (noisy segmentation)
            continue

    table = np.array(table)  # convert the table list to a numpy array
    return table, counts, one, two, three, four, five, six, coloc, rel_dist_center, rel_dist_axis, cellnum


def allowed_dist(cell1, cell2, maxcol_dist):
    """
    Output:
    This function checks if the distance between two particles is less than the maximum distance of colocalization.
    Input:
    - cell1: The coordinates of the first particle.
    - cell2: The coordinates of the second particle.
    - maxcol_dist: The maximum distance of colocalization.
    """

    dist = math.sqrt((cell1[0]-cell2[0])**2+(cell1[1]-cell2[1])**2)
    if (dist > maxcol_dist):
        return False, dist
    else:
        return True, dist


def colocalization(Input_chan, other_chan, pixel_dist, channel):
    """
    Output:
    This function calculates the colocalization rate for a channel. and returns the number of colocalizations in each cell.
    Input:
    - Input_chan: The dictionary of the coordinates of the particles per unique cell ID for the first channel.
    - other_chan: The dictionary of the coordinates of the particles per unique cell ID for the second channel.
    - pixel_dist: The maximum distance of colocalization in pixels.
    """

    result = 0
    boo = False
    candidate = float('inf')
    cols_in_each_cellIDchan1 = defaultdict(lambda: defaultdict(int))
    cols_in_each_cellIDchan2 = defaultdict(lambda: defaultdict(int))

    for key1 in Input_chan.keys():  # iterates through all the cells of the image


        if (len(Input_chan[key1]) != 0 and len(other_chan[key1]) != 0):

            for j,can1 in enumerate(Input_chan[key1]):

                smallest_dist = float('inf')
                if len(other_chan[key1]) != 0:

                    for y, can2 in enumerate(other_chan[key1]):

                        allowed, dist = allowed_dist(
                            Input_chan[key1][j],other_chan[key1][y] , pixel_dist)

                        if allowed and dist <= smallest_dist:
                            smallest_dist = dist
                            candidate = y
                            boo = True
                        else:
                            continue

                    if boo:
                        del other_chan[key1][candidate]
                        # add to the dictionary the number of colocalizations per cell
                        result += 1
                        if channel==1:
                            cols_in_each_cellIDchan1[key1][j+1] +=1
                        else:
                            cols_in_each_cellIDchan2[key1][y+1] +=1
                        
                    boo = False
                else:
                    continue

    return result, cols_in_each_cellIDchan1, cols_in_each_cellIDchan2



def size_colocalization(Input_chan, other_chan, pixel_dist, table, size, bin_edges):
    """
    Output: 
    This function calculates the colocalization rate per cell size for a channel.
    Input:
    - Input_chan: The dictionary of the coordinates of the particles per unique cell ID for the first channel.
    - other_chan: The dictionary of the coordinates of the particles per unique cell ID for the second channel.
    - pixel_dist: The maximum distance of colocalization in pixels.
    - table: The table of the particles per cell for the first channel.
    - size: The pixel size in nanometers.
    """
        # Initialize a list of lists for bins
    bins = [[] for _ in range(len(bin_edges) + 1)]


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
                        try:
                            value = df.set_index('cell_id').loc[key1, 'cell_length']
                        except KeyError:
                            continue

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

# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
####################################################################################################################################################################################################################
# ----------------------------------------------------------------------# MAIN FUNCTION FOR TWO CHANNELS ONLY (COLOCALIZATION)-------------------------------------------------------------------------------------------------#                                                   
####################################################################################################################################################################################################################
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
   

def main1(directory_seg, directory_res, Coloc_bysize , folder_integrated, parameters, channel, selected_folder_path): 
    """
    Output:
    This function creates the results folder and calls all the functions to create the results of the analysis.
    Input:
    - directory_seg: The directory of the segmentation files.
    - directory_res: The directory of the results files.
    - Coloc_bysize: The maximum distance of colocalization in pixels.
    - folder_integrated: The folder that contains the integrated intensity files created from the Integrated_intensity_analysis.py script.
    - parameter: The parameter file.
    - channel: Which channel is the colocalization analysis for.
    - selected_folder_path: The path of the results folder.
    """
    bin_edges = [0, 100, 200, 300, 400, 500, 600] # not used

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
    pixel = parameters['pixel_size']
    too_small = parameters['min_cell_size']  # this is in pixel length 
    names = ['cell_id', 'part_id', 'total_count', 'cell_length',
             'Part_X_(px)', 'Part_Y_(px)', 'IntegratedIntensity', 'Cellcenter_X(px)', 'Cellcenter_Y(px)', 'Distance_from_center', 'Distance_from_yaxis']

# ------ #total particles for each channel --------------------------------------------------------------------------------------------------------------------------------------------
    total_1 = np.empty(0)
    total_2 = np.empty(0)

# ------ Variables to make the histogram for les per cell size------------------------------------------------------------------------------------------------------------------

    tot_1, tot_2, tot_3, tot_4, tot_5, tot_6, tot_1_2, tot_2_2, tot_3_2, tot_4_2, tot_5_2, tot_6_2 = [np.empty(0) for i in range(12)]

# ------Variable to make the histogram of colocalization rate per distance size in pixels----------------------------------------------------------------------------------------------
    # Here you can change the number of distances you want by deleting or adding the desired pixel distance,
    new_dict = {0.5: 0, 1.0: 0, 1.5: 0, 2: 0, 2.5: 0, 3: 0}
    dict_size = {0.5: [], 1.0: [], 1.5: [], 2: [], 2.5: [], 3: []}   
    REL_dist_axis1, REL_dist_axis2, REL_dist_center1, REL_dist_center2 = [np.empty(0) for i in range(4)]


    i = 0
    number_of_cells=0
    for filename in sorted(os.listdir(directory_seg)):

        if not filename.endswith(".DS_Store") and not filename.startswith("._"):
            #print(filename)
            seg_file = os.path.join(directory_seg, filename)
            seg_masks = np.load(seg_file, allow_pickle=True).item()
            #masks = seg_masks['masks']
            outlines = utils.outlines_list(seg_masks['masks'])
            number_of_cells += len(outlines)

            with open(files[2*i], 'r', newline='') as file:
                file_content = file.read().splitlines()
                
            fiji_1= np.genfromtxt(file_content[1:], delimiter='\t')   # input channel 1
            
            if fiji_1.ndim==1: 
                fiji_1 = fiji_1.reshape(1, -1)
            #if file is empty create a fake header to avoid error
            if len(file_content)==0 :
                header = ['cell_id', 'part_id', 'total_count', 'cell_length',
                        'Part_X_(px)', 'Part_Y_(px)', 'IntegratedIntensity', 'Cellcenter_X(px)', 'Cellcenter_Y(px)', 'Distance_from_center', 'Distance_from_yaxis']
            else:
                header = np.genfromtxt([file_content[0]], delimiter='\t', dtype=str)
            
            # if fiji_1 is empty create a single row of zeros to avoid error
            if len(fiji_1)==0:
                fiji_1 = np.zeros((1, 13))

 
            with open(files[2*i+1], 'r', newline='') as file2:
                file_content2 = file2.read().splitlines()
                
            fiji_2= np.genfromtxt(file_content2[1:], delimiter='\t')   # input channel 1
            if fiji_2.ndim==1: 
                fiji_2 = fiji_2.reshape(1, -1)

            if len(file_content2)==0 :
                header = ['cell_id', 'part_id', 'total_count', 'cell_length',
                        'Part_X_(px)', 'Part_Y_(px)', 'IntegratedIntensity', 'Cellcenter_X(px)', 'Cellcenter_Y(px)', 'Distance_from_center', 'Distance_from_yaxis']
            else:
                header = np.genfromtxt([file_content[0]], delimiter='\t', dtype=str)
            
            # if fiji_1 is empty create a single row of zeros to avoid error
            if len(fiji_2)==0:
                fiji_2 = np.zeros((1, 13))

            i += 1

            table_1, counts_1, one, two, three, four,five, six, coloc, rel_dist_center1, rel_dist_axis1, cellnum1 = table_creation(
                fiji_1,header, outlines, pixel, too_small, bin_edges)  # channel 1

            table_2, counts_2, one2, two2, three2, four2,five2, six2, coloc2, rel_dist_center2, rel_dist_axis2, cellnum2 = table_creation(
                fiji_2,header, outlines, pixel, too_small, bin_edges)  # channel 2


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


            REL_dist_axis1 = np.concatenate((rel_dist_axis1, REL_dist_axis1),  axis=None)
            REL_dist_axis2 = np.concatenate((rel_dist_axis2, REL_dist_axis2),  axis=None)
            REL_dist_center1 = np.concatenate((rel_dist_center1, REL_dist_center1),  axis=None)
            REL_dist_center2 = np.concatenate((rel_dist_center2, REL_dist_center2),  axis=None)


# -----------------Creating a dataframe to save the results from the table_creation---------------

            df1 = pd.DataFrame(table_1, columns=names)
            file_path1 = os.path.join(
                selected_folder_path, 'First_channel'+str(i)+'.csv')
            # drop the rows where cell_ID ==0
            df1=df1[df1.cell_id != 0]
            
            df2 = pd.DataFrame(table_2, columns=names)
            file_path2 = os.path.join(
                selected_folder_path, 'Second_channel'+str(i)+'.csv')
            # drop the rows where cell_ID ==0
            df2=df2[df2.cell_id != 0]

#============================COLOCALIZATION============================
#----Iterate through the different distances of colocalization [0.5, 1.0, 1.5, 2, 2.5, 3] and calculate the colocalization rate for each distance
            for pixx in new_dict.keys(): 
                if channel==1:
                    inco = copy.deepcopy(coloc)
                    othco = copy.deepcopy(coloc2)
                    rate, dict1_cel_ID, dict2_cel_ID = colocalization(inco, othco, pixx, channel) 
                    new_dict[pixx] += rate # add the rate of colocalizations to the dictionary
                else: # if channel is 2
                    inco = copy.deepcopy(coloc2)
                    othco = copy.deepcopy(coloc)
                    rate, dict1_cel_ID, dict2_cel_ID = colocalization(inco, othco, pixx, channel)
                    new_dict[pixx] += rate

 #------ We save the number of colocalisation for each cell ID (within distance of 2 pixels (130 nm)) in the table for csv file
                if pixx == 2:  
                    if channel == 1:
                        df1['coloc'] = 0
                        for key1 in dict1_cel_ID.keys():
            
                            condition1 = df1['cell_id'] == key1

                            for id in dict1_cel_ID[key1].keys():
                                
                                condition2 = df1['part_id'] == id
                               
                                print("dict1_cel_ID[key1][id]", dict1_cel_ID[key1][id])

                                df1.loc[condition1 & condition2 , 'coloc'] = dict1_cel_ID[key1][id] 


           #same for channel 2 
                    else:
                        df2['coloc'] = 0

                        for key2 in dict2_cel_ID.keys():
                            condition1 = df2['cell_id'] == key2

                            for id in dict2_cel_ID[key2].keys():
                                condition2 = df2['part_id'] == id
                                df2.loc[condition1 & condition2, 'coloc'] = dict2_cel_ID[key2][id]
           
#------ Calculate the colocalization rate per cell size for each distance of colocalization

            for key in new_dict.keys():  # iterate through the different distances of colocalization
                inco = copy.deepcopy(coloc)
                othco = copy.deepcopy(coloc2)
                res = size_colocalization(inco, othco, key, df1, size=pixel, bin_edges=bin_edges)
                res = np.array(res)
                dict_size[key].append(res) # add the result to the dictionary of the colocalization rate per cell size

# -----------------Saving dataframes created here in CSV file ---------------
 
            df1.to_csv(file_path1, index=True, header=True, sep=',')
            df2.to_csv(file_path2, index=True, header=True, sep=',')


 # --------------Plot the histrogram of counts per cell-------------------------
    counts_histogram(total_1, total_2)

# --------------Plot histogram for each size of cell-------------------------
    size_histogram(tot_1, tot_2, tot_3, tot_4,tot_5,tot_6, 1, total_1)
    size_histogram(tot_1_2, tot_2_2, tot_3_2, tot_4_2, tot_5_2, tot_6_2, 2, total_2)

# -------------Plot histogram for colocalization--------------------------------
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
        
#--------Distance histogram from the center for both channels----------------
    distance_histogram(REL_dist_center1, total_1, channel=1, axis='center')
    distance_histogram(REL_dist_center2, total_2, channel=2, axis='center')

    distance_histogram(REL_dist_axis1, total_1, channel=1, axis='axis')
    distance_histogram(REL_dist_axis2, total_2, channel=2, axis='axis')

    integrated_int_histogram("First_channel", color_channel1) # histogram for integrated intensity
    integrated_int_histogram("Second_channel", color_channel2)   # histogram for integrated intensity
    particle_distance_against_cell_size("First_channel", color_channel1, pixel) # scatter plot for particle distance against cell size
   
    if channel == 1:
        print("First channel")
        coloc_integrated_plot("First_channel") # histogram for colocalization per integrated intensity
    else: 
        print("Second channel")
        coloc_integrated_plot("Second_channel") # histogram for colocalization per integrated intensity



# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
####################################################################################################################################################################################################################
# ----------------------------------------------------------------------# MAIN FUNCTION FOR ONE CHANNEL ONLY----------------------------------------------------------------------------------------------------------#                                                   
####################################################################################################################################################################################################################
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
   
def main2(directory_seg, directory_res, folder_integrated, parameters, selected_folder_path):
    global results_folder_path
    results_folder_path=selected_folder_path
    global number_of_cells
    bin_edges = [0, 100, 200, 300, 400, 500, 600] # not used

    pixel = parameters['pixel_size']
    too_small = parameters['min_cell_size']  # this is in pixel length  
    names = ['cell_id', 'part_id', 'total_count', 'cell_length',
             'Part_X_(px)', 'Part_Y_(px)', 'IntegratedIntensity', 'Cellcenter_X(px)', 'Cellcenter_Y(px)', 'Distance_from_center', 'Distance_from_yaxis']
    files=[]
    for filename in os.listdir(directory_res):
        f = os.path.join(directory_res, filename)
        if not filename.endswith(".DS_Store"):
            files.append(f)
    files.sort()
    print(files)

    if folder_integrated != "":
        integrated_intensity(folder_integrated)

    #--------- nunmber of total particles for each channel (total 2 = empty in this case)--------------------------------------------------------------------------------------------------------------
    total_1=np.empty(0)
    total_2=np.empty(0)
    #-------- Variables to make the histogram for particles per cell size------------------------------------------------------------------------------------------------------------------

    tot_1, tot_2, tot_3, tot_4, tot_5, tot_6 = [np.empty(0) for i in range(6)]

    #-------- Variable to make the histogram of colocalization rate per distance size in pixels----------------------------------------------------------------------------------------------

    dist1=np.empty(0)  # relative distance from the center
    dist2=np.empty(0) # relative distance from the axis

    i=0
    number_of_cells=0
    for filename in sorted(os.listdir(directory_seg)):
    
        if not filename.endswith(".DS_Store"):
            #print(filename)
        
            seg_file = os.path.join(directory_seg, filename)
            seg_masks=np.load(seg_file, allow_pickle=True).item()
            masks=seg_masks['masks']
            outlines= utils.outlines_list(seg_masks['masks'])
            number_of_cells+=len(outlines)

            with open(files[i], 'r', newline='') as file:
                file_content = file.read().splitlines()

            fiji_1= np.genfromtxt(file_content[1:], delimiter='\t')   # input channel 1
            if fiji_1.ndim==1: 
                fiji_1 = fiji_1.reshape(1, -1)
            
            #if file is empty create a fake header to avoid error
            if len(file_content)==0 :
                header = ['cell_id', 'part_id', 'total_count', 'cell_length',
                        'Part_X_(px)', 'Part_Y_(px)', 'IntegratedIntensity', 'Cellcenter_X(px)', 'Cellcenter_Y(px)', 'Distance_from_center', 'Distance_from_yaxis']
            else:
                header = np.genfromtxt([file_content[0]], delimiter='\t', dtype=str)
            
            # if fiji_1 is empty create a single row of zeros to avoid error
            if len(fiji_1)==0:
                fiji_1 = np.zeros((1, 13))

            i+=1

            table_1, counts_1, one, two, three, four, five, six, coloc, rel_distance_center, rel_distance_axis, cellnum=table_creation(fiji_1,header, outlines,pixel, too_small, bin_edges)  # channel 1    
        #------------------num of particle per cell histogram dataset-------------  

            total_1 =np.concatenate((total_1, counts_1), axis=None) 
        #--------------------particle per cell size dataset-----------------------------

            tot_1=np.concatenate((tot_1, one), axis=None) 
            tot_2=np.concatenate((tot_2, two), axis=None) 
            tot_3=np.concatenate((tot_3, three), axis=None) 
            tot_4=np.concatenate((tot_4, four), axis=None)
            tot_5=np.concatenate((tot_5, five), axis=None)
            tot_6=np.concatenate((tot_6, six), axis=None)


        #--------------distance from the center dataset--------------
            dist1=np.concatenate((dist1,rel_distance_center),  axis=None)
            dist2=np.concatenate((dist2,rel_distance_axis),  axis=None)

        #-----------------Saving results of tables created here in CSV file ---------------
            df1 = pd.DataFrame(table_1, columns=names)
            file_path = os.path.join(selected_folder_path, 'table'+str(i)+'.csv')
            df1.to_csv(file_path, index=True, header=True, sep=',')
        
    #--------------Plot the histrogram of counts per cell--------------
    counts_histogram (total_1, total_2)
    #--------------Plot histogram for each size of cell--------------
    size_histogram (tot_1, tot_2, tot_3, tot_4,tot_5, tot_6,  1, total_1)
    #------------Plot histogram for relative distance from the cell center--------------
    distance_histogram(dist1, total_1, channel=1, axis='center')
    distance_histogram(dist2, total_1, channel=1, axis='axis')

    integrated_int_histogram('table', color_channel1)   #input channel 1
    particle_distance_against_cell_size('table', color_channel1, pixel) # scatter plot for particle distance against cell size
