
import os
import pandas as pd

import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt
from cellpose import plot, utils
import statistics as s 
import sys, os
import math
import re
import csv
import scipy 
from natsort import natsorted

#pixel=int(input("Enter pixel size in nanometer:"))
analysis= input(("Enter B if you want to analyse spotsBound only. Enter L if you want to analyse spotsLifetime only. Enter BL if you want to analyse BOTH:"))
header=["tr_identifi", "tr_fram", "x_tr", "y_tr", "inten"]
dic_B={}     #dictonnarty with Key: file name (imgID_ Cell_ID) Value: ( majorAxis, minorAxis ,pix list(outlines), cellID, imageID)
dic_L={}   
# result_path=input(("enter the path where you want your results:"))
# result_path=result_path.replace('\\', '/')

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

def cell_length(outline):
    #---initialize the min, max values of the coordinates of the outline of one cell
    maxx= float('-inf')
    maxy= float('-inf')
    minx= float('inf')
    miny= float('inf')

#--This for loop finds the max and min value for x and y coordinates: maxx= max (x) pixel coord, minx= min (x) pixel coord, maxy= max (y) pixel coord, miny= "" ""
    for i in range (len(outline)):

        x=outline[i][0]
        y=outline[i][1]
        
        if x>maxx:
            maxx=x
        if y>maxy:
            maxy=y
        if y<miny:
            miny=y
        if x<minx:
            minx=x

    length=math.sqrt((maxx-minx)**2+(maxy-miny)**2) #cell length
    majoraxis=max(maxx-minx, maxy-miny) #major axis
    minoraxis=min(maxx-minx, maxy-miny) #minor axis

    return length, majoraxis, minoraxis



def fn(outlines, extracted_data, file_path, extension, imID, B_or_L, spots_or_tracks):
    extracted_df = pd.DataFrame(extracted_data)

    for i in range (len(outlines)):

        rows_to_add = []
        poly=outlines[i]
        length, majoraxis, minoraxis= cell_length(poly)
        if length>12:
            for j in range ( len(extracted_df)):

    #MAKE SURE TO FIRST OPEN THE RESULTS DATA TO SEE HOW THE COLUMNS WERE SAVED (sometimes FIJI saves the results with one extra empty column at the beginning)
                if spots_or_tracks =="spots":
                    x=extracted_data[j][2] #CHANGE THIS SO THE COLUMN IS ALIGNED WITH X_(px)
                    y=extracted_data[j][3]    #CHANGE THIS SO THE COLUMN IS ALIGNED WITH Y_(px)
                
                else: 
                    x=extracted_data[j][16] #CHANGE THIS SO THE COLUMN IS ALIGNED WITH X_(px)
                    y=extracted_data[j][17]

                if point_in_poly(x,y,poly): #verifies if the track is in the cell outline
                
                    row_data = extracted_df.iloc[j].values
                    row= np.append(row_data, i)  #add cell id at the end of the row 
                    rows_to_add.append(row)

            #create the table for each cellID and save it as a csv file
            if len(rows_to_add)==0:
                continue

            else: 


                path=file_path +"_"+str(i)+ extension
                with open( path, "w") as f:
          
                    writer = csv.writer(f)
                    writer.writerows(rows_to_add)



            if spots_or_tracks=="spots": 
        #THIS IS ONLY FOR SPOTS TO CREATE A DICTIONNARY                 
                file_name = os.path.basename(path)  # Extract the file name from the path
                file_name_without_extension = os.path.splitext(file_name)[0]  # Remove the file extension

                desired_part = file_name_without_extension.split('.tif')[0]

                if B_or_L=="B":   
                    dic_B[desired_part]=[]
                    dic_B[desired_part].extend([majoraxis, minoraxis, poly, i, imID] )
                
                else:
                    dic_L[desired_part]=[]
                    dic_L[desired_part].extend([majoraxis, minoraxis, poly, i, imID] )
            else:
                continue

        else:
            continue
            
    return 


                    
def main(folderTrack, folderSeg):
    SpotL= False
    SpotB=False
    TrL=False
    TrB=False

    files=[]
    j=0
    m=0


    new_directory = folderTrack+"/Result_tracks/"
    # Create the new directory if it doesn't exist
    if not os.path.exists(new_directory):
        os.makedirs(new_directory)


    for filename in natsorted(os.listdir(folderSeg)):
        f = os.path.join(folderSeg, filename)
        if not filename.endswith(".DS_Store"):
            files.append(f)


    for filename in natsorted(os.listdir(folderTrack)):
        m+=1
   
        if filename.endswith("DS.store"):
            continue

        if (SpotB and SpotL and TrB and TrL and analysis=="BL"):

            if f1==f2==f3==f4: 
                
                seg_file = files[j]
                print(seg_file)
                seg_masks=np.load(seg_file, allow_pickle=True).item()
                outlines= utils.outlines_list(seg_masks['masks'])
                
                fn(outlines, spotsB,  fpath_SpotB, ".tif_spotsBound.csv", j+1, "B", "spots")
                fn(outlines, spotsL, fpath_SpotL, ".tif_spotsLifetime.csv", j+1, "L", "spots")
                fn(outlines, tracksB, fpath_TrB, ".tif_tracksBound.csv", j+1, "B", "tracks")
                fn(outlines, tracksL, fpath_TrL, ".tif_tracksLifetime.csv", j+1, "L", "tracks")

                SpotB=False
                SpotL=False
                TrB=False
                TrL=False

                j+=1

  
        if filename.endswith("spotsBound.csv"):

            tr= os.path.join(folderTrack, filename)
            spotsB= np.genfromtxt(tr, delimiter=',', skip_header=0)  #reads the files 
      
            f1= re.sub('.tif_spotsBound.csv$', '', filename)
            SpotB=True
            fpath_SpotB = os.path.join(new_directory, f1)#to save the new results
   
        if (filename.endswith("spotsLifetime.csv")):

            tr= os.path.join(folderTrack, filename)
            spotsL= np.genfromtxt(tr, delimiter=',', skip_header=0)

            f2= re.sub('.tif_spotsLifetime.csv$', '', filename)
            SpotL=True
            fpath_SpotL = os.path.join(new_directory, f2) #to save the new results
  
        if filename.endswith("tracksBound.csv"):

            tr= os.path.join(folderTrack, filename)
            tracksB= np.genfromtxt(tr, delimiter=',', skip_header=0)  #reads the files 
      
            f3= re.sub('.tif_tracksBound.csv$', '', filename)
            TrB=True
            fpath_TrB = os.path.join(new_directory, f1)#to save the new results
  
        if (filename.endswith("tracksLifetime.csv")):

            tr= os.path.join(folderTrack, filename)
            tracksL= np.genfromtxt(tr, delimiter=',', skip_header=0)

            f4= re.sub('.tif_tracksLifetime.csv$', '', filename)
            TrL=True
            fpath_TrL = os.path.join(new_directory, f2) #to save the new results
             
        # elif (SpotB and analysis=="B"):
                    
        #     seg_file = files[j]
        #     seg_masks=np.load(seg_file, allow_pickle=True).item()
        #     outlines= utils.outlines_list(seg_masks['masks'])

        #     j+=1 
        #     dfB = pd.DataFrame(fn(outlines, spotsB, fpath_SpotB, ".tif_spotsBound.csv"))
        #     dfB.to_csv(fpath_SpotB, index=False, header=False, sep=',')
        #     SpotB=False
  
 
        # elif (SpotL and analysis=="L"):

        #     seg_file = files[j]
        #     seg_masks=np.load(seg_file, allow_pickle=True).item()
        #     outlines= utils.outlines_list(seg_masks['masks'])

        #     j+=1
        #     dfL = pd.DataFrame(fn(outlines, spotsL, fpath_SpotL,".tif_spotsLifetime.csv" ))
        #     dfL.to_csv(fpath_SpotL, index=False, header=False, sep=',')
        #     SpotL=False

        
        if (SpotB and SpotL and TrB and TrL and analysis=="BL"):
          
            if f1==f2==f3==f4: 
                
                seg_file = files[j]
                print(seg_file)
                seg_masks=np.load(seg_file, allow_pickle=True).item()
                outlines= utils.outlines_list(seg_masks['masks'])
                
                fn(outlines, spotsB,  fpath_SpotB, ".tif_spotsBound.csv", j+1, "B", "spots")
                fn(outlines, spotsL, fpath_SpotL, ".tif_spotsLifetime.csv", j+1, "L", "spots")
                fn(outlines, tracksB, fpath_TrB, ".tif_tracksBound.csv", j+1, "B", "tracks")
                fn(outlines, tracksL, fpath_TrL, ".tif_tracksLifetime.csv", j+1, "L", "tracks")

                SpotB=False
                SpotL=False
                TrB=False
                TrL=False

                j+=1

        else:
            continue 


if __name__ == "__main__":
    
    folderTrack= input("enter Analysis folder path:")
    folderTrack = folderTrack.replace('\\', '/')
    folderSeg= input("enter Segmentations folder path:")
    folderSeg = folderSeg.replace('\\', '/')

    main(folderTrack, folderSeg)


new_directory = folderTrack+"/Result_tracks/"
fpath_SpotB = os.path.join(new_directory, 'outputB.mat')
fpath_SpotL = os.path.join(new_directory, 'outputL.mat')
scipy.io.savemat(fpath_SpotL, dic_L)
scipy.io.savemat(fpath_SpotB, dic_B)
