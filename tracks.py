g
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

pixel=int(input("Enter pixel size in nanometer:"))
analysis= input(("Enter B if you want to analyse TracksBound only. Enter L if you want to analyse TracksLifetime only. Enter BL if you want to analyse BOTH:"))
rows_to_add=[]
header=["Cell_ID","Track_IDs", "spt_tr", "spt_widt", "mean_sp", "max_sp", "min_sp", "med_sp", "std_sp", "mean_q", "max_q_tr", "min_q_tr", "med_q_tr", "std_q_tr", "tr_dur", "tr_start", "tr_fin", "x_lc", "y_lc"]





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


def fn(outlines, extracted_data):
    extracted_df = pd.DataFrame(extracted_data)
    print(extracted_df)
    for i in range (len(outlines)):
        
        poly=outlines[i]
      
        for j in range ( len(extracted_data)):
#MAKE SURE TO FIRST OPEN THE RESULTS DATA TO SEE HOW THE COLUMNS WERE SAVED (sometimes FIJI saves the results with one extra empty column at the beginning)
            x=extracted_data[j][16]/(pixel/1000)  #CHANGE THIS SO THE COLUMN IS ALIGNED WITH X_(px)
            y=extracted_data[j][17]/(pixel/1000)    #CHANGE THIS SO THE COLUMN IS ALIGNED WITH Y_(px)
     
            if point_in_poly(x,y,poly): 
       
                row_data = extracted_df.iloc[j].values
                new_list=[i]+row_data
                rows_to_add.append(new_list)
    return rows_to_add
                    
def main(folderTrack, folderSeg):
    boL= False
    boB=False
    boS=False
    files=[]
    j=0
    
    for filename in sorted(os.listdir(folderSeg)):
        f = os.path.join(folderSeg, filename)
        if not filename.endswith(".DS_Store"):
            files.append(f)

    for filename in sorted(os.listdir(folderTrack)):
        if filename.endswith("DS.store"):
            continue

        elif filename.endswith("tracksBound.csv"):

            tr= os.path.join(folderTrack, filename)
            tracksB= np.genfromtxt(tr, delimiter=',', skip_header=0)
            f1= filename
            boB=True
            file_pathB = os.path.join(folderTrack, "trim"+filename+".csv")

        elif filename.endswith(".npy"):

            seg_file = files[i]
            seg_masks=np.load(seg_file, allow_pickle=True).item()

            # masks=seg_masks['masks']
            outlines= utils.outlines_list(seg_masks['masks'])
            #print(outlines)
            n=len(outlines)
            boS=True
            
        elif (filename.endswith("tracksDiffusive.csv")):

            tr= os.path.join(folderTrack, filename)
            tracksL= np.genfromtxt(tr, delimiter=',', skip_header=0)
            f2=filename
            boL=True
            file_pathL = os.path.join(folderTrack, "trim"+filename+".csv")
        
        elif (boB and boS and analysis=="B"):
            j+=1 
            dfB = pd.DataFrame(fn(outlines, tracksB))
            dfB.columns=header
            dfB.to_csv(file_pathB, index=True, header=True, sep=',')

            boB=False
            boS=False
 
        elif (boL and boS and analysis=="L"):
            j+=1
            dfL = pd.DataFrame(fn(outlines, tracksL))
            dfL.columns=header
            dfL.to_csv(file_pathL, index=True, header=True, sep=',')

            boS=False
            boL=False

        elif (boB and boS and boL and analysis=="BL"):
            j+=1
            f1 = re.sub('tracksBound$', '', f1)
            f2= re.sub('tracksLifetime$', '', f2)

            if f1==f2: 

                dfB = pd.DataFrame(fn(outlines, tracksB))
                dfL = pd.DataFrame(fn(outlines, tracksL))

                dfB.columns=header
                dfL.columns=header

                dfB.to_csv(file_pathB, index=True, header=True, sep=',')
                dfL.to_csv(file_pathL, index=True, header=True, sep=',')

                boB=False
                boS=False
                boL=False
            else: 
                continue 

        else:
            continue 




if __name__ == "__main__":
    
    folderTrack= input("enter Analysis folder path:")
    folderTrack = folderTrack.replace('\\', '/')
    folderSeg= input("enter Segmentations folder path:")
    folderSeg = folderSeg.replace('\\', '/')

    main(folderTrack, folderSeg)