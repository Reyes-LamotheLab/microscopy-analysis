# Steps to have an almost complete analysis of particles in your Bacteria!

## First: Do cell segmentation using Cellpose and particle detection using FIJI

### **1.1** Particles detection script: <u>particles_detection.ijm</u>
  - This file does particle detection that can be customized for each channel. Then it outputs either a csv or a txt file of the Results created by running the ComDet pluggin from FIJI and saves the results in a folder. 
### **1.2** Cellpose scripts to do cell segmentation:
**Attention**: Make sure you have Cellpose installed on your computer 
If you already have your trained segmentation model and you know how Cellpose works go to b) directly.

####  a) Install CellPose and learn how to annotate your images manually on CellPose GUI:
- Follow this [link](https://github.com/MouseLand/cellpose/tree/main) for documentation on Cellpose
- Follow this [video](https://www.youtube.com/watch?v=3Y1VKcxjNy4) to understand how to use the GUI!
#### b) Either use the default or custom trained model, or train a new model:
  - Cellpose_Train_Model.ipynb
    - Use it to train a custom model for your microscopy images of cells if the default base models provided by CellPose do not segment well.
  - Cellpose_Use_CustomModel.ipynb
    - Don't forget the instructions in the notebook and change the parameters: `diameter` and `model_type` 
    - Use it to segment your images (png, tig, jpeg). It will output the segmentations as _seg.npy files which you need to have to do the analysis of particles in cells (use the app_analysis particles.py).
      
## Second: Use app or the config.toml
### Run this script app_analysis_particles.py either on Visual Studio Code or on your terminal
**Attention:** This can be used ONLY after you ran <u>`particles_detection.ijm`</u> and <u>`Cellpose_Use_CustomModel.ipynb`</u> since it relies on their output.
  - Use this analysis script to ouput a table per image analyse and histograms that can be useful such as the following:
    - Cells Length Distribution
    - Particle's Relative Distance from the Center of the Cell 
    - Counts of Particles per Cell
    - Colocalization for any Channel
    - Colocalization per Cell Length for any Channel
    - Comparison of the Ratio of Colocalization and the Total Number of Particles against the Integrated Intensity

## Third: Use the other scripts if you want more 
You can also use the output of ***`integrated_intensity.py`*** to get a swarm plot of the integrated intensities for all cells of multiple images
### integrated_intensity.py
  - This file will use the output of the segmentation done by Cellpose as well as the TIF images in Fiji to be able to give us the Integrated intensity of each cell for images and then outputs a result file as csv format for each image.
### tracking_fullframe.py
  - Use this file if you have the tracking results from FIJI, it ouputs files outputB.mat for spotsbound and spotslifetime outputL.mat that can be analysed using Pablo's MatLab script.
    
# Output example

- ***Histograms created by the app***
![My project-1](https://github.com/Reyes-LamotheLab/microscopy-analysis/assets/83682336/fec1c7eb-dec2-4217-8be1-cdf02b20eed8)
![My project-2](https://github.com/Reyes-LamotheLab/microscopy-analysis/assets/83682336/77a5a3d8-2859-4f3c-8c41-ddb2e2abaa65)

- ***Table created for each image***

![image](https://github.com/Reyes-LamotheLab/microscopy-analysis/assets/83682336/60c2691b-1450-4411-bc2a-5a6023f95737)

