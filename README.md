# microscopy-analysis
## Description of each script
### particles_detection.ijm
  - This file does particle detection that can be customized for each channel. Then is outputs either a csv or a txt file of the Results created by running the ComDet pluggin from FIJI and saves the results in a folder. 
  - Make sure to open the results and identify how many columns there are since it is important for the Analysis_with_seg_user code
### integrated_intensity.py
  - This file will use the output of the segmentation done by Cellpose as well as the TIF images in Fiji to be able to give us the Integrated intensity of each cell for images and then outputs a result file as csv format for each image.
### Cellpose scripts to do cell segmentation
  - Cellpose_Train_Model.ipynb
    - use it to train a custom model for your microscopy images of cells if the base models provided by them do not segment well
  - Cellpose_Use_CustomModel.ipynb
    - use it to segment your images (png, tig, jpeg). It will output the segmentations as _seg.npy files which you need to have to do the analysis of particles in cells
### app_analysis_particles.py 
  - Use this analysis script to ouput a table per image analyse and histograms that can be useful such as the following:
  - This can be used ONLY after you ran ***`particles_detection.ijm`*** and ***`Cellpose_Use_CustomModel.ipynb`*** since it relies on their output.
  - You can also use the output of ***`integrated_intensity.py`*** to get a swarm plot of the integrated intensities for all cells of multiple images
    - Cells Length Distribution
    - Particle's Relative Distance from the Center of the Cell 
    - Counts of Particles per Cell
    - Colocalization for any Channel
    - Colocalization per Cell Length for any Channel
   
### Output example

- ***Histograms created by the app***
![My project-1](https://github.com/Reyes-LamotheLab/microscopy-analysis/assets/83682336/fec1c7eb-dec2-4217-8be1-cdf02b20eed8)
![My project-2](https://github.com/Reyes-LamotheLab/microscopy-analysis/assets/83682336/77a5a3d8-2859-4f3c-8c41-ddb2e2abaa65)

- ***Table created for each image***

![image](https://github.com/Reyes-LamotheLab/microscopy-analysis/assets/83682336/60c2691b-1450-4411-bc2a-5a6023f95737)

