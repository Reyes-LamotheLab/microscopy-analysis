# microscopy-analysis
## Description of each script
### particles_detection.ijm
  - This file does particle detection that can be customized for each channel. Then is outputs either a csv or a txt file of the Results created by running the ComDet pluggin from FIJI and saves the results in a folder. 
  - Make sure to open the results and identify how many columns there are since it is important for the Analysis_with_seg_user code
### integrated_intensity.py
  - This file will use the output of the segmentation done by Cellpose to be able to give us the Integrated intensity of each cell for images and then outputs a result file as csv format for each image.
### Cellpose scripts to segment 
  - Cellpose_Train_Model.ipynb
    - use it to train a custom model for your microscopy images of cells if the base models provided by them do not segment well
  - Cellpose_Use_CustomModel.ipynb
    - use it to segment your images (png, tig, jpeg). It will output the segmentations as _seg.npy files which you need to have to do the analysis of particles in cells
### app_analysis_particles.py 
  - Use this analysis script to ouput a table per image analyse and histograms that can be useful such as the following:
    - Cells Length Distribution
    - Particle's Relative Distance from the Center of the Cell 
    - Counts of Particles per Cell
    - Colocalization for any Channel
    - Colocalization per Cell Length for any Channel
![image](https://github.com/Reyes-LamotheLab/microscopy-analysis/assets/83682336/bf86997d-6648-4a23-97af-855584cd3e1d)
