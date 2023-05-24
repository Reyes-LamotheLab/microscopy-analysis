# microscopy-analysis
## Description of each script
### Particle_Detection_thereal.ijm
  - This file does particle detection that can be customized for each channel. Then is outputs either a csv or a txt file of the Results created by running the ComDet pluggin from FIJI and saves the results in a folder. 
  - Make sure to open the results and identify how many columns there are since it is important for the Analysis_with_seg_user code
### Script_intensity.py
  - This file will use the output of the segmentation done by Cellpose to be able to give us the Integrated intensity of each cell for images and then outputs a result file as csv format for each image.
