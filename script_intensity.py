import os
from ij import IJ
from ij.plugin.frame import RoiManager
from ij.gui import PolygonRoi
from ij.gui import Roi
from java.awt import FileDialog
from ij.measure import Measurements, ResultsTable
from ij.plugin.filter import Analyzer

# Image directory
image_dir = "/Users/ceciliaacosta/Desktop/data_Viri/FIJI/train+model_Viri/integrated/images/"  # Replace with the directory containing the images
image_files = sorted([file_name for file_name in os.listdir(image_dir) if file_name.endswith(".png")])

# Iterate over image files in the directory
for image_file in image_files:
	
    if image_file.endswith(".png"):
        image_path = os.path.join(image_dir, image_file)
        imp = IJ.openImage(image_path)
        imp.show()


        # Open text file
        text_file_dir = "/Users/ceciliaacosta/Desktop/data_Viri/FIJI/train+model_Viri/integrated/outlines/"  # Replace with the directory containing the text files
        text_file_prefix = image_file[:-4]  # Remove the ".tif" extension from the image file name
        text_file_suffix = "_cp_outlines.txt"  # Suffix to match the text file names
        text_file_name = [file_name for file_name in os.listdir(text_file_dir) if file_name.startswith(text_file_prefix) and file_name.endswith(text_file_suffix)]
        RM = RoiManager()

        rm = RM.getRoiManager()
        if len(text_file_name) == 1:
            text_file_path = os.path.join(text_file_dir, text_file_name[0])
            file_name=text_file_path
            print(image_path)
            print(file_name)

            # Read coordinates from text file and create ROIs
            textfile = open(file_name, "r")
            for line in textfile:
                xy = map(int, line.rstrip().split(","))
                X = list(xy)[::2]
                Y = list(xy)[1::2]
                imp.setRoi(PolygonRoi(X, Y, Roi.POLYGON))
                # IJ.run(imp, "Convex Hull", "")
                roi = imp.getRoi()
                rm.addRoi(roi)

            textfile.close()
            rm.runCommand("Associate", "true")
            rm.runCommand("Show All with labels")
            # Measure data from ROIs
            rm.runCommand(imp, "Measure")
            # Get the measurement results from the RoiManager
            rt = ResultsTable.getResultsTable()
			
            results_dir = "/Users/ceciliaacosta/Desktop/data_Viri/FIJI/train+model_Viri/integrated/results/"  # Replace with the directory to save the results
            results_file = image_file.replace(".png", ".csv")  # Assuming results file names correspond to image names
            results_path = os.path.join(results_dir, results_file)
            
            rt.saveAs(results_path)
            rm.runCommand("Delete")  # Delete all ROIs in the ROI Manager

            RM.reset()
            rt.reset()
            imp.close()


        else:
            print("Text file not found or multiple matching text files found.")


print("All images processed.")