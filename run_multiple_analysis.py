import json
from analysis_non_app import *  # Import all functions from your analysis functions file

# Load configuration
with open('config.json', 'r') as config_file:
    config = json.load(config_file)

def run_analysis(analysis_name, params, output_folder):
    """
    Executes an analysis with given parameters and stores results in the specified folder.
    """
    # Example function call - replace with actual function calls as needed
    
    main1(directory_seg, directory_res, Coloc_bysize , folder_integrated, parameter, channel, selected_folder_path): 
    # Assuming your function accepts keyword arguments
    # Save result to the specified output folder
    # This could involve saving files, plots, etc.
    save_result(output_folder, result)

def save_result(output_folder, result):
    # Implement saving results to files or directories here
    pass

# Iterate over each analysis in the configuration and run them
for analysis_name, details in config.items():
    run_analysis(analysis_name, details['params'], details['output_folder'])