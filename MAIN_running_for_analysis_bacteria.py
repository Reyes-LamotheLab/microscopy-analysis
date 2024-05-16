import toml
from analysis_particles_Fn_only import *  # Make sure to define your analysis functions here

# Load TOML configuration
with open('config.toml', 'r') as config_file:
    config = toml.load(config_file)

# Extracting specific sections into variables for easier access
path_config = config['path']
coloc_config = config['coloc']
parameters = config['parameters']


def perform_analysis(path_config, coloc_config, parameters):
    """
    Executes the analysis based on the provided configuration.
    """
    # Example of using the path configuration
    segmentation_path = path_config['segmentation']
    particles_path = path_config['particles']

    output_folder = path_config['output_folder']
    bins=parameters["bins"] 
    # Example of handling colocalization settings
    if coloc_config['colocalization']:
        
        if coloc_config['colocalization_per_size']:
        # Implement colocalization analysis here
            main1(segmentation_path, particles_path, True , coloc_config["folder_integrated"] , parameters, parameters["channel"], output_folder)
        else:
            main1(segmentation_path, particles_path, False , coloc_config["folder_integrated"] , parameters, parameters["channel"], output_folder)
    else: 
        print("Colocalization is not enabled")
        main2(segmentation_path, particles_path,coloc_config["folder_integrated"], parameters, output_folder)



# Implement the analysis functions here

def save_result(output_folder, config):
    # Implement saving tolm here
    with open(output_folder + '/result.toml', 'w') as result_file:
        toml.dump(config, result_file)

# Perform the analysis
perform_analysis(path_config, coloc_config, parameters)
save_result(path_config['output_folder'], config)