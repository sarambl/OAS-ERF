import os
import socket
from useful_scit.util.make_folders import make_folders
import pandas as pd
from oas_erf.project_root import get_project_base

##################################################################
######################## EDIT HERE: ##############################
# SPECIFY DIRECTORY WHERE PROJECT IS PLACED :
# Replace next line by
# user_spec_base_path = '/path/to/dir/'
user_spec_base_path = None
# Specify folder where model output is located
raw_data_path_NorESM = None
##################################################################


project_name = 'OAS-ERF'
hostname = socket.gethostname()
# Raw data path:
nird_project_code = 'NS9066K'
_data_folder = 'model_output/archive/'

if user_spec_base_path is None:
    project_base_path = get_project_base(hostname, nird_project_code)
else:
    project_base_path = user_spec_base_path
if raw_data_path_NorESM is None:
    raw_data_path_NorESM = project_base_path + _data_folder

pathdic_raw_data = {'NorESM': raw_data_path_NorESM}  # [file_source]}


def get_input_datapath(model='NorESM'):
    return pathdic_raw_data[model]


# Plots path:
path_plots = project_base_path + '/Plots_' + project_name + '/'

paths_plotsave = dict(maps=path_plots + 'maps/',
                      comparison=path_plots + 'global_comparison/',
                      one_value=path_plots + 'one_value/',
                      lineprofiles=path_plots + 'lineprofiles/',
                      sizedist=path_plots + 'sizedistribution/',
                      sizedist_time=path_plots + 'sizedist_time/',
                      levlat=path_plots + 'levlat/',
                      eusaar=path_plots + 'eusaar/'
                      )


def get_plotpath(key):
    """

    :param key:
    :return:
    """
    if key in paths_plotsave:
        return paths_plotsave[key]

    else:
        return path_plots + '/' + key + '/'


path_eusaar_data = project_base_path + '/EUSAAR_data'

path_EBAS_data = project_base_path + '/EBAS_data'

# Output data:

path_outdata = project_base_path + '/Output_data_' + project_name + '/'

path_eusaar_outdata = path_eusaar_data + '/EUSAAR_data/'

latlon_path = path_outdata + 'latlon.nc'


def get_outdata_base():
    return path_outdata


outpaths = dict(
    pressure_coords=path_outdata + '/Fields_pressure_coordinates',
    original_coords=path_outdata + '/computed_fields_ng',
    computed_fields_ng=path_outdata + '/computed_fields_ng',  # native grid computed fields
    pressure_coords_converstion_fields=path_outdata + '/Pressure_coordinates_conversion_fields',
    pressure_density_path=path_outdata + '/Pressure_density',
    masks=path_outdata + '/means/masks/',
    area_means=path_outdata + '/means/area_means/',
    map_means=path_outdata + '/means/map_means/',
    levlat_means=path_outdata + '/means/levlat_means/',
    profile_means=path_outdata + '/means/profile_means/',
    sizedistrib_files=path_outdata + '/sizedistrib_files',
    collocated=path_outdata + '/collocated_ds/',
    eusaar=path_outdata + '/eusaar/',
)


def get_outdata_path(key):
    if key in outpaths:
        return outpaths[key]
    else:
        print('WARNING: key not found in outpaths, constants.py')
        return path_outdata + '/' + key


make_folders(path_outdata)

# data info
proj_lc = project_name.lower().replace('-', '_')

path_data_info = project_base_path + f'{project_name}/{proj_lc}/data_info/'
# output locations:
path_locations_file = path_data_info + 'locations.csv'
if os.path.isfile(path_locations_file):
    collocate_locations = pd.read_csv(path_locations_file, index_col=0)
else:
    _dic = dict(Hyytiala={'lat': 61.51, 'lon': 24.17},
                Melpitz={'lat': 51.32, 'lon': 12.56},
                Amazonas={'lat': -3., 'lon': -63.},
                Beijing={'lat': 40, 'lon': 116})
    collocate_locations = pd.DataFrame.from_dict(_dic)
    collocate_locations.to_csv(path_locations_file)

