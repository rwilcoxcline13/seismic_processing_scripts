import numpy as np
import os as os
import glob as glob
import pandas as pd
import shutil as shutil

directory_of_new_runs = 'New_Runs_052317'
final_data_directory_name = 'Group_Data_New_Runs_052317'
path = '/home/rwcline/axisem-9f0be2f/SOLVER/{}' .format(directory_of_new_runs)
final_data_directory_path = '{}/{}' .format(path, final_data_directory_name)
os.chdir(path)
#os.mkdir(final_data_directory_name)

#Store directories of runs in an array
Runs = os.listdir(path)
Runs = np.sort(Runs)


#Load station names
stations = []


with open('receiver_gll_locs_mtt_m_mpp.dat') as f:

    next(f)

    for line in f:
        station_name = line.split('_II')[0]
        stations.append(station_name)

    f.close()

station_location = np.round(np.loadtxt('receiver_gll_locs_mtt_m_mpp.dat', skiprows=1, usecols=(1,3)))
paths_to_data = []

#move processed data to new directory

for r in Runs:

    if r != 'receiver_gll_locs_mtt_m_mpp.dat' and r != 'Group_Data_New_Runs_052317':
        path_to_data = '{}/{}/Data_Postprocessing/SEISMOGRAMS' .format(path, r)
        print(path_to_data)
        paths_to_data.append(path_to_data)

for p in paths_to_data:

    for station in stations:

        station_files = '{}/{}' .format(p, station)
        os.chdir(station_files)
        station_data = os.listdir(station_files)

        for data in station_data:

            name = data.split('_')[0]
            component = data.split('_')[-1]
            run_name = p.split('/')[6]
            new_file_name = '{}_{}_{}' .format(name, run_name, component)
            final_destination = '{}/{}' .format(final_data_directory_path, new_file_name)
            shutil.move("{}" .format(data), "{}" .format(final_destination) )
