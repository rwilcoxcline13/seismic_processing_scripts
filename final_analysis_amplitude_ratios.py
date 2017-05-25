import numpy as np
import os as os
from glob import glob
from scipy.signal import correlate
from scipy.signal import convolve
from obspy.taup import TauPyModel
import matplotlib.pyplot as plt
from amplitude_ratio_calculator import amplitude_correlation, amplitude_ratio

path_to_data = '/home/rwcline/axisem-9f0be2f/SOLVER/New_Runs_052317/Group_Data_New_Runs_052317'
path_to_script = '/home/rwcline/Desktop/Research/Amplitude_Ratios/Python_scripts'
os.chdir(path_to_data)

#Load station names
stations = []
station_location = np.round(np.loadtxt('receiver_gll_locs_mtt_m_mpp.dat', skiprows=1, usecols=(1,3)))
event_depth = [350]


with open('receiver_gll_locs_mtt_m_mpp.dat') as f:

    next(f)

    for line in f:
        station_name = line.split('_II')[0]
        stations.append(station_name)

    f.close()

PREM_model_name = 'PREM_real_ISO_051717'

desired_stations = stations[::5]
desired_stations = desired_stations[::-1]
desired_location = station_location[ ::5, 0]
desired_location = desired_location[::-1]



desired_component = ['E']

PREM_files = []
Perturbation_files = []
S_arrival_times = np.zeros(len(desired_location))
SS_arrival_times = np.zeros(len(desired_location))


os.chdir(path_to_script)

for distance in range(len(desired_location)):



    Tau_p_times = TauPyModel(model="prem_axisem")
    arrivals = Tau_p_times.get_travel_times(source_depth_in_km=350, distance_in_degree= desired_location[distance], phase_list=['S', 'SS', 'Sdiff'])

    if len(arrivals) == 2:

        S_arrival = arrivals[1].time
        SS_arrival = arrivals[0].time

    else:

        S_arrival = arrivals[0].time
        SS_arrival = arrivals[1].time

    S_arrival_times[distance] = S_arrival
    SS_arrival_times[distance] = SS_arrival

os.chdir(path_to_data)



#specify type of perfurbation

type_of_perturbation = 'Gradient'
perturbation_amount = '02'
perturbation_name = '{}{}' .format(type_of_perturbation, perturbation_amount)

perturbation_depth = np.arange(100000, 800000, 50000)
desired_depths =(3480000+perturbation_depth)/(10000)

SS_S = np.zeros((len(desired_location), len(desired_depths)))
travel_time_delay_S = np.zeros((len(desired_location), len(desired_depths)))
travel_time_delay_SS = np.zeros((len(desired_location), len(desired_depths)))

for i in range(len(desired_stations)):

    print(desired_stations[i])
    PREM_file_names = glob('{}_{}*_{}.dat' .format(desired_stations[i], PREM_model_name, desired_component))[0]
    prem_data = np.loadtxt(PREM_file_names)
    prem_time = prem_data[:, 0]
    prem_amplitude = prem_data[:, 1]

    for d in range(len(desired_depths)):

        for component in desired_component:

            file_name = '{}_{}_Depth{}_{}.dat' .format(desired_stations[i], perturbation_name, str(desired_depths[d])[:3], component)
            print(file_name)
            #Perturbation_files.append(file_name)
            model_data = np.loadtxt(file_name)
            model_time = model_data[:, 0]
            model_amplitude = model_data[:, 1]


            #Calculate Amplitude Ratio
            window_S = np.where((prem_time >=  S_arrival_times[i]-40) & (prem_time <= S_arrival_times[i] + 40))
            window_SS = np.where((prem_time >=  SS_arrival_times[i]-40) & (prem_time <=  SS_arrival_times[i] + 40))

            prem_time_window_S = prem_time[window_S]
            prem_amplitude_window_S = prem_amplitude[window_S]
            model_amplitude_window_S = model_amplitude[window_S]

            prem_time_window_SS = prem_time[window_SS]
            prem_amplitude_window_SS = prem_amplitude[window_SS]
            model_amplitude_window_SS = model_amplitude[window_SS]

            tau_shift = amplitude_correlation(prem_amplitude_window_S, prem_time_window_S, model_amplitude_window_S)
            tau_shift_SS = amplitude_correlation(prem_amplitude_window_SS, prem_time_window_SS, model_amplitude_window_SS)

            travel_time_delay_S[i, d] = tau_shift[0]
            travel_time_delay_SS[i, d] = tau_shift_SS[0]

            ratio = amplitude_ratio(prem_amplitude_window_S, model_amplitude_window_S, prem_amplitude_window_SS, model_amplitude_window_SS, prem_time_window_S, prem_time_window_SS)
            SS_S[i, d] = ratio
