import numpy as np
import matplotlib.pyplot as plt
from glob import glob
from obspy.taup import TauPyModel
import os
from scipy.signal import correlate
from scipy.signal import convolve
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)

#Functions to calculate amplitude ratios

def amplitude_correlation(prem_A, prem_T, model_A):

    amplitude_corr = correlate(prem_A, model_A, mode = 'full')
    xcorr = np.arange(amplitude_corr.size)
    lags = xcorr - (prem_A.size-1)
    distancePerLag = (prem_T[-1] - prem_T[0])/float(prem_T.size)
    offsets = -lags*distancePerLag
    maximum_correlation_indx = np.argmax(amplitude_corr)
    tau = offsets[maximum_correlation_indx]
    maximum_correlation = amplitude_corr[maximum_correlation_indx]

    return tau, maximum_correlation

def amplitude_ratio(prem_A_S, model_A_S, prem_A_SS, model_A_SS, prem_T_S, prem_T_SS):

    A1_S = amplitude_correlation(prem_A_S,prem_T_S, model_A_S)[1]/amplitude_correlation(model_A_S,prem_T_S, model_A_S)[1]
    A2_S = amplitude_correlation(prem_A_S,prem_T_S, prem_A_S)[1]/amplitude_correlation(prem_A_S,prem_T_S, model_A_S)[1]

    A1_SS = amplitude_correlation(prem_A_SS,prem_T_SS, model_A_SS)[1]/amplitude_correlation(model_A_SS,prem_T_SS, model_A_SS)[1]
    A2_SS = amplitude_correlation(prem_A_SS,prem_T_SS, prem_A_SS)[1]/amplitude_correlation(prem_A_SS,prem_T_SS, model_A_SS)[1]

    A1_SS_S = np.min([A1_SS, A2_SS])/np.max([A1_S, A2_S])
    A2_SS_S = np.max([A1_SS, A2_SS])/np.min([A1_S, A2_S])

    Ratio = (1/2)*(A1_SS_S+A2_SS_S)

    return Ratio


#Identify directory where grouped seismic data is located
os.chdir('/home/rwcline/Desktop/5_second_runs_2_percent_gradient')
dirs = np.sort(glob('*'))


#Create empty lists for station_names, PREM seismic models for the corresponding station
#Create empty lists for pertrubed models and the maximum amplitude from each model (used for stacking)
names = []
prem_models = []
models = []
maximum_amplitudes = []
model_names = []
colors = ['red', 'green', 'blue', 'orange', 'cyan', 'purple']
desired_depth = [358, 368, 378, 388]
desired_depth_labels = ['100km','200km','300km', '400km']
desired_gradients = ['02']
dictionary = {}


for d in dirs:
    os.chdir('/home/rwcline/Desktop/5_second_runs_2_percent_gradient/{}' .format(d))
    names.append(d)


#Read in station data and distance from event (deg)
os.chdir('/home/rwcline/Desktop')
station_distance = np.zeros(len(names))

with open('receiver_gll_locs_mtt_m_mpp.dat', 'r') as f:


    lines=f.readlines()
    lines=lines[1:]
    for line in lines:
        parts = line.split()
        stat_name = parts[0].split('_')[0]
        if stat_name in names:
            print(stat_name)
            idx = names.index(stat_name)
            station_distance[idx] = parts[1]

#Group stations and models based on distance from event


sorted_idx = np.argsort(station_distance)
sorted_distances = []
sorted_names = []
sorted_prem = []

for i in sorted_idx:
    sorted_names.append(names[i])
    sorted_distances.append(station_distance[i])

for d in sorted_names:
    os.chdir('/home/rwcline/Desktop/5_second_runs_2_percent_gradient/{}' .format(d))
    files = np.sort(np.array(glob('Gradient*')))
    dictionary[str(d)] = {}
    for f in files:
        file_name = str(f)
        gradient_split = file_name.split('Gradient')[1]
        gradient = gradient_split.split('_')
        if len(gradient) == 5:

            gradient =gradient[1:]
        gradient_name = gradient[0]
        model_names.append(gradient_name)
        dictionary[str(d)][str(gradient_name)] = file_name



models = []
travel_time_delays = np.zeros((len(sorted_names), len(desired_depth)))
travel_time_delays_SS = np.zeros((len(sorted_names), len(desired_depth)))
SS_S_ratio = np.zeros((len(sorted_names), len(desired_depth)))





for i in range(len(sorted_names)):

    print(sorted_names[i])

    for j in desired_gradients:
        os.chdir('/home/rwcline/Desktop/5_second_runs_2_percent_gradient/{}' .format(sorted_names[i]))
        models.append(np.loadtxt(dictionary[sorted_names[i]][j]))

    sorted_prem.append(np.loadtxt(glob('PREM*.dat')[0]))


for n  in range(len(sorted_names)):

    print(n)

    for d in range(len(desired_depth)):

        os.chdir('/home/rwcline/Desktop/5_second_runs_2_percent_gradient/{}' .format(sorted_names[n]))

        name = 'Gradient{}_Depth{}*' .format(desired_gradients[0], desired_depth[d])
        model_name = glob(name)[0]
        print(model_name)
        model = np.loadtxt(model_name)
        model_time = model[:, 0]
        model_amplitude = model[:, 1]

        prem_file_name = glob('*PREM*')[0]
        prem = np.loadtxt(prem_file_name)
        print(prem_file_name)
        prem_time = prem[:, 0]
        prem_amplitude = prem[:, 1]




        Tau_p_times = TauPyModel(model="prem_axisem")
        arrivals = Tau_p_times.get_travel_times(source_depth_in_km=350, distance_in_degree= sorted_distances[n], phase_list=['S', 'SS', 'Sdiff'])
        #print(arrivals)


        if len(arrivals) == 2:

            S_arrival = arrivals[0].time
            SS_arrival = arrivals[1].time

        else:

            S_arrival = arrivals[0].time
            SS_arrival = arrivals[1].time

        window_S = np.where((prem_time >=  S_arrival-40) & (prem_time <=  S_arrival + 40))
        window_SS = np.where((prem_time >=  SS_arrival-40) & (prem_time <=  SS_arrival+ 40))

        prem_time_window_S = prem_time[window_S]
        prem_amplitude_window_S = prem_amplitude[window_S]
        model_amplitude_window_S = model_amplitude[window_S]

        prem_time_window_SS = prem_time[window_SS]
        prem_amplitude_window_SS = prem_amplitude[window_SS]
        model_amplitude_window_SS = model_amplitude[window_SS]

        tau_shift = amplitude_correlation(prem_amplitude_window_S, prem_time_window_S, model_amplitude_window_S)
        tau_shift_SS = amplitude_correlation(prem_amplitude_window_SS, prem_time_window_SS, model_amplitude_window_SS)
        travel_time_delays[n, d] = tau_shift[0]
        travel_time_delays_SS[n, d] = tau_shift_SS[0]

        ratio = amplitude_ratio(prem_amplitude_window_S, model_amplitude_window_S, prem_amplitude_window_SS, model_amplitude_window_SS, prem_time_window_S, prem_time_window_SS)
        dummy_time = np.arange(np.shape(window_SS)[1])
        SS_S_ratio[n, d] = ratio

        print(desired_depth[d])

        plt.plot(prem_time, model_amplitude, color = colors[d], label = r'${}$' .format(desired_depth_labels[d]))

    plt.plot(prem_time, prem_amplitude, 'k')

#Plot single seismogram with multiple models

    plt.xlabel(r'$time (sec)$')
    plt.ylabel(r'$Amplitude$')
    plt.xlim([520, 3600])
    plt.title(r'${}$ $Degrees$' .format(np.round(sorted_distances[n])))

    plt.legend(loc = 0)
    plt.savefig( '/home/rwcline/Desktop/{}_{}_{}.png' .format(desired_gradients[0], desired_depth_labels[d], np.round(sorted_distances[n])))
    plt.cla()
    plt.clf()

for i in range(len(desired_depth)):

    #plt.plot(travel_time_delays[:, i], '-o', label = '{}' .format(desired_depth_labels[i]))
    plt.plot(sorted_distances, travel_time_delays[:, i], '-o', color = colors[i], label = r'${}$' .format(desired_depth_labels[i]), alpha = 0.7, markeredgewidth=0.0)


#Plot delay times

plt.xlabel(r'$Station$ $Location$')
plt.ylabel(r'$\tau_{max}$')
plt.title(r'$2$ $percent$ $Perturbation$')
plt.legend(loc = 0, numpoints=1)
plt.xlim([25, 125])
plt.savefig('/home/rwcline/Desktop/tau_shift_grad02.png')

#Plot amplitude ratios

for i in range(len(desired_depth)):

    plt.plot(sorted_distances, SS_S_ratio[:, i], '-o', color = colors[i], label = r'${}$' .format(desired_depth_labels[i]))

plt.title(r'$2$ $percent$ $Perturbation$')
plt.xlabel(r'$Station$ $Location$')
plt.ylabel(r'$\frac{SS}{S}$ $Ratio$')
plt.legend(loc = 0,  numpoints=1)
plt.xlim([25, 125])
plt.savefig('/home/rwcline/Desktop/SS_S_2percent.png')

plt.show()

mean_SS = np.zeros(len(desired_depth))

for i in range(len(desired_depth)):
    log_data = np.log(SS_S_ratio[i, :])
    mean = np.mean(log_data)
    geometric_mean = np.exp(mean)
    mean_SS[i] = geometric_mean
