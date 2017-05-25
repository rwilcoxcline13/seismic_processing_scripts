import numpy as np
from scipy.signal import correlate
from scipy.signal import convolve
from obspy.taup import TauPyModel


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
