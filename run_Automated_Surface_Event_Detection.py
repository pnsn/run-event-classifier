#!/home/ahutko/miniconda3/envs/surface/bin/python

# import std packages
from zenodo_get import zenodo_get
import pandas as pd
import os
import sys

# import third party packages
from zenodo_get import zenodo_get
from joblib import dump, load
from obspy.clients.fdsn import Client
client = Client('IRIS')
import obspy

# import event classifier specific packges
from db.get_event_info import unix_to_true_time
from db.get_event_info import get_event_info

#----- get evid from input arguments
try:
    evid = sys.argv[1]
except:
    print("Usage: thisscript.py evid")
    print("evid must be a valid int event_id")
    sys.exit(1)


##### User Defined Parameters for specific case. stride/before/dur in seconds
stride = 5 # Increment/step size for running the model.
before = 90
dur = 210 # Total duration of waveform.

#
orid, ordate, lat, lon, dep, mag, mindist, maxdist, netstas, dists_km, analyst_class = get_event_info(evid)
strordate = ordate.strftime("%Y-%m-%dT%H:%M:%S.%f")[:-4]
stations_id = netstas
event_info = [ evid, orid, ordate, lat, lon, dep, mag, mindist, maxdist, netstas, dists_km, analyst_class ]
starttime = obspy.UTCDateTime(ordate) - before
location = '*'
dists_km = [round(x*10)/10. for x in dists_km]

# Download .joblib models and scalar_params*.csv files. Only do this once.
#doi = '10.5281/zenodo.12629637'
#files = zenodo_get([doi])

##### Specify the model
filename = 'P_50_100_F_1_10_50' # 150s model (default)
#filename = 'P_10_100_F_1_10_50' # 110s model
#filename = 'P_10_30_F_1_15_50' # 40s model

##### These parameters will be automatically defined from the model name, so no need to edit this

mindist = min(dists_km)
if mindist > 50:
    offset = mindist/6. # rough estimate of earliest arrival w.r.t. origin time
    starttime = starttime + ( mindist / 6. )
lowcut = 1  # lower limit of bandpass 
highcut = int(filename.split('_')[5]) # higher limit of bandpass
## window length for moving across the trace (s)
win = int(filename.split('_')[1]) + int(filename.split('_')[2]) # measurement window length
samp_freq = 100
fs = 100
num_corners = 4
original_sr = 100
new_sr = int(filename.split('_')[-1])
scaler_params = pd.read_csv('scaler_params_'+filename+'.csv')
best_model = load('best_rf_model_all_features_'+filename+'.joblib')

##### 
import pandas as pd
import numpy as np
from glob import glob
import obspy

import matplotlib.pyplot as plt
from tqdm import tqdm

from obspy import UTCDateTime
from dateutil import parser
from pytz import timezone
import obspy
from obspy.clients.fdsn.mass_downloader import CircularDomain, \
    Restrictions, MassDownloader

import yaml


# importing the dependencies. 

import scipy as sc
from scipy import signal
import h5py

from obspy.signal.filter import envelope

import tsfel
import random
from datetime import timedelta
import calendar
from tsfel import time_series_features_extractor
from sklearn.preprocessing import StandardScaler

from scipy import stats

#%config InlineBackend.figure_format = "png"

#from Feature_Extraction import compute_hibert

import warnings

# Ignore all warnings
warnings.filterwarnings("ignore")

# displaying all columns from pandas dataframe
# Set display options to show all columns
pd.set_option('display.max_columns', None)

from sklearn.ensemble import RandomForestClassifier
from joblib import dump, load


import time
from datetime import datetime
import seaborn as sns

import sys
sys.path.append('../Common_Scripts')
import seis_feature_new


from common_processing_functions import apply_cosine_taper
from common_processing_functions import butterworth_filter
from common_processing_functions import resample_array

import matplotlib.lines as mlines

import tsfel
cfg_file = tsfel.get_features_by_domain()

import json
import os
from zenodo_get import zenodo_get

from obspy.clients.fdsn import Client
client = Client('IRIS')

def event_detection(starttime = starttime, stations_id = stations_id, dur = dur, stride = stride, new_sr = new_sr, 
                           original_sr = original_sr, lowcut = lowcut, highcut = highcut, num_corners = num_corners, scaler_params = scaler_params, best_model = best_model, location = location, samp_freq = samp_freq, win = win
):
    st_data_full = []
    result_stns = []
    index_stns = []
    prob_stns = []

    st_overall = []
    st_overall_data = []
    st_overall_times = []
    sncls = []

#    for stn_id in tqdm(stations_id): # if you want to see progress bars
    for stn_id in stations_id:

        ## Extracting the station
        stn = stn_id.split('.')[1]

        ## Extracting the network
        network = stn_id.split('.')[0]

        st = []

        # starttime set the starttime. 
        lenst = len(st)
        for channel in [ 'HHZ', 'BHZ', 'EHZ', 'HNZ', 'ENZ' ]:
            if len(st) == lenst:
                try:
                    st += client.get_waveforms(starttime=starttime, endtime=starttime + dur, station=stn,
                                      network=network, channel=channel, location=location)
                    #print("DOING: ", starttime, dur, stn, network, channel, location, sncls)
                    sncl = network + '.' + stn + '.--.' + channel
                    sncls.append(sncl)
                except:
                    pass

        try: 
            
            # detrending and resampling all the data to 100 Hz since thats 
            st = obspy.Stream(st)
            st = st.resample(samp_freq)
            st.detrend()

            st_data_full = []

            ## Ideally the data should come in single stream (len(st) == 1) but if it comes in multiple streams
            ## We will take the first stream. 
            times = st[0].times()

            for i in range(len(st)):
                st_data_full = np.hstack([st_data_full, st[i].data])


                if i+1 < len(st):
                    diff = st[i+1].stats.starttime - st[i].stats.endtime

                    times = np.hstack((times, st[i+1].times()+times[-1]+diff))




                    #print('Final times')

            # st_overall_data will store the full data array for each station. 
            st_overall_data.append(st_data_full)

            # st_overall will store the streams from all stations. 
            st_overall.append(st)

            # so clearly there can be significant differences between the length of st_overall and st_overall_data
            # if data comes in multiple streams. 

            # storing the times. length of times would be equal to number of st_overall_data and length of stations. 
            st_overall_times.append(times)

            # storing the trace data and times. 
#            trace_data = [st_data_full[i:i+int(win*samp_freq)] for i in tqdm(range(0, len(st_data_full), stride*samp_freq)) if len(st_data_full[i:i+int(win*samp_freq)]) == int(win*samp_freq)]
#            trace_times = [times[i] for i in tqdm(range(0, len(st_data_full), stride*samp_freq)) if len(st_data_full[i:i+int(win*samp_freq)]) == int(win*samp_freq)]
            trace_data = [st_data_full[i:i+int(win*samp_freq)] for i in range(0, len(st_data_full), stride*samp_freq) if len(st_data_full[i:i+int(win*samp_freq)]) == int(win*samp_freq)]
            trace_times = [times[i] for i in range(0, len(st_data_full), stride*samp_freq) if len(st_data_full[i:i+int(win*samp_freq)]) == int(win*samp_freq)]


            ## applying processing to the trace. 
            trace_data = np.array(trace_data)
            tapered = apply_cosine_taper(trace_data, 10)
            filtered = butterworth_filter(tapered, lowcut, highcut, original_sr, num_corners, filter_type='bandpass')

            # Applying the normalization.
            norm = filtered / np.max(abs(np.stack(filtered)), axis=1)[:, np.newaxis]

            # Applying resampling
            if new_sr != original_sr:
                norm = np.array([resample_array(arr, original_sr, new_sr) for arr in norm])

     
            result = []
            prob = []
            time = []
            index = []

            for i in range(len(norm)):

                try:
                    tsfel_features = time_series_features_extractor(cfg_file, norm[i], fs= new_sr, verbose = 0)
                    physical_features = seis_feature_new.FeatureCalculator(norm[i],  fs = new_sr).compute_features()

                    final_features = pd.concat([tsfel_features, physical_features], axis=1)
                    columns = scaler_params['Feature'].values

                    features = final_features.loc[:, columns]
                    scaler_params.index = scaler_params['Feature'].values


                    for j in range(len(columns)):
                        features.loc[:, columns[j]] = (features.loc[:, columns[j]] - scaler_params.loc[columns[j],'Mean'])/scaler_params.loc[columns[j], 'Std Dev']

                    features['hod'] = (starttime).hour - 8
                    features['dow'] = (starttime).weekday
                    features['moy'] = (starttime).month



                    #features['E_20_50'] = 0.001
                    # extracting the results.
                    result.append(best_model.predict(features))
                    prob.append(best_model.predict_proba(features))
                    index.append(i)

                except:
                    pass
           




            result_stns.append(result)
            index_stns.append(index)
            prob_stns.append(prob)
            #print("RESULTS: ", result, index, prob, len(prob_stns), len(prob_stns[0]) )
        except:
            pass

    
    
    return result_stns, index_stns, prob_stns, st_overall, st_overall_data, st_overall_times, sncls



st_overall_data = []
st_overall_times = []
st_overall = []
result_stns = []
index_stns = []
prob_stns = []


def plot_detection_results(st_overall_data = st_overall_data, st_overall_times = st_overall_times, st_overall = st_overall, result_stns = result_stns, index_stns = index_stns, prob_stns = prob_stns, xlim = [0,dur], ev_markers = [before,dur], shift = stride, win = win, filename = filename, dists_km = dists_km, event_info = event_info ):

    [ evid, orid, ordate, lat, lon, dep, mag, mindist, maxdist, netstas, dists_km, analyst_class ] = event_info

    plt.rcParams['xtick.labelsize'] = 16  # Font size for xtick labels
    plt.rcParams['ytick.labelsize'] = 20  # Font size for ytick labels

    fig, axs = plt.subplots(len(index_stns)+1, 1, figsize=(15, 3*(len(index_stns)+1)))

    for k in range(len(index_stns)):

        ## This is plotting the normalized data

        trace_data = np.array(st_overall_data[k])
        filtered = butterworth_filter([trace_data], 0.33, 10., samp_freq, num_corners, filter_type='bandpass')
        norm = filtered / np.max(abs(np.stack(filtered)), axis=1)[:, np.newaxis]
        norm = norm[0]

        axs[k].plot(st_overall_times[k], norm, color = 'g' )
#        axs[k].plot(st_overall_times[k], st_overall_data[k] / np.max(abs(st_overall_data[k])))

        ## Setting the title of the plot
        #axs[k].set_title(st_overall[k][0].id, fontsize=20)
        eqprob = np.max(np.array(prob_stns)[k,:,0,0])
        exprob = np.max(np.array(prob_stns)[k,:,0,1])
        noprob = np.max(np.array(prob_stns)[k,:,0,2])
        suprob = np.max(np.array(prob_stns)[k,:,0,3])
        if eqprob > exprob and eqprob > suprob:
            etype = 'EQ'
            etypemax = eqprob
        elif exprob > eqprob and exprob > suprob:
            etype = 'EX'
            etypemax = exprob
        elif suprob > eqprob and suprob > exprob:
            etype = 'SU'
            etypemax = suprob
        else:
            etype = 'NOISE'
            etypemax = noprob
        axs[k].set_title(st_overall[k][0].id + '   EQ: ' + str(eqprob) + '  EX: ' + str(exprob) + '  SU: ' + str(suprob) + '   Prediction: ' + etype + '    Dist= ' + str(dists_km[k]) + ' km', fontsize=20)


        ## These are the colors of detection window. 
        colors = ['black', 'blue', 'white', 'red']
        for i in range(len(index_stns[k])):
            axs[k].axvline(shift * index_stns[k][i] + (win/2), ls='--', color=colors[int(result_stns[k][i])], alpha = 0.6)
            
        # Plot circles on top of the line plot
        for i in range(len(index_stns[k])):
            if result_stns[k][i] == 3:
                axs[k].scatter(shift * np.array(index_stns[k])[i] + (win/2), np.array(prob_stns[k])[:, :, 3][i], ec='k', marker='o', c='red', s=100, zorder=5)
            elif result_stns[k][i] == 0:
                axs[k].scatter(shift * np.array(index_stns[k])[i] + (win/2), np.array(prob_stns[k])[:, :, 0][i], ec='k', marker='o', c='black', s=100, zorder=5)
            elif result_stns[k][i] == 1:
                axs[k].scatter(shift * np.array(index_stns[k])[i] + (win/2), np.array(prob_stns[k])[:, :, 1][i], ec='k', marker='o', c='blue', s=100, zorder=5)
            else:
                axs[k].scatter(shift * np.array(index_stns[k])[i] + (win/2), np.array(prob_stns[k])[:, :, 2][i], ec='k', marker='o', c='white', s=100, zorder=5)

        # Create custom legend for circular markers
        legend_elements = [
           mlines.Line2D([], [], marker='o', color='red', mec = 'k', label='Prob (Su)', markersize=10),
            mlines.Line2D([], [], marker='o', color='k', mec = 'k', label='Prob (Eq)', markersize=10), 
            mlines.Line2D([], [], marker='o', color='white', mec = 'k',  label='Prob (No)', markersize=10),
             mlines.Line2D([], [], marker='o', color='blue', mec = 'k',  label='Prob (Exp)', markersize=10)
        ]
        axs[k].legend(handles=legend_elements, loc='upper right', fontsize=12)

        axs[k].set_xlabel('Time(s) since ' + str(starttime).split('.')[0], fontsize=20)
        axs[k].set_xlim(xlim[0], xlim[1])  # Set x-axis limits if needed
        axs[k].set_ylim(-0.999,1)
        axs[k].axvline(ev_markers[0], ls = '-', c = 'k', lw = 2)
        axs[k].axvline(ev_markers[1], ls = '-', c = 'k', lw = 2)
        ydash = [ 0, 0.25, 0.5, 0.75, 1 ]
        for y_val in ydash:
            axs[k].hlines(y_val, xmin=xlim[0], xmax=xlim[1], colors='gray', linestyles='dotted', alpha = 0.5, linewidth = 0.4 )


    # Finding the size of the biggest array
    max_size = max(max(arr) for arr in index_stns)+1

    # initializing different lists. 
    eq_probs = np.zeros([len(index_stns), max_size])
    exp_probs = np.zeros([len(index_stns), max_size])
    su_probs = np.zeros([len(index_stns), max_size])
    no_probs = np.zeros([len(index_stns), max_size])



    # saving the probabilities. 

    for i in range(len(index_stns)):
        for j in range(len(index_stns[i])):
            
            try:
                #print(f"Processing i={i}, j={j}")
                #print(f"index_stns[i][j] = {index_stns[i][j]}")
                #print(f"prob_stns[i][j] = {prob_stns[i][j]}")

                eq_probs[i, index_stns[i][j]] = prob_stns[i][j][0][0]
                exp_probs[i, index_stns[i][j]] = prob_stns[i][j][0][1]
                no_probs[i, index_stns[i][j]] = prob_stns[i][j][0][2]
                su_probs[i, index_stns[i][j]] = prob_stns[i][j][0][3]

            except IndexError as e:
                            print(f"IndexError: {e}")
                            #print(f"eq_probs shape: {eq_probs.shape}")
                            #print(f"index_stns[i]: {index_stns[i]}")
                            #print(f"prob_stns shape: {len(prob_stns)}, {len(prob_stns[i])}, {len(prob_stns[i][j][0])}")         
            
            
    mean_eq_probs = np.mean(eq_probs, axis = 0)
    mean_exp_probs = np.mean(exp_probs, axis = 0)
    mean_no_probs = np.mean(no_probs, axis = 0)
    mean_su_probs = np.mean(su_probs, axis = 0)
    
    axs[-1].scatter(shift*np.arange(max_size)+(win/2), mean_eq_probs, marker = 'o', c = 'k', s = 100, ec = 'k')
    axs[-1].scatter(shift*np.arange(max_size)+(win/2), mean_exp_probs, marker = 'o', c = 'b', s = 100, ec = 'k')
    axs[-1].scatter(shift*np.arange(max_size)+(win/2),mean_no_probs, marker = 'o', c = 'white', s = 100, ec = 'k')
    axs[-1].scatter(shift*np.arange(max_size)+(win/2),mean_su_probs, marker = 'o', c = 'r', s = 100, ec= 'k')
    axs[-1].set_xlim(xlim[0], xlim[1])
    axs[-1].set_ylim(0, 1)
    axs[-1].set_ylabel('Mean Probability', fontsize = 20)
    axs[-1].legend(handles=legend_elements, loc='upper right', fontsize=12)
    ydash = [ 0, 0.25, 0.5, 0.75, 1 ]
    for y_val in ydash:
        axs[-1].hlines(y_val, xmin=xlim[0], xmax=xlim[1], colors='gray', linestyles='dotted', alpha = 0.5, linewidth = 0.4 )

    # Title at top of figure
    axs[0].text(0,1.5,'Event classification model: '+filename + ' evid: ' + str(evid) + ' ' + mag + ' z=' + str(int(dep)) + ' Analyst: ' + analyst_class, fontsize = 21)
   
    eqmean, exmean, sumean, nomean = 0,0,0,0
    for i in range(len(st_overall)):
        eqprob = np.max(np.array(prob_stns)[i,:,0,0])
        exprob = np.max(np.array(prob_stns)[i,:,0,1])
        noprob = np.max(np.array(prob_stns)[i,:,0,2])
        suprob = np.max(np.array(prob_stns)[i,:,0,3])
        eqmean += eqprob
        exmean += exprob
        sumean += suprob
        nomean += noprob
        if eqprob > exprob and eqprob > suprob:
            etype = 'EQ'
            etypemax = eqprob
        elif exprob > eqprob and exprob > suprob:
            etype = 'EX'
            etypemax = exprob
        elif suprob > eqprob and suprob > exprob:
            etype = 'SU'
            etypemax = suprob
        else:
            etype = 'NOISE1'
            etypemax = noprob
        # Print station-level results to screen
        print("{:<9s} {:5s} {:6.3f} {:8.3f} {:5.1f} {:<4d} {:15s} {:5.3f} {:5.3f} {:5.3f} {:5.3f} {:>6s} {:5.3f}".format(evid, mag, lat, lon, dep, i, sncls[i], eqprob, exprob, noprob, suprob, etype, etypemax ))

    eqmean = round(eqmean*1000/len(st_overall))/1000.
    exmean = round(exmean*1000/len(st_overall))/1000.
    sumean = round(sumean*1000/len(st_overall))/1000.
    if eqmean >= exmean and eqmean >= sumean:
        etype = 'EQ'
    elif exmean > eqmean and exmean > sumean:
        etype = 'EX'
    elif sumean > eqmean and sumean > exmean:
        etype = 'SU'
    if eqmean < 0.5 and exmean < 0.5 and sumean < 0.5:
        etype = 'NOISE2'
    axs[-1].set_title('evid: ' + str(evid) + '  ' + str(mag) + '   Mean prob:   EQ: ' + str(eqmean) + '  EX: ' + str(exmean) + '  SU: ' + str(sumean) + '   Prediction: ' + etype, fontsize = 20)

    # Print event-level results to screen
    print(str(evid) + '  ' + strordate + ' ' + str(mag) + '   Mean probabilities:   EQ: ' + str(eqmean) + '  EX: ' + str(exmean) + '  SU: ' + str(sumean) + '   Prediction: ' + etype + '  Analyst: ' + analyst_class )
    print('')
    plt.tight_layout()  # Adjust subplots to avoid overlap
    plt.savefig('event_classification_' + str(evid) + '_filttraces.png')


##### Run the model.
result_stns, index_stns, prob_stns, st_overall, st_overall_data, st_overall_times, sncls = event_detection(starttime = starttime, stations_id = stations_id, dur = dur, stride = stride, new_sr = new_sr, original_sr = original_sr, lowcut = lowcut, highcut = highcut, num_corners = num_corners, scaler_params = scaler_params, best_model = best_model, location = location, samp_freq = samp_freq, win = win)

##### Plot the results.
plot_detection_results(st_overall_data = st_overall_data, st_overall_times = st_overall_times, 
                       st_overall = st_overall, result_stns = result_stns, index_stns = index_stns, 
                       prob_stns = prob_stns, xlim = [0,dur], ev_markers = [before,dur], shift = stride, win = win, filename = filename, dists_km = dists_km, event_info = event_info )


