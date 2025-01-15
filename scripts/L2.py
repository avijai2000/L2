import argparse
import NuRadioReco
import matplotlib.pyplot as plt
from NuRadioReco.modules.io.RNO_G.readRNOGDataMattak import readRNOGData
import pandas as pd
import numpy as np
from NuRadioReco.utilities import units
from NuRadioReco.modules import channelBandPassFilter
from NuRadioReco.detector import detector
import datetime
from NuRadioReco.modules import sphericalWaveFitter
from NuRadioReco.modules import channelAddCableDelay
import snr 
import rpr
import reco
import csw 
import hilbert 
import impulsivity 
import time 
import utils
import reco_utils
import itertools
from NuRadioReco.utilities import fft as fft_reco
import config_files
import data_loading 
import glitch_removal
import dedisperse_new 

parser = argparse.ArgumentParser(description='L2')
parser.add_argument('--file', type=str, required=True)
parser.add_argument('--stat', type=int, required=True)
parser.add_argument('--year', type=int, required=True)
args = parser.parse_args()
filename = args.file
station_id = args.stat
year = args.year

detectorpath = "/data/i3store/users/avijai/RNO_season_2023.json"
channels_PA = [0,1,2,3]
channels_PS = [0,1,2,3,5,6,7]
channels_all = [0,1,2,3,5,6,7,9,10,22,23]
do_envelope = True
res = 100
solution = "direct_ice"

list_of_root_files = [filename]
run_no = filename.split("/")[-1]
rpr = rpr.RPR()
csw = csw.CSW()
reco = reco.Reco()
snr = snr.SNR()
hilbert = hilbert.Hilbert()
impulsivity = impulsivity.Impulsivity()
data_loading = data_loading.DataLoading()
glitch = glitch_removal.GlitchFinder()

readRNOGData = NuRadioReco.modules.io.RNO_G.readRNOGDataMattak.readRNOGData()
readRNOGData.begin(list_of_root_files, mattak_kwargs = {"backend":"pyroot"})

mappath = reco.build_travel_time_maps(detectorpath, station_id, channels_all) 
ttcs = utils.load_ttcs(mappath, channels_all)
csw_channels = [channels_PA, channels_PS, channels_all]
csw_chan_names = ["PA", "PS", "ALL"]
freqs = np.arange(50, 1001) / 1000
dedisperse = dedisperse_new.Dedisperse()
phase_splines = dedisperse.load_phase_response_as_spline()
analysis_config = year 

L2 = {}
for i_event, event in enumerate(readRNOGData.run()): 
    L2[event.get_id()] = {}
    
    station = event.get_station(station_id)
    is_bad = glitch.run(event, station)

    data_loading.run(event, station, station_id, analysis_config, phase_splines)

    maxcorr_point, maxcorr, score, t_ab, surf_corr_ratio, max_surf_corr = reco.run(event, station, detectorpath, station_id, channels_all, do_envelope, res, mappath, ttcs)
    
    avg_snr, rms = snr.run(event, station)
    avg_rpr = rpr.run(event, station, rms)
    csw_rpr = {}
    csw_snr = {}
    csw_hilbert_snr = {}
    csw_impulsivity = {}
    
    for i in range(len(csw_channels)):
        chans = csw_channels[i]
        csw_times, csw_values = csw.run(event, station, detectorpath, station_id, chans, solution, ttcs, maxcorr_point, maxcorr, score, t_ab)
        csw_snr[csw_chan_names[i]] = snr.get_snr_single(csw_times, csw_values)
        csw_rpr[csw_chan_names[i]] = rpr.get_single_rpr(csw_times, csw_values)
        csw_hilbert_snr[csw_chan_names[i]] = hilbert.hilbert_snr(csw_values)
        csw_impulsivity[csw_chan_names[i]] = impulsivity.calculate_impulsivity_measures(csw_values)
    
    L2[event.get_id()]["glitch"] = is_bad 
    L2[event.get_id()]["max_correlation"] = maxcorr_point 
    L2[event.get_id()]["reco_vars"] = maxcorr 
    L2[event.get_id()]["avg_snr"] = avg_snr
    L2[event.get_id()]["avg_rpr"] = avg_rpr
    L2[event.get_id()]["csw_snr"] = csw_snr
    L2[event.get_id()]["csw_rpr"] = csw_rpr
    L2[event.get_id()]["csw_hilbert_snr"] = csw_hilbert_snr
    L2[event.get_id()]["csw_impulsivity"] = csw_impulsivity
    L2[event.get_id()]["surf_corr_ratio"] = surf_corr_ratio
    L2[event.get_id()]["max_surf_correlation"] = max_surf_corr 
    

df = pd.DataFrame(L2).transpose()
df.to_csv(f'/data/condor_shared/users/avijai/RNO_reco/rno_dep/source/NuRadioMC/NuRadioReco/examples/RNO_data/read_data_example/L2_files_improved/L2_{station_id}_{run_no}.csv')



    
