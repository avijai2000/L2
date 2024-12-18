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
import dedisperse_new
import csw 
import hilbert 
import impulsivity 
import time 
import utils
import reco_utils
import itertools

parser = argparse.ArgumentParser(description='L2')
parser.add_argument('--file', type=str, required=True)
parser.add_argument('--stat', type=int, required=True)
args = parser.parse_args()
filename = args.file
station_id = args.stat

#station_id = 12
detectorpath = "/data/i3store/users/avijai/RNO_season_2023.json"
channels_PA = [0,1,2,3]
channels_PS = [0,1,2,3,5,6,7]
channels_all = [0,1,2,3,5,6,7,9,10,22,23]
do_envelope = True
res = 100
solution = "direct_ice"

#list_of_root_files = ['/data/i3store/users/avijai/station_runs/station12/run1611']
list_of_root_files = [filename]
run_no = filename.split("/")[-1]
rpr = rpr.RPR()
csw = csw.CSW()
reco = reco.Reco()
snr = snr.SNR()
dedisperse = dedisperse_new.Dedisperse()
hilbert = hilbert.Hilbert()
impulsivity = impulsivity.Impulsivity()

readRNOGData = NuRadioReco.modules.io.RNO_G.readRNOGDataMattak.readRNOGData()
readRNOGData.begin(list_of_root_files, mattak_kwargs = {"backend":"pyroot"})

mappath = reco.build_travel_time_maps(detectorpath, station_id, channels_all) 
ttcs = utils.load_ttcs(mappath, channels_all)
#mappath2 = reco.build_travel_time_maps(detectorpath, station_id, channels_PS)
#ttcs2 = utils.load_ttcs(mappath, channels_PS)
csw_channels = [channels_PA, channels_PS, channels_all]
csw_chan_names = ["PA", "PS", "ALL"]
freqs = np.arange(50, 1001) / 1000
dedisperse = dedisperse_new.Dedisperse()
phase_splines = dedisperse.load_phase_response_as_spline()
 
L2 = {}
for i_event, event in enumerate(readRNOGData.run()): 
    tstart = time.time()
    L2[event.get_id()] = {}
    station_id = event.get_station_ids()[0]
    station = event.get_station(station_id)
    maxcorr_point, maxcorr, score, t_ab = reco.run(event, station, detectorpath, station_id, channels_all, do_envelope, res, mappath, ttcs)
    tstop = time.time()
    #print(tstop-tstart, "total")

    
    avg_snr, rms = snr.run(event, station)
    avg_rpr = rpr.run(event, station, rms)
    #can change channels for CSW
    csw_rpr = {}
    csw_snr = {}
    csw_hilbert_snr = {}
    csw_impulsivity = {}
    
    tstart = time.time()
    for i in range(len(csw_channels)):
        #tstart = time.time()
        chans = csw_channels[i]
        csw_times, csw_values = csw.run(event, station, detectorpath, station_id, chans, solution, ttcs, maxcorr_point, maxcorr, phase_splines, score, t_ab)
        #tstop = time.time()
        #print(tstop-tstart, "csw", csw_chan_names[i])
        csw_snr[csw_chan_names[i]] = snr.get_snr_single(csw_times, csw_values)
        csw_rpr[csw_chan_names[i]] = rpr.get_single_rpr(csw_times, csw_values)
        csw_hilbert_snr[csw_chan_names[i]] = hilbert.hilbert_snr(csw_values)
        csw_impulsivity[csw_chan_names[i]] = impulsivity.calculate_impulsivity_measures(csw_values)
        
    
    tstop = time.time()
    #print(tstop-tstart, "total")

    
    L2[event.get_id()]["max_correlation"] = maxcorr_point 
    L2[event.get_id()]["reco_vars"] = maxcorr 
    L2[event.get_id()]["avg_snr"] = avg_snr
    L2[event.get_id()]["avg_rpr"] = avg_rpr
    L2[event.get_id()]["csw_snr"] = csw_snr
    L2[event.get_id()]["csw_rpr"] = csw_rpr
    L2[event.get_id()]["csw_hilbert_snr"] = csw_hilbert_snr
    L2[event.get_id()]["csw_impulsivity"] = csw_impulsivity
    

df = pd.DataFrame(L2)
df.to_csv(f'/data/condor_builds/users/avijai/RNO_reco/rno_dep/source/NuRadioMC/NuRadioReco/examples/RNO_data/read_data_example/L2_files/L2_{station_id}_{run_no}.csv')



    
