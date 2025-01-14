import ROOT
import os
import logging
import importlib.resources as pkg_resources
import yaml
import config_files
import numpy as np
import array
import matplotlib.pyplot as plt

ROOT.gSystem.Load(os.environ.get("RNO_G_DEPS_INSTALL_DIR")+"/lib/libRootFftwWrapper.so")
ROOT.gInterpreter.Declare(f'#include "{os.environ.get("RNO_G_DEPS_INSTALL_DIR")}/include/FFTtools.h"')

class CW_Filter:

    def run(self, event, station, station_id, analysis_config):

        file = pkg_resources.open_text(config_files, 
                                   "analysis_configs.yaml")
        file_content = yaml.safe_load(file)

        cw_filters = {}

        try:
            this_station_config = file_content[f"station{station_id}"][f"config{analysis_config}"]
        except:
            logging.error(f"Could not find station {station_id}, config {analysis_config} in the cw config file")
            raise
    
        if this_station_config["filters"] is not None:
            for filter_name, config_settings in this_station_config["filters"].items():
                the_filter = ROOT.FFTtools.SineSubtract(3,
                                                   config_settings["min_power_ratio"],
                                                   False)
                the_filter.setVerbose(False)
                the_filter.setFreqLimits(config_settings["min_freq"], 
                                    config_settings["max_freq"])
                ROOT.SetOwnership(the_filter, True) # give python full ownership
                cw_filters[filter_name] = the_filter

        file.close() 
        
        all_times = {}
        all_volts = {}

        for channel in station.iter_channels():
            volts = channel.get_trace()
            times = channel.get_times()
            graph = ROOT.TGraph(len(times), times, volts)

            for filter_i, filt in cw_filters.items():
                new_graph = filt.subtractCW(graph, -1)
                all_times[channel.get_id()] = new_graph.GetX()
                all_volts[channel.get_id()] = new_graph.GetY()
            
            channel.set_trace(np.array(new_graph.GetY()), channel.get_sampling_rate())

        #return all_times, all_volts


class Bandpass:

    def __init__(self):

        self.minfreq = 0.05 #GHz
        self.maxfreq = 0.7
        self.interp_tstep = 0.2 #ns

    def run(self, event, station):

        dummy = ROOT.FFTtools.ButterworthFilter(ROOT.FFTtools.LOWPASS, 2, 100)
        ROOT.SetOwnership(dummy, True) # give python full ownership

        nyquist = 1./(2.*self.interp_tstep*1E-9) # interp speed in seconds
        freq_lopass = 700E6/nyquist # lowpass fitler at 700 MHz, in units of nyquist
        freq_hipass = 50E6/nyquist # highpass filter at 50 MHz, in units of nyquist

        lp = ROOT.FFTtools.ButterworthFilter(ROOT.FFTtools.LOWPASS, 5, freq_lopass)
        hp = ROOT.FFTtools.ButterworthFilter(ROOT.FFTtools.HIGHPASS, 5, freq_hipass)
        ROOT.SetOwnership(lp, True) # give python full ownership
        ROOT.SetOwnership(hp, True) # give python full ownership
        
        for channel in station.iter_channels():
            volts = channel.get_trace()
            times = channel.get_times()
            volts_lp = array.array("d", [0]*len(volts))
            volts_hp = array.array("d", [0]*len(volts))

            lp.filterOut(len(times), volts, volts_lp)
            hp.filterOut(len(times), volts_lp, volts_hp)

            channel.set_trace(np.array(volts_hp), channel.get_sampling_rate())



class Interpolate:

    def __init__(self):
        self.interp_tstep = 0.2 #ns 

    def run(self, event, station):
        all_times = {}
        for channel in station.iter_channels():

            volts = channel.get_trace()
            times = channel.get_times()
            all_times[channel.get_id()] = times 

            graph = ROOT.TGraph(len(times), times, volts)

            interp_wave = ROOT.FFTtools.getInterpolatedGraph(graph,self.interp_tstep)
            ROOT.SetOwnership(interp_wave, True) # give python full ownership
            
            volts_interp = np.array(interp_wave.GetY())
            time_interp = np.array(interp_wave.GetX())
            if (len(volts_interp) % 2 != 0):
                volts_interp = np.append(volts_interp, 0)
                time_interp = np.append(time_interp, time_interp[-1] + self.interp_tstep)
                
            channel.set_trace(volts_interp, channel.get_sampling_rate())
            channel.set_trace_start_time(0)

