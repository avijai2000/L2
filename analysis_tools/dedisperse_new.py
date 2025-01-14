import csv
import numpy as np
from scipy import interpolate
import waveform_utilities as wu
import time 
import matplotlib.pyplot as plt
import math
from scipy.signal import savgol_filter
from NuRadioReco.utilities import units
from NuRadioReco.detector import detector
import logging
import datetime
import pandas as pd

class Dedisperse:

    def __init__(self):

        self.path = "resp_felix.csv"
        self.channels = [0,1,2,3,5,6,7,9,10,22,23]
        self.freq = np.linspace(50, 1000, 1000) * units.MHz 

    def load_phase_response_as_spline(self):
        resp = {}
        for ch in self.channels:
            resp[ch] = []
        
        count = 0
        with open(self.path) as f:
            reader = csv.reader(f)
            for row in reader:
                if (count > 0):
                    for i in range(1, len(row)):
                        resp[self.channels[i-1]].append(np.angle(complex((row[i]))))
                count += 1
    

        phase_splines = {}
        
        for ch in self.channels:
            phs_unwrapped = np.unwrap(resp[ch]) # unwrapped phase in radians
            
            
            the_phase_spline = interpolate.Akima1DInterpolator(
                self.freq, phs_unwrapped,
                method="makima",
            )
            
            the_phase_spline.extrapolate = False
            
            phase_splines[ch] = the_phase_spline
            
        return phase_splines

    def eval_splined_phases(self, freqs_to_evaluate, phase_splines, ch):
        phase_spline = phase_splines[ch]
        these_phases = savgol_filter(phase_spline(freqs_to_evaluate), 40, 2)
        these_phases = np.nan_to_num(these_phases) # convert nans to zeros
        return these_phases

    def run(self, event, station, phase_splines):
        for channel in station.iter_channels():
            if (channel.get_id() in self.channels):
                volts = channel.get_trace()
                times = channel.get_times()
                
                if len(times) != len(volts):
                    raise Exception("The time and volts arrays are mismatched in length. Abort.")

                # first thing to do is get the frequency domain representation of the trace
                
                freqs, spectrum = wu.time2freq(times, volts)
                
                # interpolate the *unwrapped phases* to the correct frequency base
                phased_interpolated = self.eval_splined_phases(freqs, phase_splines, channel.get_id())
            
                # convert these into a complex number
                phased_rewrapped = np.exp((0 + 1j)*phased_interpolated)
                
                # do complex division to do the dedispersion

                spectrum /= phased_rewrapped

                # back to the time domain
                times, volts = wu.freq2time(times, spectrum)                    

                channel.set_trace(volts, channel.get_sampling_rate())




