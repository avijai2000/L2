import dedisperse_new
from NuRadioReco.modules.RNO_G import channelBlockOffsetFitter
import filters
import calibrate

class DataLoading:

    def run(self, event, station, station_id, analysis_config, phase_splines):

        dedisperse = dedisperse_new.Dedisperse()
        cw_filter = filters.CW_Filter()
        bp_filter = filters.Bandpass()
        interpolate = filters.Interpolate()
        calib = calibrate.Calibrate()
        block_offsets = channelBlockOffsetFitter.channelBlockOffsets()

        block_offsets.remove_offsets(event, station)
        calib.run(event, station)
        interpolate.run(event, station)
        dedisperse.run(event, station, phase_splines)
        cw_filter.run(event, station, station_id, analysis_config)
        bp_filter.run(event, station)

