import argparse 
import glob
import os
import csv
import numpy as np

parser = argparse.ArgumentParser(description='pick_files')
parser.add_argument('--stat', type=int, required=True)
args = parser.parse_args()
stat = args.stat

ts_path = f"/data/condor_builds/users/avijai/RNO/tutorials-rnog/get_daqstatus/trigger_timestamps/station{stat}"
main_path = f"/data/condor_builds/users/avijai/RNO/tutorials-rnog/get_daqstatus/trigger_thresholds/station{stat}"

week = 604800

def map_ts_run(indir):
    runs = []
    files = sorted(glob.glob(os.path.join(indir, "*")))
    ts_run_map = {}
    for f in files:
        run = f.split("/")[-1].split(".")[0]
        run_no = run_num(run)
        ts = np.load(f)
        ts_conv = ts[0]/(10**(6))*118
        if (ts_conv > 0):
            ts_run_map[ts_conv] = run_no
    return ts_run_map 

def run_num(run):
    nums = []
    for char in run:
        if (char.isdigit() == True):
            nums.append(char)
    num = 0

    for i in range(len(nums)):
        num += float(nums[i])*(10**(len(nums)-i-1))
    return int(num)

ts_run_map = map_ts_run(ts_path)
sorted_ts = sorted(ts_run_map.keys())

runs = []
ts_1 = sorted_ts[0]
for ts in sorted_ts:
    if ts <= ts_1 + week:
        print(ts)
        runs.append(ts_run_map[ts])



np.save(f"/data/condor_builds/users/avijai/RNO_reco/rno_dep/source/NuRadioMC/NuRadioReco/examples/RNO_data/read_data_example/run_lists/runs_{stat}.npy", np.array(runs))
