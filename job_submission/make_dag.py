#!/usr/env python 
import numpy as np
import os

dag = ""
main_dir = "/data/condor_builds/users/avijai/RNO_reco/rno_dep/source/NuRadioMC/NuRadioReco/examples/RNO_data/read_data_example/run_lists/runs_11.npy"

stations = [11]
runs = list(np.load(main_dir))

for station in stations:
    main_path = f"/data/i3store/users/avijai/station_runs/station{station}"
    for run in runs:
        run_path = main_path + "/run" + str(run)
        job_name = f"{station}_{run}_L2"
        cmd = f"JOB {job_name} L2.sub\n"
        cmd += f"VARS {job_name} job_name=\"{job_name}\""
        cmd += f" cmd=\"'python3 /data/condor_builds/users/avijai/RNO_reco/rno_dep/source/NuRadioMC/NuRadioReco/examples/RNO_data/read_data_example/L2.py "
        cmd += f"--file {run_path} "
        cmd += f"--stat {station} "
        cmd += "'"
        cmd += f"\"\n"
        dag += cmd

open("L2.dag", 'w').write(dag)



