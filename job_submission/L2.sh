#!/bin/bash
  
cmd=$1

source /data/condor_builds/users/avijai/RNO_reco/rno_dep/setup.sh
export PYTHONPATH=/data/condor_builds/users/avijai/RNO_reco/rno_dep/source/NuRadioMC:$PYTHONPATH

echo $LD_LIBRARY_PATH

$cmd

