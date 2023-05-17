#!/bin/bash

# to be set by the user
#_________________________
ENABLED0=false
ENABLEJPSI=false
ENABLEJPSID0=true
INPUTDATA=input_data.txt
TIMELIMIT=180
CONFIGJSON=configuration_dq_hf.json
#_________________________
COMMONCONFIGS="-b --configuration json://$CONFIGJSON --aod-memory-rate-limit 1000000000 --shm-segment-size 7500000000 --time-limit $TIMELIMIT"

# HF + HELPERS
TRACKTOCOLLASS=""
PIDTOFBASEWORKFLOW=""
PIDTPCBASEWORKFLOW=""
PIDTOFWORKFLOW=""
PIDTPCWORKFLOW=""
HFSKIMCREATORWOFRFLOW=""
CREATORD0WOFRFLOW=""
SELECTORD0WOFRFLOW=""
TASKD0WOFRFLOW=""
TASKJPSID0WOFRFLOW=""
# DQ + HELPERS
FWDTRACKEXT=""
DQTABLEMAKER=""
DQTABLEREADER=""
if [ "$ENABLED0" = true ] || [ "$ENABLEJPSID0" = true ]; then
    TRACKTOCOLLASS="o2-analysis-track-to-collision-associator $COMMONCONFIGS |"
    PIDTOFBASEWORKFLOW="o2-analysis-pid-tof-base $COMMONCONFIGS |"
    PIDTPCBASEWORKFLOW="o2-analysis-pid-tpc-base $COMMONCONFIGS |"
    PIDTOFWORKFLOW="o2-analysis-pid-tof-full $COMMONCONFIGS |"
    PIDTPCWORKFLOW="o2-analysis-pid-tpc-full $COMMONCONFIGS |"
    HFSKIMCREATORWOFRFLOW="o2-analysis-hf-track-index-skim-creator $COMMONCONFIGS |"
    CREATORD0WOFRFLOW="o2-analysis-hf-candidate-creator-2prong $COMMONCONFIGS |"
    SELECTORD0WOFRFLOW="o2-analysis-hf-candidate-selector-d0 $COMMONCONFIGS |"
fi
if [ "$ENABLEJPSI" = true ] || [ "$ENABLEJPSID0" = true ]; then
    FWDTRACKEXT="o2-analysis-fwdtrackextension $COMMONCONFIGS |"
    DQTABLEMAKER="o2-analysis-dq-table-maker $COMMONCONFIGS |"
    DQTABLEREADER="o2-analysis-dq-table-reader $COMMONCONFIGS |"
fi
if [ "$ENABLED0" = true ]; then
    TASKD0WOFRFLOW="o2-analysis-hf-task-d0 $COMMONCONFIGS |"
fi
if [ "$ENABLEJPSID0" = true ]; then
    TASKD0WOFRFLOW="o2-analysis-dq-task-jpsi-hf $COMMONCONFIGS |"
fi

WORKFLOW="$TRACKTOCOLLASS $PIDTPCBASEWORKFLOW $PIDTOFBASEWORKFLOW $PIDTPCWORKFLOW $PIDTOFWORKFLOW $HFSKIMCREATORWOFRFLOW $CREATORD0WOFRFLOW $SELECTORD0WOFRFLOW $TASKD0WOFRFLOW $TASKJPSID0WOFRFLOW $FWDTRACKEXT $DQTABLEMAKER $DQTABLEREADER o2-analysis-trackselection $COMMONCONFIGS | o2-analysis-track-propagation $COMMONCONFIGS | o2-analysis-multiplicity-table $COMMONCONFIGS | o2-analysis-timestamp $COMMONCONFIGS | o2-analysis-zdc-converter $COMMONCONFIGS | o2-analysis-event-selection $COMMONCONFIGS | o2-analysis-timestamp $COMMONCONFIGS --aod-file @$INPUTDATA --resources-monitoring 2 --fairmq-ipc-prefix ."

# RUN!
echo $WORKFLOW > run.sh
chmod +x run.sh
./run.sh
rm run.sh
rm localhost*
rm dpl-config.json