#!/bin/bash

# to be set by the user
#_________________________
ENABLED0=false
ENABLEJPSI=false
ENABLEJPSID0=true
INPUTDATA=input_data_523792.txt
TIMELIMIT=18000000000
CONFIGJSON=configuration_dq_hf.json
NTOT=516
N=51
OUTPUTDIR=LHC22m_pass4_523792
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

NFILESPERCHUNK=$((NTOT / N))
PARTFILE=${INPUTDATA/.txt/.part}
split -d --lines=$NFILESPERCHUNK $INPUTDATA $PARTFILE

WORKFLOW="$TRACKTOCOLLASS $PIDTPCBASEWORKFLOW $PIDTOFBASEWORKFLOW $PIDTPCWORKFLOW $PIDTOFWORKFLOW $HFSKIMCREATORWOFRFLOW $CREATORD0WOFRFLOW $SELECTORD0WOFRFLOW $TASKD0WOFRFLOW $TASKJPSID0WOFRFLOW $FWDTRACKEXT $DQTABLEMAKER $DQTABLEREADER o2-analysis-trackselection $COMMONCONFIGS | o2-analysis-track-propagation $COMMONCONFIGS | o2-analysis-zdc-converter $COMMONCONFIGS | o2-analysis-timestamp $COMMONCONFIGS | o2-analysis-event-selection $COMMONCONFIGS | o2-analysis-ft0-corrected-table $COMMONCONFIGS | o2-analysis-timestamp $COMMONCONFIGS --resources-monitoring 2 --fairmq-ipc-prefix ."
#| o2-analysis-zdc-converter $COMMONCONFIGS

# RUN!

echo "if [ ! -d \"$OUTPUTDIR\" ]; then" > run_workflow.sh
echo "    mkdir $OUTPUTDIR" >> run_workflow.sh
echo "fi" >> run_workflow.sh
echo "declare -a ARRAY" >> run_workflow.sh
echo "for ((i=0; i<=$((N - 1)); i++)); do" >> run_workflow.sh
echo "    formattedi=\$(printf \"%02d\" \$i)" >> run_workflow.sh
echo "    ARRAY[\$i]=\$formattedi" >> run_workflow.sh
echo "    if [ ! -d \"$OUTPUTDIR/\$formattedi\" ]; then" >> run_workflow.sh
echo "        mkdir $OUTPUTDIR/\$formattedi" >> run_workflow.sh
echo "    fi" >> run_workflow.sh
echo "    cp ${INPUTDATA/.txt/.part}\$formattedi $OUTPUTDIR/\$formattedi/" >> run_workflow.sh
echo "    cp $CONFIGJSON $OUTPUTDIR/\$formattedi/$CONFIGJSON" >> run_workflow.sh
echo "    sed -i \"10 s/.*/        \\\"aod-file\\\": \\\"@$PARTFILE\${formattedi}\\\"/\" $OUTPUTDIR/\$formattedi/$CONFIGJSON" >> run_workflow.sh
echo "done" >> run_workflow.sh
echo "parallel -j $N \"cd $OUTPUTDIR/{}; ${WORKFLOW}\" ::: \${ARRAY[@]}" >> run_workflow.sh
echo "rm input_data.part*" >> run_workflow.sh
echo "for ((i=0; i<=$((N - 1)); i++)); do" >> run_workflow.sh
echo "    formattedi=\$(printf \"%02d\" \$i)" >> run_workflow.sh
echo "    rm $OUTPUTDIR/\$formattedi/localhost*" >> run_workflow.sh
echo "done" >> run_workflow.sh
echo "rm localhost*" >> run_workflow.sh
chmod +x run_workflow.sh
