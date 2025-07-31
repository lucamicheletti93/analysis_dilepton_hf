#!/bin/bash

NTIGHT=0
NLOOSE=0
NTRKCUTS=15

# test tighter selections
if (( NTIGHT > 0 )); then
    for number in $(seq 0 $((NTIGHT-1))); do
        echo -e "\033[32mRUNNING SYSTEMATICS FOR TIGHT CUT ${number}\033[0m"
        python3 compute_efficiencies.py -c systematics/config_dzero_mb_y06_tight$number.yml -p -e
        python3 extract_raw_yields.py -c systematics/config_dzero_mb_y06_tight$number.yml -p -t -f
        python3 compute_fraction_cutvar.py systematics/config_cutvar_dzero_mb_y06_tight$number.json
        python3 compute_xsection.py -ir systematics/rawyields_nocut_dzero_LHC24_JPsiD_y06_tight$number.root -ie systematics/efficiencies_nocutnp_dzero_LHC24k3_trackTuner_ptSmearing1p5_JPsiD_y06_tight$number.root -ic systematics/promptfrac_dzero_LHC24k3_trackTuner_ptSmearing1p5_JPsiD_y06_tight$number.root -o systematics -s _y06_tight$number
    done
fi

# test looser selections
if (( NLOOSE > 0 )); then
    for number in $(seq 0 $((NLOOSE-1))); do
        echo -e "\033[32mRUNNING SYSTEMATICS FOR LOOSE CUT ${number}\033[0m"
        python3 compute_efficiencies.py -c systematics/config_dzero_mb_y06_loose$number.yml -p -e
        python3 extract_raw_yields.py -c systematics/config_dzero_mb_y06_loose$number.yml -p -t -f
        python3 compute_fraction_cutvar.py systematics/config_cutvar_dzero_mb_y06_loose$number.json
        python3 compute_xsection.py -ir systematics/rawyields_nocut_dzero_LHC24_JPsiD_y06_loose$number.root -ie systematics/efficiencies_nocutnp_dzero_LHC24k3_trackTuner_ptSmearing1p5_JPsiD_y06_loose$number.root -ic systematics/promptfrac_dzero_LHC24k3_trackTuner_ptSmearing1p5_JPsiD_y06_loose$number.root -o systematics -s _y06_loose$number
    done
fi

# test track selections
if (( NTRKCUTS > 0 )); then
    for number in $(seq 0 $((NTRKCUTS-1))); do
        echo -e "\033[32mRUNNING SYSTEMATICS FOR TRACK CUT ${number}\033[0m"
        python3 compute_efficiencies.py -c systematics/config_dzero_mb_y06_trkcut$number.yml -p -e
        python3 extract_raw_yields.py -c systematics/config_dzero_mb_y06_trkcut$number.yml -p -t -f
        python3 compute_fraction_cutvar.py systematics/config_cutvar_dzero_mb_y06_trkcut$number.json
        python3 compute_xsection.py -ir systematics/rawyields_nocut_dzero_LHC24_JPsiD_y06_trkcut$number.root -ie systematics/efficiencies_nocutnp_dzero_LHC24k3_trackTuner_ptSmearing1p5_JPsiD_y06_trkcut$number.root -ic systematics/promptfrac_dzero_LHC24k3_trackTuner_ptSmearing1p5_JPsiD_y06_trkcut$number.root -o systematics -s _y06_trkcut$number -in ../../data_shared/luminosity_dzero_LHC24_minBias_sampled_trksys.root
    done
fi