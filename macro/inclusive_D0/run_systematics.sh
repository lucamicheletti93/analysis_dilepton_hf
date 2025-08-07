#!/bin/bash

NTIGHT=4
NLOOSE=10
NTRKCUTS=15
declare -a CUTVARIATION=(_altstep1 _altstep2 _altstep3 _widerleft1 _widerleft2 _widerright1 _widerright2 _widerboth _strictleft1 _strictleft2 _strictright1 _strictright2 _strictboth)
NTRIALSRAWY=324

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

# test prompt fraction
if (( ${#CUTVARIATION[@]} > 0 )); then
    for cutvar in "${CUTVARIATION[@]}"; do
        echo -e "\033[32mRUNNING SYSTEMATICS FOR CUT-VARIATION METHOD ${cutvar}\033[0m"
        python3 compute_fraction_cutvar.py systematics/config_cutvar_dzero_mb_y06$cutvar.json
        python3 compute_xsection.py -ir rawyields/rawyields_nocut_dzero_LHC24_JPsiD_y06.root -ie efficiencies/efficiencies_nocutnp_dzero_LHC24k3_trackTuner_ptSmearing1p5_JPsiD_y06.root -ic systematics/promptfrac_dzero_LHC24k3_trackTuner_ptSmearing1p5_JPsiD_y06$cutvar.root -o systematics -s _y06$cutvar -in ../../data_shared/luminosity_dzero_LHC24_minBias_sampled.root
    done
fi

# test signal extraction
if (( NTRIALSRAWY > 0 )); then
    for number in $(seq 0 $((NTRIALSRAWY-1))); do
        echo -e "\033[32mRUNNING SYSTEMATICS FOR RAW-YIELD EXTRACTION ${number}\033[0m"
        python3 compute_xsection.py -ir systematics/multitrial/rawyields_dzero_multitrial.root -sry _${number} -ie efficiencies/efficiencies_nocutnp_dzero_LHC24k3_trackTuner_ptSmearing1p5_JPsiD_y06.root -ic cutvariation/promptfrac_dzero_pp13dot6tev_LHC24_JPsiD_y06.root -o systematics/multitrial -s _y06_$number -in ../../data_shared/luminosity_dzero_LHC24_minBias_sampled.root
    done
fi
