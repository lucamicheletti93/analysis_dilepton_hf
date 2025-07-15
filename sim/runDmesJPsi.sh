NJOBS=160
NCORES=80
NEVENTSPERJOB=5000000
MODE=kMonash
PROCESS=kSoftQCD
SEEDSTART=0
SEEDEND=$(($SEEDSTART + $NJOBS))
parallel -j $NCORES "root -l -q -b 'simulateDmesJPsi.cc($NEVENTSPERJOB, $MODE, $PROCESS, false, 13600, {}, \"jpsidmes/pythia8_jpsidmes_${MODE}_${PROCESS}_seed{}.root\")' > log_jpsidmes_${MODE}_seed{}.txt" ::: $(seq $SEEDSTART $SEEDEND)
