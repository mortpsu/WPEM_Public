#!/bin/bash

JOBCOUNT=$1
WORKDIR=$2
RESULTSDIR=$3

set -e 

echo "MATLAB post-run script started for " $JOBCOUNT $WORKDIR $RESULTSDIR

UCINITDIR=$SLURM_SUBMIT_DIR/UC_Init_Solution

cd $WORKDIR
cp week$JOBCOUNT.csv $RESULTSDIR
cp SolverStatus$JOBCOUNT.csv $RESULTSDIR
cp MarginalCost$JOBCOUNT.csv $RESULTSDIR
cp MultipliedDemand$JOBCOUNT.csv $RESULTSDIR
cp UnservedEnergywhole$JOBCOUNT.csv $RESULTSDIR
cp Genperhour$JOBCOUNT.csv $RESULTSDIR
cp NewEnergywholeunserved$JOBCOUNT.csv $RESULTSDIR
cp NewDemandwhole$JOBCOUNT.csv $RESULTSDIR
cp ShedEnergywhole$JOBCOUNT.csv $RESULTSDIR
cp CheckJob$JOBCOUNT.csv $RESULTSDIR
cp results.gdx $RESULTSDIR/results10$JOBCOUNT.gdx
cp results2.gdx $RESULTSDIR/results20$JOBCOUNT.gdx
cp results3.gdx $RESULTSDIR/results30$JOBCOUNT.gdx
cp MarginalElectricityPrice.gdx $RESULTSDIR/MarginalElectricityPrice$JOBCOUNT.gdx
cp TransmissionImpacts.gdx $RESULTSDIR/TransmissionImpacts$JOBCOUNT.gdx
cp GenByBus.gdx $RESULTSDIR/GenByBus$JOBCOUNT.gdx
cp OutageByBus.gdx $RESULTSDIR/OutageByBus$JOBCOUNT.gdx
cp UnitCommitment_p1.gdx $UCINITDIR/UnitCommitment_p1$JOBCOUNT.gdx

cp results4.gdx $RESULTSDIR/results40$JOBCOUNT.gdx
cp results5.gdx $RESULTSDIR/results50$JOBCOUNT.gdx
cp results6.gdx $RESULTSDIR/results60$JOBCOUNT.gdx

echo "MATLAB post-run script complete for " $JOBCOUNT

