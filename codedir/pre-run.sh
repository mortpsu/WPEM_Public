#!/bin/bash

#debug info

#echo "SLURM_SUBMIT_DIR:" $SLURM_SUBMIT_DIR
#echo "Job arg 1 and 2:" $1 $2

set -e

JOBCOUNT=$1
WORKDIR=$2

echo "Input script started for " $JOBCOUNT $WORKDIR

#COALDIR=$/Coal_Cali_Data
UCINITDIR=$SLURM_SUBMIT_DIR/UC_Init_Solution

cp -r $SLURM_SUBMIT_DIR/Gams_data/* $WORKDIR
#cp $COALDIR/CoalBinaryCali$JOBCOUNT.gdx $COALDIR/CoalGenCali$JOBCOUNT.gdx $WORKDIR
cp $UCINITDIR/UnitCommitment_p1$JOBCOUNT.gdx $WORKDIR
cp $WORKDIR/UnitCommitment_p1$JOBCOUNT.gdx $WORKDIR/UnitCommitment_p1.gdx
#cp $WORKDIR/CoalBinaryCali$JOBCOUNT.gdx $WORKDIR/CoalBinaryCali.gdx
#cp $WORKDIR/CoalGenCali$JOBCOUNT.gdx $WORKDIR/CoalGenCali.gdx

echo "Input script complete for " $JOBCOUNT

