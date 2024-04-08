#! /bin/bash

IterationCount=$1

echo "Starting convergence run for iteration Count $IterationCount" 

cd $SLURM_SUBMIT_DIR

module load gams/44.3-s
export LD_PRELOAD=/storage/icds/RISE/sw8/gams/gams44.3_linux_x64_64_sfx/libstdc++.so.6

#cd results
#mkdir Iteration${IterationCount}
cd $SLURM_SUBMIT_DIR
RESULTSDIR=$SLURM_SUBMIT_DIR/Results_Each_Iter/Iteration$IterationCount
CGEDIR=$SLURM_SUBMIT_DIR/modelcge


if [ ! -d "results/Iteration$IterationCount" ]
then
 echo "Creating directory results/Iteration$IterationCount"
 mkdir -p results/Iteration$IterationCount
fi

CONRESDIR=$SLURM_SUBMIT_DIR/results/Iteration$IterationCount

matlab -nodisplay -nosplash -r "Trial_Coupling($IterationCount);exit" > output200$IterationCount.txt
rc=$?

echo "Convergence job complete with return code of $rc"

cp CheckConvgJob.csv $RESULTSDIR
mv output200$IterationCount.txt $RESULTSDIR
mv mat_output*.txt $RESULTSDIR


#Name the CGE directory to copy CGEtoM.gdx and CGEtoM2.gdx files
cd $CGEDIR
mv CGEtoM.gdx $CONRESDIR
mv CGEtoM2.gdx $CONRESDIR
mv MtoCGE.gdx $CONRESDIR
mv MtoCGE2.gdx $CONRESDIR

cd $SLURM_SUBMIT_DIR

cp $SLURM_SUBMIT_DIR/DemandChangeData_zones.mat $CONRESDIR
cp $SLURM_SUBMIT_DIR/DemandChangeData_econr.mat $CONRESDIR

if [ $IterationCount -eq 1 ]
then
 echo "Copying Unitcommitment_5000.gms"
 cp $SLURM_SUBMIT_DIR/data/programs_uc/UnitCommitment_5000.gms $SLURM_SUBMIT_DIR/Gams_data/UnitCommitment.gms
fi
if [ $IterationCount -eq 2 ]
then
 echo "Copying Unitcommitment_4000.gms"
 cp $SLURM_SUBMIT_DIR/data/programs_uc/UnitCommitment_4000.gms $SLURM_SUBMIT_DIR/Gams_data/UnitCommitment.gms
fi
if [ $IterationCount -eq 3 ]
then
 echo "Copying Unitcommitment_3000.gms"
 cp $SLURM_SUBMIT_DIR/data/programs_uc/UnitCommitment_3000.gms $SLURM_SUBMIT_DIR/Gams_data/UnitCommitment.gms
fi
if [ $IterationCount -ge 4 ]
then
 echo "Copying Unitcommitment_2000.gms"
 cp $SLURM_SUBMIT_DIR/data/programs_uc/UnitCommitment_2000.gms $SLURM_SUBMIT_DIR/Gams_data/UnitCommitment.gms
fi



 
echo "Ending convergence run for iteration Count $IterationCount" 

exit $rc
