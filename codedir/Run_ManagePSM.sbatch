#!/bin/bash

### #SBATCH --job-name=TEST1-run1
#SBATCH --output=scen.%j.out
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=13
#SBATCH --mem-per-cpu=13GB
#SBATCH --time=35:00:00
#SBATCH --constraint=sc
#SBATCH --account=pches
#SBATCH --partition=sla-prio
### #SBATCH --partition=open
#SBATCH --mail-type=ALL,TIME_LIMIT_90
#SBATCH --mail-user=jmp629@psu.edu


#=============================================================================#
#     user Defined Parameters for the FORTRAN code
#=============================================================================#

export I_MANAGE_PSM_ITER_START=1                        # set start iteration   - if continuing a iteration, change start and max values
export I_MANAGE_PSM_ITER_MAX=35                         # set maximum iteration - Note: iter_max > iter_start
export I_MANAGE_PSM_GAMS_THREADS=1                      # set the number threads: 1 to 8 -- CPLEX options threads must match -- not working?
export I_MANAGE_PSM_DEBUG_LVL=3                         # set debut level: == 0 minimum output, ==1, ==2, up to 3 == verbose
export I_MANAGE_PSM_JOBDIR="/scratch/jmp629/run-"        # set scratch directory must be quotations: Example - "/scratch/jmp629/run-"

#=========================================================================================#
#     set the following parameters to completely specify a case to run
#=========================================================================================#

# MODEL        : GFDL-CM3  ; CCSM4
# YEAR          : 50, 64, 89; 72, 94, 93 
# RETIRECASE    : 27, 25, 28
# SHOCK         : == 0 baseline, == 1 Water (WBM), == 2 Retirement, == 3 both Water and Retire
# DEMANDRESPONSE: == 0 disable demand response, == 1 enable  demand response

# MODEL=GFDL-CM3
# YEAR=50
# RETIRECASE=27

# SHOCK=1
# DEMANDRESPONSE=0

#=========================================================================================#
#     define directory locations 
#=========================================================================================#

DIR_PSM=Gams_data          # psm model directory
DIR_REM=modelcge           # rem model directory
DIR_RES1=results           # main coupled model results
DIR_RES2=Results_Each_Iter # detailed psm results

#=========================================================================================#
#     get started
#=========================================================================================#

echo "Job started on `hostname` at `date`"

## Load intel and intel MPI
module purge
module load intel
module load impi

## Go to the correct place

cd $SLURM_SUBMIT_DIR

#=============================================================================#
#     copy PSM model econmic region files
#=============================================================================#

# define the non-electricity productive multipliers (Sullivan multipliers)
cp data/files_economic_regions/${ECONR}_SullivanParameters.mat  SullivanParameters.mat

# define generator data file for Parallizing PSM
cp data/files_economic_regions/${ECONR}_Gendata2010.mat         ${DIR_PSM}/Gendata2010.mat

# define economic region crosswalk for PSM data import
cp data/files_economic_regions/${ECONR}_trialimport_econr.gms   ${DIR_PSM}/trialimport_econr.gms 

#=============================================================================#
#     copy PSM model starting files
#=============================================================================#

if [ ${I_MANAGE_PSM_ITER_START} -eq 1 ]
then

	# start iterations with a demand multiplier of zero
	cp data/files_economic_regions/${ECONR}_Input_data_DC_zeros.csv ${DIR_PSM}/Input_data_DC.csv
	cp data/files_economic_regions/${ECONR}_Input_data_UC_zeros.csv ${DIR_PSM}/Input_data_UC.csv

	# start iterations with longer runtime for UC model to find a good solution
	cp data/programs_uc/UnitCommitment_7000.gms Gams_data/UnitCommitment.gms

fi

#=============================================================================#
#     copy PSM scenario data files for the specific shock case requested
#=============================================================================#

# define retirement scenario file
cp data/scenarios_RET/gretire${RETIRECASE}.gms                                      ${DIR_PSM}/trialimport_gretire.gms

# define WBM gcm and year scenario file
cp data/scenarios_WBM/Year${YEAR}/PowerOutagesYear20New${YEAR}${MODEL}_rcp85.mat    ${DIR_PSM}/PowerOutages.mat

# define demand response gcm file
cp data/scenarios_DRM/load_scale_factors_${MODEL}_rcp85_294_busses_x_8760_hours.mat ${DIR_PSM}/load_scale_factors.mat

#=============================================================================#
#     copy TrialImportData with the specific scenario requested switch
#=============================================================================#

if [[ ${SHOCK} -eq 0 ]]    # baseline scenario 
then
cp data/files_trialimport_scenarios/TrialImportData_ShockNONE.gms   ${DIR_PSM}/TrialImportData.gms	
elif [[ ${SHOCK} -eq 1 ]]  # WBM scearnio
then
cp data/files_trialimport_scenarios/TrialImportData_ShockWBM.gms    ${DIR_PSM}/TrialImportData.gms	
elif [[ ${SHOCK} -eq 2 ]]  # retiement scenario
then
cp data/files_trialimport_scenarios/TrialImportData_ShockRET.gms    ${DIR_PSM}/TrialImportData.gms	
else
cp data/files_trialimport_scenarios/TrialImportData_ShockWBMRET.gms ${DIR_PSM}/TrialImportData.gms	
fi

#=============================================================================#
#     copy ParallelizingTrial and Trial_Coupling with or without demand response based
#     on the on/off switch
#=============================================================================#

if [ ${DEMANDRESPONSE} -eq 0 ]
then
cp data/files_parallelizing/ParallelizingTrial_REF.m  ParallelizingTrial.m 
cp data/files_coupling/Trial_Coupling_REF.m           Trial_Coupling.m
else
cp data/files_parallelizing/ParallelizingTrial_DR.m   ParallelizingTrial.m
cp data/files_coupling/Trial_Coupling_DR.m            Trial_Coupling.m
fi

#=============================================================================#
#     run the job
#=============================================================================#

echo "Job run starting at `hostname` at `date`" 
# cat /proc/cpuinfo
time srun ./Manage_PSM $STARTITER $ENDITER
echo "Job run ending at `hostname` at `date`" 

#=============================================================================#
#     finish up
#=============================================================================#

echo "Job Ended at `date`"
