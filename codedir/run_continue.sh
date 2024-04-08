#!/bin/bash
#set -x

#############################################################################################
# Script to continue running the PSM-DREM coupled model. The specification is based on the 
# run_single_start.sh. 
#
# No input from the user is required. THE USER MUST CHANGE Run_ManagePSM.sbatch WITH THE NEW
# ITERATION START AND STOP VALUES!!!!!!!!!!!!!
#
# Created by Joseph Perla (jmp629) and based on script from Vijay Kumar (vuk19). 
# Created: 2021-06-02
# Updated: 2021-06-02
#
# ## informtion on variable
# REGION        : 14
# MODEL         : GFDL-CM3  ; CCSM4
# YEAR          : 50, 64, 89; 72, 94, 93 
# RETIRECASE    : 27, 25, 28
# SHOCK         : == 0 baseline, == 1 Water (WBM), == 2 Retirement, == 3 both Water and Retire
# DEMANDRESPONSE: == 0 disable demand response, == 1 enable  demand response
###############################################################################################

echo
echo "===Continue PSM-DREM Coupled Model Run===" 
echo

#=========================================================================================#
#     USE DEFINED VARIABLES
#=========================================================================================#

echo '----User defined variables from run_single_start.sh----'
echo

# define the number of regions to use # currently only 14 region files are setup
REGION=$(grep 'REGION=' ./run_start.sh | cut -d '=' -f2 | tr -d '[:space:]')                 
# define global climate model to use
GCM=$(grep 'GCM=' ./run_start.sh | cut -d '=' -f2 | tr -d '[:space:]')
# define water scenario years      
WS=$(grep 'WS=' ./run_start.sh | cut -d '=' -f2 | tr -d '[:space:]')       
# define retirement scenario      
RS=$(grep 'RS=' ./run_start.sh | cut -d '=' -f2 | tr -d '[:space:]')    
# define which scenario to runtime 
RUN_SCEN=$(grep 'RUN_SCEN=' ./run_start.sh | cut -d '=' -f2 | tr -d '[:space:]')   
# define whether the demand response is on/off 
RUN_DR=$(grep 'RUN_DR=' ./run_start.sh | cut -d '=' -f2 | tr -d '[:space:]')       

echo "Economic Region: ${REGION[@]}"
echo "Scenario Cases: ${RUN_SCEN[@]}"
echo "Demand Resonse Cases: ${RUN_DR[@]}"
echo "Global Climate Models: ${GCM}"
echo "Water Scenarios: ${WS[@]}"
echo "Retirement Scenarios: ${RS[@]}"

#=========================================================================================#
#     create economic region name for copying files
#=========================================================================================#

# define economic region name
ECONR="${REGION}econr"

echo "Economic Region Name: ${ECONR}"
echo 

#=========================================================================================#
#     define variable list for sbatch 
#=========================================================================================#

echo "----Define Coupled Model Specification----"
echo

if [ ${RUN_SCEN} -eq 1 ] && [ ${RUN_DR} -eq 0 ]                           # scenario: WBM only 
then

	# define variable list for sbatch
	varlist="ECONR=${ECONR}, SHOCK=1, DEMANDRESPONSE=0, MODEL=${GCM}, YEAR=${WS}, RETIRECASE=27"        
						
	# define job name
	DIR_SCEN=${GCM}_w${WS}
	
	# magic sbatch call line
	sbatch -J ${DIR_SCEN} -v "${varlist}" Run_ManagePSM.sbatch  ## we need to test this

	echo						
	echo "sbatch Run_ManagePSM.sbatch for ${DIR_SCEN}: econonomic region set to ${ECONR}; GCM set to ${GCM}; scenario set to 1; demand response set to 0; water shock year set to ${WS}; ret shock set to 27 (NOT USED)"
	echo			
		
elif [ ${RUN_SCEN} -eq 2 ] && [ ${RUN_DR} -eq 0 ]                         # scenario: retirements only 
then

	# define variable list for sbatch        
	varlist="ECONR=${ECONR}, SHOCK=2, DEMANDRESPONSE=0, MODEL=GFDL-CM3, YEAR=64, RETIRECASE=${RS}"                              

	# define job name
	DIR_SCEN=r${RS}
	
	# magic sbatch call line              
	sbatch -J ${DIR_SCEN} -v "${varlist}" Run_ManagePSM.sbatch

	echo						
	echo "sbatch Run_ManagePSM.sbatch for ${DIR_SCEN}: econonomic region set to ${ECONR}; GCM set to ${GCM} (NOT USED); scenario set to 2; demand response set to 0; water shock year set to 64 (NOT USED); ret shock set to ${RS}"
	echo	
	
elif [ ${RUN_SCEN} -eq 3 ] && [ ${RUN_DR} -eq 0 ]                         # scenario: WBM and retirements 
then

	# define variable list for sbatch        
	varlist="ECONR=${ECONR}, SHOCK=3, DEMANDRESPONSE=0, MODEL=${GCM}, YEAR=${WS}, RETIRECASE=${RS}"
			
	# create a temp directory in scratch to run the code
	DIR_SCEN=${GCM}_w${WS}r${RS}
	
	# magic sbatch call line
	sbatch -J ${DIR_SCEN} -v "${varlist}" Run_ManagePSM.sbatch
	
	echo
	echo "sbatch Run_ManagePSM.sbatch for ${DIR_SCEN}: econonomic region set to ${ECONR}; GCM set to ${GCM}; scenario set to 3; demand response set to 0; water shock year set to ${WS}; ret shock set to ${RS}"
	echo
			
elif [ ${RUN_SCEN} -eq 1 ] && [ ${RUN_DR} -eq 1 ]                         # scenario: Demand response with WBM
then

	# define variable list for sbatch        
	varlist="ECONR=${ECONR}, SHOCK=1, DEMANDRESPONSE=1, MODEL=${GCM}, YEAR=${WS}, RETIRECASE=27"

	# create a temp directory in scratch to run the code
	DIR_SCEN=d${GCM}_w${WS}
	
	# magic sbatch call line
	sbatch -J ${DIR_SCEN} -v "${varlist}" Run_ManagePSM.sbatch

	echo						
	echo "sbatch Run_ManagePSM.sbatch for ${DIR_SCEN}: econonomic region set to ${ECONR}; GCM set to ${GCM}; scenario set to 1; demand response set to 1; water shock year set to ${WS}; ret shock set to 27 (NOT USED)"
	echo					

elif [ ${RUN_SCEN} -eq 2 ] && [ ${RUN_DR} -eq 1 ]                         # scenario: Demand response with retirements 
then

	# define variable list for sbatch        
	varlist="ECONR=${ECONR}, SHOCK=2, DEMANDRESPONSE=1, MODEL=${GCM}, YEAR=64, RETIRECASE=${RS}"

	# define job name
	DIR_SCEN=d${GCM}_r${RS}

	# magic sbatch call line
	sbatch -J ${DIR_SCEN} -v "${varlist}" Run_ManagePSM.sbatch
	
	echo                    
	echo "sbatch Run_ManagePSM.sbatch for ${DIR_SCEN}: econonomic region set to ${ECONR}; GCM set to ${GCM}; scenario set to 2; demand response set to 1; water shock year set to 64 (NOT USED); ret shock set to ${RS}"
	echo	
			
elif [ ${RUN_SCEN} -eq 3 ] && [ ${RUN_DR} -eq 1 ]                         # scenario: Demand response with WBM and retirements 
then

	# define variable list for sbatch        
	varlist="ECONR=${ECONR}, SHOCK=3, DEMANDRESPONSE=1, MODEL=${GCM}, YEAR=${WS}, RETIRECASE=${RS}"

	# define job name
	DIR_SCEN=d${GCM}w${WS}r${RS}
	
	# magic sbatch call line
	sbatch -J ${DIR_SCEN} -v "${varlist}" Run_ManagePSM.sbatch

	echo						
	echo "sbatch Run_ManagePSM.sbatch for ${DIR_SCEN}: econonomic region set to ${ECONR}; GCM set to ${GCM}; scenario set to 3; demand response set to 1; water shock year set to ${WS}; ret shock set to ${RS}"
	echo
	
fi # of if statement
