#!/bin/bash
# set -x

#############################################################################################
# Script to run PSM-DREM coupled model. This script allows the user to specific which
# scenerios, global climate model, demand response, WBM year, and reitrement case to run on
# ACI.
#
# There is currently only a single region specification. If necessary, I can add this into a
# future version
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
echo "===Run PSM-DREM Coupled Model under a Single Specification===" 
echo

#=========================================================================================#
#     USE DEFINED VARIABLES
#=========================================================================================#

echo '----User defined variables----'
echo

# define the number of regions to use # currently only 14 region files are setup
REGION=14                    
# define global climate model to use
GCM=GFDL-CM3
# define water scenario years      
WS=64        
# define retirement scenario      
RS=27     
# define which scenario to runtime 
RUN_SCEN=3    
# define whether the demand response is on/off 
RUN_DR=1        

echo "Economic Region: ${REGION}"
echo "Scenario Cases: ${RUN_SCEN}"
echo "Demand Resonse Cases: ${RUN_DR}"
echo "Global Climate Models: ${GCM}"
echo "Water Scenarios: ${WS}"
echo "Retirement Scenarios: ${RS}"

#=========================================================================================#
#     create economic region name for copying files
#=========================================================================================#

# define economic region name
ECONR=${REGION}econr

echo "Economic Region Name: ${ECONR}"
echo

#=========================================================================================#
#     define directory locations 
#=========================================================================================#

echo "----Define Model and Result Directories in codedir----"
echo

DIR_PSM=Gams_data          # psm model directory
DIR_REM=modelcge           # rem model directory
DIR_RES1=results           # main coupled model results
DIR_RES2=Results_Each_Iter # detailed psm results

#=========================================================================================#
#     remove old files from perviious run
#=========================================================================================#

echo "----Remove Old Result and Output Files (if they exist)----"
echo

# Delete all results subfolders to start clean
#     MAKE SURE PREVIOUS CASE HAS BEEN COPIED!!!!!!!!
rm -rf ${DIR_RES1}/Iteration*
rm -rf ${DIR_RES2}/Iteration*
rm -rf ${DIR_REM}/listings/*.*
rm -rf ${DIR_REM}/logs/*.*
rm -rf ${DIR_REM}/output/csv/*.*
rm -rf ${DIR_REM}/output/gdx/*.*
rm -rf ${DIR_REM}/output/coupling/drm/*.*

rm CheckConvgJob*
rm outf*
rm -rf DemandChangeData*
# rm *o.*

#=========================================================================================#
#     define directory names for the multiple 
#=========================================================================================#

echo "----Define Coupled Model Specification----"
echo

if [ ${RUN_SCEN} -eq 1 ] && [ ${RUN_DR} -eq 0 ]                           # scenario: WBM only 
then

	# define variable list for qsub
	varlist="ECONR=${ECONR}, SHOCK=1, DEMANDRESPONSE=0, MODEL=${GCM}, YEAR=${WS}, RETIRECASE=27"        
						
	# define job name
	DIR_SCEN=${GCM}_w${WS}
	
	# magic qsub call line
	qsub -N ${DIR_SCEN} -v "${varlist}" Run_ManagePSM.pbs

	echo						
	echo "qsub Run_ManagePSM.pbs for ${DIR_SCEN}: econonomic region set to ${ECONR}; GCM set to ${GCM}; scenario set to 1; demand response set to 0; water shock year set to ${WS}; ret shock set to 27 (NOT USED)"
	echo			
		
elif [ ${RUN_SCEN} -eq 2 ] && [ ${RUN_DR} -eq 0 ]                         # scenario: retirements only 
then

	# define variable list for qsub        
	varlist="ECONR=${ECONR}, SHOCK=2, DEMANDRESPONSE=0, MODEL=GFDL-CM3, YEAR=64, RETIRECASE=${RS}"                              

	# define job name
	DIR_SCEN=r${RS}
	
	# magic qsub call line              
	qsub -N ${DIR_SCEN} -v "${varlist}" Run_ManagePSM.pbs

	echo						
	echo "qsub Run_ManagePSM.pbs for ${DIR_SCEN}: econonomic region set to ${ECONR}; GCM set to ${GCM} (NOT USED); scenario set to 2; demand response set to 0; water shock year set to 64 (NOT USED); ret shock set to ${RS}"
	echo	
	
elif [ ${RUN_SCEN} -eq 3 ] && [ ${RUN_DR} -eq 0 ]                         # scenario: WBM and retirements 
then

	# define variable list for qsub        
	varlist="ECONR=${ECONR}, SHOCK=3, DEMANDRESPONSE=0, MODEL=${GCM}, YEAR=${WS}, RETIRECASE=${RS}"
			
	# create a temp directory in scratch to run the code
	DIR_SCEN=${GCM}_w${WS}r${RS}
	
	# magic qsub call line
	qsub -N ${DIR_SCEN} -v "${varlist}" Run_ManagePSM.pbs
	
	echo
	echo "qsub Run_ManagePSM.pbs for ${DIR_SCEN}: econonomic region set to ${ECONR}; GCM set to ${GCM}; scenario set to 3; demand response set to 0; water shock year set to ${WS}; ret shock set to ${RS}"
	echo
			
elif [ ${RUN_SCEN} -eq 1 ] && [ ${RUN_DR} -eq 1 ]                         # scenario: Demand response with WBM
then

	# define variable list for qsub        
	varlist="ECONR=${ECONR}, SHOCK=1, DEMANDRESPONSE=1, MODEL=${GCM}, YEAR=${WS}, RETIRECASE=27"

	# create a temp directory in scratch to run the code
	DIR_SCEN=d${GCM}_w${WS}
	
	# magic qsub call line
	qsub -N ${DIR_SCEN} -v "${varlist}" Run_ManagePSM.pbs

	echo						
	echo "qsub Run_ManagePSM.pbs for ${DIR_SCEN}: econonomic region set to ${ECONR}; GCM set to ${GCM}; scenario set to 1; demand response set to 1; water shock year set to ${WS}; ret shock set to 27 (NOT USED)"
	echo					

elif [ ${RUN_SCEN} -eq 2 ] && [ ${RUN_DR} -eq 1 ]                         # scenario: Demand response with retirements 
then

	# define variable list for qsub        
	varlist="ECONR=${ECONR}, SHOCK=2, DEMANDRESPONSE=1, MODEL=${GCM}, YEAR=64, RETIRECASE=${RS}"

	# define job name
	DIR_SCEN=d${GCM}_r${RS}

	# magic qsub call line
	qsub -N ${DIR_SCEN} -v "${varlist}" Run_ManagePSM.pbs
	
	echo                    
	echo "qsub Run_ManagePSM.pbs for ${DIR_SCEN}: econonomic region set to ${ECONR}; GCM set to ${GCM}; scenario set to 2; demand response set to 1; water shock year set to 64 (NOT USED); ret shock set to ${RS}"
	echo	
			
elif [ ${RUN_SCEN} -eq 3 ] && [ ${RUN_DR} -eq 1 ]                         # scenario: Demand response with WBM and retirements 
then

	# define variable list for qsub        
	varlist="ECONR=${ECONR}, SHOCK=3, DEMANDRESPONSE=1, MODEL=${GCM}, YEAR=${WS}, RETIRECASE=${RS}"

	# define job name
	DIR_SCEN=d${GCM}w${WS}r${RS}
	
	# magic qsub call line
	qsub -N ${DIR_SCEN} -v "${varlist}" Run_ManagePSM.pbs

	echo						
	echo "qsub Run_ManagePSM.pbs for ${DIR_SCEN}: econonomic region set to ${ECONR}; GCM set to ${GCM}; scenario set to 3; demand response set to 1; water shock year set to ${WS}; ret shock set to ${RS}"
	echo
	
fi # of if statement
