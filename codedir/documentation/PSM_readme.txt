PSM Model Info:

The PSM Model is in 3 files Coal_Cali_Data, Gams_data and UC_Init_Solution. 
Coal_Cali_data is a file which has California plants coal generation. This is a redundant file and can be delted in next version.
Deletion of this file will require extensive code changes on the ACI system and hence has not been performed till now.

Gams_data is the main PSM folder and more information of each fileis given below.

UC_Init_Solution has a previous solution for each Unit Commitment model which can be used as a warm start for the next run.



Following are the files in Gams_data:

1. Bus_Final.gpr is a project file for the PSM.
2. cplex.opt is an option file to run the UnitCommitment and OPF models.
3. DCOPF.gms is the file with Optimal Power Flow model in GAMS.
4. Demandfor2010.csv is the demand data for year 2010.
5. Gendata2010.mat is the generation data and the renewable energy profile for each hour of the year.
6. input_data.csv is a file which has the demand multiplier for California and Rowecc obtained from IMPLAN.
7. PowerOutagesYear2090GFDL-CM3_rcp85.mat is the list of generators turned off when high temperatures are present. This file will change with different WBM outputs.
8. ReshapeMatrix is an accessory function to reshape the values obtained from GAMS. Can be redundant as we can use the 'full' option in rgdx.
9. TrialImportData.gms is a folder which has all the sets and other parameter definitions required to run UnitCommitment and DCOPF model.s
10. UnitCommitment.gms is the Unit commitment model of the system.
11. WECC_line_data2.txt is the line data for the WECC system. 


ParallelizingTrial.m has to be changed in the main folder( original working directory) and not inside Gams_data.   