#!/bin/bash
# 

#source ~/Cloud/Github/DeepFields.SuperDeblending/Softwares/SETUP


# 
# Check
# 
if [[ ! -d "result_datatable_per_redshift_bin" ]]; then
    echo "Error! \"result_datatable_per_redshift_bin\" was not found! Please first run \"a_dzliu_code_Plot_SMF_dzliu_model.bash\"!"; exit
fi


# 
# Step 1
# 
echo "macro read a_dzliu_code_Simu_Galaxies_dzliu_model.sm simu_Galaxies" | sm | tee "log_Simu_Galaxies_dzliu_model_Step_1.txt"


