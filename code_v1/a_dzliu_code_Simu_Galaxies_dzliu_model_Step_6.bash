#!/bin/bash
# 

#source ~/Cloud/Github/DeepFields.SuperDeblending/Softwares/SETUP


# 
# Check
# 
if [[ ! -f "result_simu_galaxies_radec.txt" ]]; then
    echo "Error! \"result_simu_galaxies_radec.txt\" was not found! Please first run \"a_dzliu_code_Simu_Galaxies_dzliu_model_Step_4.bash\"!"; exit
fi
if [[ ! -f "result_simu_galaxies_list.txt" ]]; then
    echo "Error! \"result_simu_galaxies_list.txt\" was not found! Please first run \"a_dzliu_code_Simu_Galaxies_dzliu_model_Step_[1234].bash\"!"; exit
fi
if [[ ! -d "result_simu_galaxies" ]]; then
    echo "Error! \"result_simu_galaxies\" was not found! Please first run \"a_dzliu_code_Simu_Galaxies_dzliu_model_Step_[1234].bash\"!"; exit
fi


# 
# Step 6 -- get z_Mstar_SFR
# 
echo "macro read a_dzliu_code_Simu_Galaxies_dzliu_model.sm print_Galaxy_z_Mstar_SFRs" | sm | tee "log_Simu_Galaxies_dzliu_model_Step_6.txt"


