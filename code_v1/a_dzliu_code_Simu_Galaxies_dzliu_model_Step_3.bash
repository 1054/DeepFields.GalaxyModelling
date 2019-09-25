#!/bin/bash
# 

#source ~/Cloud/Github/DeepFields.SuperDeblending/Softwares/SETUP


# 
# Check
# 
if [[ ! -f "result_simu_galaxies_list.txt" ]]; then
    echo "Error! \"result_simu_galaxies_list.txt\" was not found! Please first run \"a_dzliu_code_Simu_Galaxies_dzliu_model_Step_2.bash\"!"; exit
fi
if [[ ! -d "result_simu_galaxies" ]]; then
    echo "Error! \"result_simu_galaxies\" was not found! Please first run \"a_dzliu_code_Simu_Galaxies_dzliu_model_Step_[12].bash\"!"; exit
fi


# 
# Step 3
# 
echo "macro read a_dzliu_code_Simu_Galaxies_dzliu_model.sm generate_Galaxy_SEDs" | sm | tee "log_Simu_Galaxies_dzliu_model_Step_3.txt"


