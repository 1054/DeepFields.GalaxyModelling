#!/bin/bash
# 

#source ~/Cloud/Github/DeepFields.SuperDeblending/Softwares/SETUP


# 
# Check
# 
if [[ ! -f "log_Simu_Galaxies_dzliu_model_Step_1.txt" ]]; then
    echo "Error! \"log_Simu_Galaxies_dzliu_model_Step_1.txt\" was not found! Please first run \"a_dzliu_code_Simu_Galaxies_dzliu_model_Step_1.bash\"!"; exit
fi
if [[ ! -d "result_simu_galaxies" ]]; then
    echo "Error! \"result_simu_galaxies\" was not found! Please first run \"a_dzliu_code_Simu_Galaxies_dzliu_model_Step_1.bash\"!"; exit
fi


# 
# Step 2
# 
cat "log_Simu_Galaxies_dzliu_model_Step_1.txt" | grep "# CrabPhotMonteCarlo: writing to" | sed -e 's/# CrabPhotMonteCarlo: writing to //g' > "result_simu_galaxies_list.txt"
date +"# %Y-%m-%d %Hh%Mm%Ss %Z" > "log_Simu_Galaxies_dzliu_model_Step_2.txt"
echo "Output to \"result_simu_galaxies_list.txt\"!" 
echo "Output to \"result_simu_galaxies_list.txt\"!" >> "log_Simu_Galaxies_dzliu_model_Step_2.txt"


