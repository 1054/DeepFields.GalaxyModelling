#!/bin/bash
# 

source ~/Cloud/Github/DeepFields.SuperDeblending/Softwares/SETUP


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
# Check computer
# 
if [[ $(hostname) != *"aida42198"* ]]; then 
    echo "Error! Unable to work on current computer \"$(hostname)\"!"
fi


# 
# Step 5
# 
RunFolder=$(readlink -f $(dirname ${BASH_SOURCE[0]}))
iSTART=1
iSTOP=699262
iBand="Ks" # 870

for (( i=$iSTART; i<=$iSTOP; i++ )); do
    echo $i
    for band in ${iBand[@]}; do
        if [[ ! -f "result_simu_galaxies_flux_at_$band/$i.txt" ]]; then
            echo "Errro! \"result_simu_galaxies_flux_at_$band/$i.txt\" was not found!"
            exit 1
        fi
        if [[ $i == 1 ]]; then
            cat "result_simu_galaxies_flux_at_$band/$i.txt" > "result_simu_galaxies_flux_at_$band.txt"
        else
            cat "result_simu_galaxies_flux_at_$band/$i.txt" >> "result_simu_galaxies_flux_at_$band.txt"
        fi
    done
done


# head -n 32 result_simu_galaxies_flux_at_1250.txt


