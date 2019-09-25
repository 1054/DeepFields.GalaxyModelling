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
iSTOP=699262 # 100 # 699262
iBand="Ks"

for (( i=$iSTART; i<=$iSTOP; i++ )); do
    if [[ ! -f "result_simu_galaxies_flux_at_$iBand/$i.txt" ]]; then
        #echo "Current computer \"$(hostname)\""
        #echo "Current directory \"$RunFolder\""
        echo "Counting screen task "$(screen -ls | grep "fgmod_" | wc -l)
        while       [[   30  -le    $(screen -ls | grep "fgmod_" | wc -l) ]]; do
        sleep             3
        echo "Counting screen task "$(screen -ls | grep "fgmod_" | wc -l)
        done
        RunScreen="fgmod_$(date +%s)_$i"
        echo "Adding screen task \"$RunScreen\" with command \"echo macro read a_dzliu_code_Simu_Galaxies_dzliu_model.sm generate_Galaxy_Fluxes $i | sm\""
        screen -d -S "$RunScreen" -m bash -c "cd $RunFolder; echo macro read a_dzliu_code_Simu_Galaxies_dzliu_model.sm generate_Galaxy_Fluxes $i | sm"
        echo Sleeping 0.15s
        sleep 0.15
    fi
done


