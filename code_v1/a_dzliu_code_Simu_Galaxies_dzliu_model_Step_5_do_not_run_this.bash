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
# Step 5
# 
echo "macro read a_dzliu_code_Simu_Galaxies_dzliu_model.sm generate_Galaxy_Fluxes" | sm | tee "log_Simu_Galaxies_dzliu_model_Step_5.txt"


sm << EOF
foreach var {Ks irac1 irac2 irac3 irac4 24 100 160 250 350 450 500 850 1100 1200 1250 1300 10cm 20cm} \\
{data result_simu_galaxies_flux_at_\$var.txt read f\$var 4.f}
print "result_simu_galaxies_flux_at_selected_bands.txt" \\
'%15.6e %15.6e %15.6e %15.6e %15.6e %15.6e %15.6e %15.6e %15.6e %15.6e %15.6e %15.6e %15.6e %15.6e %15.6e %15.6e %15.6e %15.6e %15.6e\n' \\
{fKs firac1 firac2 firac3 firac4 f24 f100 f160 f250 f350 f450 f500 f850 f1100 f1200 f1250 f1300 f10cm f20cm}
EOF


