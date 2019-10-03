#!/bin/bash
# 

set -e


script_dir=$(dirname "${BASH_SOURCE[0]}")
echo "script_dir = $script_dir"


echo "mkdir"
if [[ ! -d "$script_dir/python/lib/crab" ]]; then
    mkdir -p "$script_dir/python/lib/crab"
fi
if [[ ! -d "$script_dir/python/lib/a3cosmos-gas-evolution" ]]; then
    mkdir -p "$script_dir/python/lib/a3cosmos-gas-evolution"
fi


echo "cp"
cp -r ~/Cloud/Github/Crab.Toolkit.Python/lib/crab/crabgalaxy "$script_dir/python/lib/crab/"
cp -r ~/Cloud/Github/Crab.Toolkit.Python/lib/crab/crabtable "$script_dir/python/lib/crab/"

#cp ~/Cloud/GitLab/AlmaCosmos/Plot/Common_Python_Code/calc_galaxy_stellar_mass_function.py "$script_dir/python/lib/a3cosmos-gas-evolution/"
#cp ~/Cloud/GitLab/AlmaCosmos/Plot/Common_Python_Code/calc_galaxy_main_sequence.py "$script_dir/python/lib/a3cosmos-gas-evolution/"
#cp ~/Cloud/GitLab/AlmaCosmos/Plot/Common_Python_Code/calc_cosmic_star_formation_rate_density.py "$script_dir/python/lib/a3cosmos-gas-evolution/"
#cp ~/Cloud/GitLab/AlmaCosmos/Plot/Common_Python_Code/apply_cosmology.py "$script_dir/python/lib/a3cosmos-gas-evolution/"
#cp ~/Cloud/GitLab/AlmaCosmos/Plot/Common_Python_Code/setup_matplotlib.py "$script_dir/python/lib/a3cosmos-gas-evolution/"


echo "Done!"

