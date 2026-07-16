#!/bin/bash

# create a new experiment directory by copying relevant files across

# Source (import) the new_expt file into this script environment
source "new_expt.sh"

# Control experiments, PERTURBATION is appended to each base experiment name
expts=("zstar" "hycom1" "alt-hycom" "adapt")
PERTURBATION="-weaker-KD"
# OVERRIDE is appended to the MOM_override file in the new experiment
OVERRIDE="KD = 1E-6"
# Call the function, passing the new expt control directory and OVERRIDE
# new_perturbation_expt "$1" "$PERTURBATION"
# for e in "${expts[@]}"; do
#     p_expt="${e}${PERTURBATION}"
#     echo "Creating a new control directory: $p_expt"
#     new_perturbation_expt "$p_expt" "$OVERRIDE"
# done
new_perturbation_expt "adapt-test" "$OVERRIDE"
