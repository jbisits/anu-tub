#!/bin/bash

# create a new experiment directory by copying relevant files across

# Source (import) the new_expt file into this script environment
source "new_expt.sh"

# Control experiments, PERTURBATION is appended to each base experiment name
# OVERRIDE is appended to the MOM_override file in the new PETURBATION control directory
expts=("zstar" "hycom1" "alt-hycom" "adapt")
PERTURBATION="-weaker-KD"
OVERRIDE="KD = 1E-6"
# now create new control directories, copy restart files, setup payu and run!
for e in "${expts[@]}"; do
    p_expt="${e}${PERTURBATION}"
    echo "Creating a new control directory: $p_expt"
    new_perturbation_expt "$p_expt" "$OVERRIDE"
done
