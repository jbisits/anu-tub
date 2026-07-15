#!/bin/bash

# create a new experiment directory by copying relevant files across

# 1. Source (import) the new_expt file into this script environment
source "new_expt.sh"

# 2. Call the function, passing the directory argument along
new_perturbation_expt "$1" "$2"
