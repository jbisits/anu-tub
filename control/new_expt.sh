# copying the config.yaml, MOM_override from the control version of the same experiment.
# E.g. zstar-dec-diff would create a new experiment directory based off zstar. This is to
# avoid manually doing this when there are many new experiments to create.
new_perturbation_expt() {
    local new_expt="$1"
    local override="$2"

    # 1. Error if no directory is passed
    if [ -z "$new_expt" ]; then
        echo "Error: Please specify a destination directory." >&2
        return 1
    fi
    # 2. If the directory doesn't exist, create it (including parent folders)
    if [ -d "$new_expt" ]; then
        echo "Error: Directory '$new_expt' already exists. Aborting to protect data." >&2
        return 1
    else
        echo "Experiment directory '$new_expt' does not exist. Creating it now..."
        mkdir -p "$new_expt"
    fi

    # 3. Extract the first 3 letters of the new experiment and other directory info
    local control_dir=$(pwd)
    local common_dir="${control_dir}/common"
    local expt_dir_name=$(basename "$new_expt")
    local copy_expt="${expt_dir_name:0:3}"
    local control_expt=""
    local output_dir=""

    # match the new experiment to a control
    if [[ "$copy_expt" == "zst" ]]; then
	control_expt="zstar"
	output_dir="zstar"
    elif [[ "$copy_expt" == "hyc" ]]; then
	control_expt="hycom1"
	output_dir="hycom"
    elif [[ "$copy_expt" == "ada" ]]; then
	control_expt="adapt"
	output_dir="AG"
    elif [[ "$copy_expt" == "alt" ]]; then
	control_expt="alt-hycom"
	output_dir="alt-hycom"
    fi

    echo "Copying config.yaml and MOM_override from $control_expt to $new_expt"
    local CFG="config.yaml"
    local MO="MOM_override"
    cp "${control_dir}/${control_expt}/${CFG}" "$new_expt"
    cp "${control_dir}/${control_expt}/${MO}" "$new_expt"

    echo "Copying sync script"
    local SYNC="sync_output_to_gdata.sh"
    cp "${control_dir}/${control_expt}/${SYNC}" "$new_expt"

    echo "Copying set_default from $common_dir to $new_expt"
    local SD="set_default.sh"
    cp "${common_dir}/${SD}" "$new_expt"
 
    echo "Set up default inputs and update $control_expt to $new_expt where appropriate and add override if present"
    pushd "$new_expt" > /dev/null
    ./set_default.sh
    sed -i "s|${output_dir}|${new_expt}|g" "$SYNC"
    sed -i "s|one-deg-tub-${output_dir,,}|${new_expt}|g" "$CFG"
    if [ -z "$override" ]; then
        echo "Nothing appended to MOM_override so expt is unchanged."
    else
        echo "appending ${override} to $MO"
        sed -i "\$a #override ${override}" "$MO"
    fi
    payu_setup_new_expt
    popd > /dev/null
}
# After a new experiment directory has been, set up a payu experiment and run it.
payu_setup_new_expt(){
    module use /g/data/vk83/modules
    module load payu
    payu setup
    local work_dir_target=$(readlink "work")
    local exptname=$(basename "$work_dir_target")
    copy_restart_to_archive "$exptname"
    git init
    payu sweep
    # payu setup
    # payu sweep
    # payu run
}

# Copy the most restart file based on the vertical coordinate to the `archive`.
copy_restart_to_archive() {
    local exptname=$(basename "$1")
    local archivepath="/scratch/e14/jb2381/mom6/archive/$exptname"
    local restart_path="/g/data/e14/jb2381/one-degree-anu-tub/restarts"
    local vc="${exptname:0:3}"
    local verticalcoord=""
    local restart_directory=""

    if [[ "$vc" == "zst" ]]; then
        verticalcoord="zstar"
	restart_directory="restart051"
    elif [[ "$vc" == "ada" ]]; then
        verticalcoord="AG"
	restart_directory="restart051"
    elif [[ "$vc" == "hyc" ]]; then
        verticalcoord="hycom"
	restart_directory="restart051"
    elif [[ "$vc" == "alt" ]]; then
        verticalcoord="alt-hycom"
	restart_directory="restart050"
    fi

    echo "Copying ${restart_path}/${verticalcoord}/${restart_directory} to $archivepath"
    cp -r "${restart_path}/${verticalcoord}/${restart_directory}" "$archivepath"
}
