#!/bin/bash

CONDA_DIR=/a/data/tehanu/$USER/miniconda3

dry_run=0

print_usage() {
    printf "Usage: setup_pipeline [-d] [-p <CONDA_PATH>]
            -d  Do a dry run (do not install)
            -p  Install miniconda to the specified path
            \n"
}

while getopts 'dp:' flag; do
    case "${flag}" in
        d) dry_run=1 ;;
        p) CONDA_DIR="${OPTARG}" ;;
        *) print_usage
           exit 1 ;;
    esac
done

if [[ $dry_run == 1 ]]
then
    echo "$FERMIPIPE/miniconda3.sh -b -p $CONDA_DIR"
    exit 0
fi

# Set up the conda environment manager
bash $FERMIPIPE/Miniconda3-latest-Linux-x86_64.sh -b -p $CONDA_DIR
source $CONDA_DIR/bin/activate

# Save the initialization text for other users to copy into their .bashrc
conda init -v -d | grep "^+" | grep -v "^+++" | cut -c2- > $FERMIPIPE/conda_init.txt

# Initialize conda
conda init
source ~/.bashrc
conda config --set auto_activate_base false

# Download and install the Fermitools and Fermipy
conda env create -f $FERMIPIPE/environment.yml

exit 0
