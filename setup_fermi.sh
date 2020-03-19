#!/bin/bash

module unload root
module unload python
unset LD_LIBRARY_PATH
unset PYTHONPATH
export PYTHONUSERSITE=True
conda activate fermi
ulimit -n 2048
