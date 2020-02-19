#!/bin/bash

module unload root
module unload python
unset LD_LIBRARY_PATH
unset PYTHONPATH
conda activate fermi
