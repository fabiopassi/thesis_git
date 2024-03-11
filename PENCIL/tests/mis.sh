#!/bin/bash

# Script to show how to launch the program (just a template)

file_gro="../input_files/misfolded_optimize_40_k/first_frame.gro"
excogito_dir="../input_files/misfolded_optimize_40_k"
code="misfolded"
num_propre="175"

PYTHON_DIR="../python_scripts"

if [ -d ../output_files/$code ]; then
        rm -r ../output_files/$code
fi
mkdir ../output_files/$code

python $PYTHON_DIR/SA.py -g $file_gro -d $excogito_dir -c $code -n $num_propre

python $PYTHON_DIR/post_processing.py -d misfolded -n $num_propre
