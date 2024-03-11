#!/bin/bash

# Script to show how to launch the program (just a template)

file_gro="../input_files/WT_optimize_40_k/first_frame.gro"
excogito_dir="../input_files/WT_optimize_40_k"
code="WT"
num_propre="258"

PYTHON_DIR="../python_scripts"

if [ -d ../output_files/$code ]; then
        rm -r ../output_files/$code
fi
mkdir ../output_files/$code

python $PYTHON_DIR/SA.py -g $file_gro -d $excogito_dir -c $code -n $num_propre

python $PYTHON_DIR/post_processing.py -n $num_propre -d WT
