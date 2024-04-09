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

for i in {0..2}; do

	python $PYTHON_DIR/SA.py -g $file_gro -d $excogito_dir -c ${code}_$i -n $num_propre

	if [ ! $i -eq 0 ]; then

		count=0
		for file in $( ls ../output_files/${code}_$i ); do
			if [[ $file == *".log" ]]; then
				count=$(($count + 1))
			fi
		done
		count=$(($count - 1))

		num_files=$(eval echo {0..$count})
		for j in $num_files; do
			new_j=$( echo "$j + $i * ($count + 1)" | bc )
			mv ../output_files/${code}_${i}/SA_${j}.log ../output_files/${code}_${i}/SA_${new_j}.log
			mv ../output_files/${code}_${i}/SA_data_${j}.txt ../output_files/${code}_${i}/SA_data_${new_j}.txt
		done

	fi

	cp ../output_files/${code}_$i/* ../output_files/$code/
	rm -r ../output_files/${code}_$i

done

python $PYTHON_DIR/post_processing.py -n $num_propre -d WT
