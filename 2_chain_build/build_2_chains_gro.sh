#!/bin/bash

file_chain_1_gro="misfolded_chain_1.gro"
file_chain_2_gro="misfolded_chain_2.gro"

file_chain_1_pdb="misfolded_chain_1.pdb"
file_chain_2_pdb="misfolded_chain_2.pdb"

file_double_chain_gro="2_chains.gro"
file_double_chain_pdb="2_chains.pdb"

single_chains_dir="single_chains"
double_chain_dir="double_chain"

if [ -f $file_chain_1_pdb ]; then rm $file_chain_1_pdb; fi
if [ -f $file_chain_2_pdb ]; then rm $file_chain_2_pdb; fi
if [ -f $file_chain_2_gro ]; then rm $file_chain_2_gro; fi
if [ -f $file_double_chain_gro ]; then rm $file_double_chain_gro; fi
if [ -f $file_double_chain_pdb ]; then rm $file_double_chain_pdb; fi
if [ -f topol.top ]; then rm topol.top; fi
if [ -f posre.itp ]; then rm posre.itp; fi


./translate_gro.py

gmx editconf -f $file_chain_1_gro -o $file_chain_1_pdb
gmx editconf -f $file_chain_2_gro -o $file_chain_2_pdb

last_line=$( wc -l $file_chain_1_pdb | awk '{print $1}' )
sed -i "${last_line}d" $file_chain_1_pdb

sed -i '1,4d' $file_chain_2_pdb

cat $file_chain_1_pdb $file_chain_2_pdb > $file_double_chain_pdb
gmx pdb2gmx -f $file_double_chain_pdb -o $file_double_chain_gro -merge all -chainsep ter

if [ ! -d $single_chains_dir ]; then
    mkdir $single_chains_dir
else 
    rm -r ${single_chains_dir}/*
fi

mv $file_chain_1_pdb $file_chain_2_pdb $file_chain_2_gro $single_chains_dir

if [ ! -d $double_chain_dir ]; then
    mkdir $double_chain_dir
else 
    rm -r ${double_chain_dir}/*
fi

mv $file_double_chain_gro $file_double_chain_pdb topol.top posre.itp $double_chain_dir
mv $single_chains_dir $double_chain_dir
