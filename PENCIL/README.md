# NAME OF THE PROGRAM

## Table of contents

- [Introduciton](#introduction)
- [Installation and Requirements](#installation-and-requirements)
- [Usage](#usage)
- [main_numba.py](#main_numbapy)
- [exe](#exe)
- [post_processing.py](#post_processingpy)
- [Example](#example)
- [Bugs (unfixed)](#bugs-unfixed)
- [Contacts](#contacts)

## Introduction

Write why this program exists and what the various programs do.

## Installation and Requirements

The programs needed to correctly install and run NAME_OF_THE_PROGRAM are the following:

- `python3` : python3 is a very common and powerful high-level object-oriented programming language. On linux platforms it should be already installed, whereas on Windows it must be installed manually. If you do not know how to do it, plenty of tutorials are available on the web.

- `vmd` : vmd is a widespread visualization tool used to visualize protein structures. This tool is optional, and it should be used to visualize the protein structure associating the probability of retaining a residue to each residue (see section [main_numba](#main_numbapy) for more info).

- `gcc` : gcc is the most common and used (family of) compiler(s) for the C programming language. On Linux platforms, it is already installed. If you need it on Windows, it can be installed via mingw.

- `make` : make is a build tool used to build executables form source code in an automatized way. It should be already installed on Linux machines (if not, it should be easily obtainable using the distro's package manager). If you have Windows, you can find tutorials of how to install make online.

> **NOTE**: This program can also be run with just `python3` installed, as it is explained in section [usage](#usage). `gcc` and `make` are required only if you plan to use the C_routine, which is recommended for big (i.e. > 2000 residues) structures.

To run the python scripts, you need the following libraries:

- mdtraj
- numpy
- pandas
- numba
- matplotlib

All these libraries can be easily obtained using tools like Anaconda (or miniconda). You can setup your environment to run NAME_OF_THE_PROGRAM launching the following commands on your Linux terminal (with Anaconda/miniconda installed) :

```bash
conda create --name NAME_OF_THE_PROGRAM
conda activate NAME_OF_THE_PROGRAM
conda install -c conda-forge numba
conda install -c conda-forge matplotlib
conda install -c conda-forge mdtraj
```

As fas as the C routine is concerned, the only required library to run in parallel mode is `openmp`, which is (as far as I know, but do not take it for granted) already installed on Linux machines.  
Compiling the code is very simple, you just have to type these few commands in your terminal (assuming that your initial directory is the one containing this README) :

```bash
cd C_routine/
make
```

After this, you should find ax executable called `exe` in the *C_routine/build* directory.

## Usage

The program is composed of 3 parts:

- The first python script `main_numba.py`
- The C executable `exe`
- The second python script `post_processing.py`

There are 2 possible pipelines, one involving only python and another one which makes use of the C executable.

***

### Just-python pipeline

This pipeline just involves two scripts and they must be run in the following order :

1. `main_numba.py`
2. `post_processing.py`

To use this pipeline, you **must supply the -n flag** to the [main_numba.py](#main_numbapy) script.

To understand how to use the single scripts, please see the sections below.

***

### Mixed pipeline

This pipeline involves all the 3 scripts/executable and they must be run in the following order :

1. `main_numba.py`
2. `exe`
3. `post_processing.py`

To use this pipeline, you **should NOT supply the -n flag** to the [main_numba.py](#main_numbapy) script.

To understand how to use the single scripts, please see the sections below.

## main_numba.py

### Command

The command to launch the main_numba.py script is:

```bash
python3 main_numba.py -c <prot_code> -d <excogito_file_directory> -g <.gro file> [-i <T_0>] [-n <num_propre>] [-p <parallel_bool>] [-r <radius_canvas>] [-s num_steps] [-t decay_time]
```

The arguments can also be passed with long double-dashed labels, which can be found in the Arguments section hereafter.

### Arguments

Here you can find a detailed list of the input arguments :

- `prot_code` : Mandatory input argument. It must be a string and the output of the program will be stored in the directory *output_files/prot_code*. Long flag: `--code`

- `excogito_file_directory` : Mandatory input argument. It must be a string and it must contain the path to the directory where the excogito output files excogito are located. Long flag: `--drectory`

- `.gro file` : Mandatory input argument. It must be a string and it must contain the path to the .gro file of the protein which you analyzed with excogito. Long flag: `--gro`

- `T_0` : Optional input argument. It must be a positive float and represents the starting temperature for the Simulated Annealing procedure. The **default value** is calculated by the program such that the initial acceptance probability is 75%. Long flag: `--initial-T`

- `num_propre` : Optional input argument. This argument must be a positive integer. If this flag is not supplied, the program stops before the simulated annealing optimization and just analyzes the excogito output. If you supply this flag, then you perform the SA optimization and you can directly skip to the execution of post_processing.py, otherwise you have to use the [exe](#exe) program to perform the SA optimization.  
If you performed the PROPRE analysis on your trajectory (and probably used the number of atoms given in output also for the EXCOGITO analysis), then you are strongly encouraged to use that result as value for this flag. If you did not perform that analysis, simply supply the number of atoms which were retained in the excogito mappings. There is **no default value**. Long flag: `--num-propre`

- `parallel_bool` : Optional input argument. This argument must be 0 or 1. The value 0 corresponds to serial execution, hence only one optimization is performed. The option 1 corresponds to parallel execution, and the program performs a number of optimizations equal to the available number of cores. The **default value** is 1. Long flag: `--parallel`

- `radius_canvas` : Optional input argument. This argument must be a positive float and represents the value (in nm) of the radius of the MG region desired in the CANVAS model (choice 2) around each all-atom residue. The **default value** is 0.7. Long flag: `--radius`

- `num_steps` : Optional input argument. It must be a positive integer and represents the number of Monte Carlo steps performed during the simulated annealing procedure. The **default value** is 40'000. Long flag: `--steps`

- `decay_time` : Optional input argument. It must be a positive float and represents the decay time for the exponential decay of the temperature parameter during the simulated annealing optimization. The **default value** is calculated such that the final temperature is 0.1% of the starting temperature. Long flag: `--decay-time`

> **SUGGESTIONS**:  
>
>1. Please, **pay attention to the pipeline** you want to follow! Remember to put the -n flag if you plan to use the just-python pipeline.
>
>2. For the input files, create a directory with the name of your protein (e.g. "4AKE") in the directory *input_files* and store inside *input_files/4AKE* all the excogito output files and the .gro file of the structure. In this way you can just use `-c 4AKE -d ?/input_files/4AKE -g ?/input_files/4AKE/my_file.gro` as mandatory flags (of course the path ? depends on the directory from which you execute the .py script) and the output will be stored in *output_files/4AKE*.  
Nice and clean!
>
>3. If you want to perform a number of optimizations higher than the number of your cores, see the [example](#example) section.

### Output

The output of this script depends on whether the -n flag is used. Some files are produced in any case, whereas other files are ONLY produced with the -n flag active.

#### Always-present output files

The output files that do not depend on the -n flag are the following:

#### -n flag exclusive output files

## exe

### Command

### Arguments

### Output

## post_processing.py

### Command

### Arguments

### Output

## Example

...

## Bugs (unfixed)

No known bugs are present at the moment.

## Contacts

<fabio.passi24@gmail.com> or <fabio.passi@studenti.unitn.it>
