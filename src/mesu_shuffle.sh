#!/bin/bash
#PBS -q alpha
#PBS -l select=1:ncpus=1
#PBS -N ams_shuffle

#load appropriate modules
#module load python-3.6

# Define scratch space scratch alpha for UV
SCRATCH=/scratchalpha/$USER/myprogram_scratch_space
PROJECT=ALLOGENOMICS
mkdir $SCRATCH
mkdir $SCRATCH/$PROJECT

# activate conda environnment
source activate conda_env

# copy ams directory to working directory
#cp -r AMS_workflow $PBS_O_WORKDIR/AMS_workflow