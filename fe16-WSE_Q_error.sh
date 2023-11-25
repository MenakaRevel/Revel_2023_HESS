#! /bin/bash

### SET "mool PBS" @ IIS U-Tokyo
#PBS -q G12
#PBS -l select=1:ncpus=12:mem=10gb
#PBS -l place=scatter
#PBS -j oe
#PBS -m ea
#PBS -M menaka@rainbow.iis.u-tokyo.ac.jp
#PBS -V
#PBS -N Figure

#source ~/.bashrc
# import virtual environment
# source ~/.bashrc
# source ~/.bash_conda

# source activate pydef

which python

# OMP Settings
NCPUS=12
export OMP_NUM_THREADS=$NCPUS

# got to working dirctory
# cd $PBS_O_WORKDIR
cd "/work/a06/menaka/Revel_2023_HESS/"

mkdir -p pdfdoc

python src/NSEAI_All_WSE_Q_error.py

wait

# conda deactivate