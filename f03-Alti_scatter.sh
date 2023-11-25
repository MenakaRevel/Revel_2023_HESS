#! /bin/bash

### SET "mool PBS" @ IIS U-Tokyo
#PBS -q E10
#PBS -l select=1:ncpus=10:mem=100gb
#PBS -l place=scatter
#PBS -j oe
#PBS -m ea
#PBS -M menaka@rainbow.iis.u-tokyo.ac.jp
#PBS -V
#PBS -N Figure3

#source ~/.bashrc
# import virtual environment
# source ~/.bashrc
# source ~/.bash_conda

# source activate pydef

which python

# OMP Settings
NCPUS=10
export OMP_NUM_THREADS=$NCPUS

# got to working dirctory
# cd $PBS_O_WORKDIR
cd "/cluster/data7/menaka/Reveletal2022"

mkdir -p figures
mkdir -p data

#CaMA-Flood directory
CaMa_dir="/cluster/data6/menaka/CaMa-Flood_v4"

# Map
mapname="amz_06min"

# syear, smon, sday
syear=2009
smon=1
sday=1

# eyear, emon, eday
eyear=2014
emon=12
eday=31

obslist="/cluster/data6/menaka/HydroDA/dat/HydroWeb_alloc_amz_06min_QC0.txt"

obslist_sim="/cluster/data6/menaka/HydroDA/dat/HydroWeb_alloc_amz_06min_QC0_simulation.txt"

obslist_val="/cluster/data6/menaka/HydroDA/dat/HydroWeb_alloc_amz_06min_QC0_validation.txt"

figname="fig03-altimetry_map"

python src/Alti_scatter.py  $syear $eyear $CaMa_dir $mapname $obslist_sim $obslist_val $figname $NCPUS &

wait

# conda deactivate