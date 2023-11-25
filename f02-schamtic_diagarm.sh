#! /bin/bash

### SET "mool PBS" @ IIS U-Tokyo
#PBS -q E10
#PBS -l select=1:ncpus=10:mem=100gb
#PBS -l place=scatter
#PBS -j oe
#PBS -m ea
#PBS -M menaka@rainbow.iis.u-tokyo.ac.jp
#PBS -V
#PBS -N Figure2

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
# CaMa_dir="/cluster/data6/menaka/CaMa-Flood_v396a_20200514"
CaMa_dir="/cluster/data6/menaka/CaMa-Flood_v4"

# Map
mapname="amz_06min"
# mapname="glb_06min"

# syear, smon, sday
syear=2009
smon=1
sday=1

# eyear, emon, eday
eyear=2014
emon=12
eday=31

figname="fig02-schemetic_diagram"

# bias
# a R_AMAZONAS_AMAZONAS_KM0842 
namea="R_AMAZONAS_JARI_KM0529"

# amplitude differnce
# b R_AMAZONAS_GRANDE_KM3669 
# R_AMAZONAS_JAPURA_KM2981
# R_AMAZONAS_JARI_KM0241
# R_AMAZONAS_JARI_KM0336
# R_AMAZONAS_JAVARI_KM3655
# R_AMAZONAS_JAVARI_KM3717
# nameb="R_AMAZONAS_JUTAI_KM3325"
nameb="R_AMAZONAS_JUTAI_KM3182"

indir="/cluster/data6/menaka/ensemble_org/CaMa_out/AMZE2O003"
# indir="/home/yamadai/data/CaMa_v400_simulations/VIC_BC_3h_06min"

python src/bias_amplitude_wse.py $syear $eyear $namea $nameb $CaMa_dir $mapname $figname $indir $NCPUS &

wait

# conda deactivate