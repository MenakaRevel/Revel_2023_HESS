#! /bin/bash

### SET "mool PBS" @ IIS U-Tokyo
#PBS -q E10
#PBS -l select=1:ncpus=10:mem=10gb
#PBS -l place=scatter
#PBS -j oe
#PBS -m ea
#PBS -M menaka@rainbow.iis.u-tokyo.ac.jp
#PBS -V
#PBS -N Figure9

#source ~/.bashrc
# import virtual environment
source ~/.bashrc
source ~/.bash_conda

source activate pydef

which python

# OMP Settings
NCPUS=40
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

# syear, smon, sday
syear=2009
smon=1
sday=1

# eyear, emon, eday
eyear=2014
emon=12
eday=31

figname="fig09-scatter_stat_wse"


#*** 0. experiment list
EXLIST="./Fig09-experiment_list.nam"
rm -r $EXLIST
cat >> ${EXLIST} << EOF
DIR: DIR_WSE_E2O_HWEB_001
ANO: ANO_WSE_E2O_HWEB_004
NOM: NOM_WSE_E2O_HWEB_004
EOF

# Exp 1: DIR_WSE_E2O_HWEB_001
# Exp 2: ANO_WSE_E2O_HWEB_004
# Exp 3: NOM_WSE_E2O_HWEB_004

python src/scatter_stat_wse.py $syear $eyear $CaMa_dir $mapname $figname $NCPUS &

wait

rm -r $EXLIST