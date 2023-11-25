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
# source ~/.bashrc
# source ~/.bash_conda

# source activate pydef

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
eyear=2009
emon=12
eday=31

# figname="fig07-NSE_boxplot"
figname="figs9-dis_max_min_scatter"

#*** 0. experiment list
EXLIST="./Figs9-experiment_list.nam"
rm -r $EXLIST
cat >> ${EXLIST} << EOF
NOM: NOM_WSE_E2O_HWEB_004
NOM_2: NOM_WSE_E2O_HWEB_006
EOF

# Exp 3b: NOM_WSE_E2O_HWEB_007
# Exp 1a: DIR_WSE_E2O_HWEB_001
# Exp 1b: DIR_WSE_E2O_HWEB_002
# Exp 2a: ANO_WSE_E2O_HWEB_004
# Exp 2b: ANO_WSE_E2O_HWEB_005

python src/discharge_ann_max_min_log.py $syear $eyear $CaMa_dir $mapname $EXLIST $figname $NCPUS &

wait

rm -r $EXLIST

# conda deactivate