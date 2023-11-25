#! /bin/bash

### SET "mool PBS" @ IIS U-Tokyo
#PBS -q E10
#PBS -l select=1:ncpus=10:mem=10gb
#PBS -l place=scatter
#PBS -j oe
#PBS -m ea
#PBS -M menaka@rainbow.iis.u-tokyo.ac.jp
#PBS -V
#PBS -N FigureS6

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

# syear, smon, sday
syear=2009
smon=1
sday=1

# eyear, emon, eday
eyear=2014
emon=12
eday=31

figname="figs6-NRMSE_ann_max_map"

#*** 0. experiment list
EXLIST="./Figs6-experiment_list.nam"
rm -r $EXLIST
cat >> ${EXLIST} << EOF
Exp 3a: NOM_WSE_E2O_HWEB_004
Exp 3b: NOM_WSE_E2O_HWEB_006
EOF


python src/discharge_ann_max_map.py $syear $eyear $CaMa_dir $mapname $figname $NCPUS &

wait

rm -r $EXLIST

rm -r $GRDCLIST
# conda deactivate