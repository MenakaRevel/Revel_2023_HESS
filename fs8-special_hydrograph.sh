#! /bin/bash

### SET "mool PBS" @ IIS U-Tokyo
#PBS -q E10
#PBS -l select=1:ncpus=10:mem=10gb
#PBS -l place=scatter
#PBS -j oe
#PBS -m ea
#PBS -M menaka@rainbow.iis.u-tokyo.ac.jp
#PBS -V
#PBS -N FigureS8

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

figname="figs8-hydrograph_compare"

#*** 0. experiment list
EXLIST="./Figs8-experiment_list.nam"
rm -r $EXLIST
cat >> ${EXLIST} << EOF
Open loop: NOM_WSE_E2O_HWEB_004
Assimimated: NOM_WSE_E2O_HWEB_004
EOF

#*** 1. GRDC list
GRDCLIST="./Figs8-GRDClist.txt"
rm -r $GRDCLIST
cat >> ${GRDCLIST} << EOF
3626000
3629001
EOF
# 3620000
# 3627035

#3625000
#3627035

python src/hydrograph_com.py $syear $eyear $CaMa_dir $mapname $figname $EXLIST $GRDCLIST $NCPUS &

wait

rm -r $EXLIST

rm -r $GRDCLIST
# conda deactivate