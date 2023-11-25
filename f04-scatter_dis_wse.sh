#! /bin/bash

### SET "mool PBS" @ IIS U-Tokyo
#PBS -q E10
#PBS -l select=1:ncpus=10:mem=10gb
#PBS -l place=scatter
#PBS -j oe
#PBS -m ea
#PBS -M menaka@rainbow.iis.u-tokyo.ac.jp
#PBS -V
#PBS -N Figure4

#source ~/.bashrc
# import virtual environment
# source ~/.bsash_conda
# source activate pydef

which python

# OMP Settings
NCPUS=10
export OMP_NUM_THREADS=$NCPUS

# cd $PBS_O_WORKDIR
cd "/cluster/data7/menaka/Reveletal2022"

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

mkdir -p figures

figname="fig04-EXP1_dis_wse_hydrograph"

#*** 0. experiment list
EXLIST="./Fig04-experiment_list.nam"
rm -r $EXLIST
cat >> ${EXLIST} << EOF
Exp 1: DIR_WSE_E2O_HWEB_001
EOF

EXPNAME="DIR_WSE_E2O_HWEB_001"

#*** 1. GRDC list
GRDCLIST="./Fig04-GRDClist.txt"
rm -r $GRDCLIST
cat >> ${GRDCLIST} << EOF
3625340
3624121
3620000
EOF
# 3625340 - PURUS, RIO; LABREA
# 3624121 - JURUA, RIO; SANTOS DUMONT
# 3620000 - AMAZON; SANTO ANTONIO DO ICA



# 3625350
# 3629001
# 3627010
# 3618000
#======
# python src/scatter_dis_wse.py $syear $eyear 2 "DIR_WSE_E2O_HWEB_001" "DIR_WSE_E2O_HWEB_002" $CaMa_dir $mapname $figname &
# python src/scatter_dis_wse.py $syear $eyear $CaMa_dir $mapname $figname $EXLIST $GRDCLIST $NCPUS &
python src/map_dis_wse_hydrograph.py $syear $eyear $CaMa_dir $mapname $figname $EXPNAME $GRDCLIST $NCPUS &

wait

rm -r $EXLIST

rm -r $GRDCLIST

# conda deactivate