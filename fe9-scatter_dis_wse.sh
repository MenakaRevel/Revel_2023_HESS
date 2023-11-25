#! /bin/bash

### SET "mool PBS" @ IIS U-Tokyo
#PBS -q F10
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

figname="fig09-EXP3b_dis_wse_hydrograph"

#*** 0. experiment list
EXLIST="./Fig09-experiment_list.nam"
rm -r $EXLIST
cat >> ${EXLIST} << EOF
Exp 3b: NOM_WSE_E2O_HWEB_006
EOF

EXPNAME="NOM_WSE_E2O_HWEB_006"

#*** 1. GRDC list
GRDCLIST="./Fig09-GRDClist.txt"
rm -r $GRDCLIST
cat >> ${GRDCLIST} << EOF
3627035
3625000
3625320
EOF

# 3627035 - MADEIRA, RIO; HUMAITA
# 3625000 - AMAZON; ITAPEUA
# 3625320 - PURUS, RIO; CANUTAMA

#======
# python src/scatter_dis_wse.py $syear $eyear 2 "NOM_WSE_E2O_HWEB_004" "NOM_WSE_E2O_HWEB_006" $CaMa_dir $mapname $figname &
# python src/scatter_dis_wse.py $syear $eyear $CaMa_dir $mapname $figname $EXLIST $GRDCLIST $NCPUS &
python src/map_dis_wse_hydrograph.py $syear $eyear $CaMa_dir $mapname $figname $EXPNAME $GRDCLIST $NCPUS &

wait

rm -r $EXLIST

rm -r $GRDCLIST

conda deactivate