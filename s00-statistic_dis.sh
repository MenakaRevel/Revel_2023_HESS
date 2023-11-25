#! /bin/bash

### SET "mool PBS" @ IIS U-Tokyo
#PBS -q F10
#PBS -l select=1:ncpus=10:mem=100gb
#PBS -l place=scatter
#PBS -j oe
#PBS -m ea
#PBS -M menaka@rainbow.iis.u-tokyo.ac.jp
#PBS -V
#PBS -N statistic

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

printf '%s%10s%10s%10s%10s' id r NSE KGE sharpness > tmp.txt
python src/statistic_dis.py $syear $eyear $expname $CaMa_dir $mapname >> tmp.txt

mkdir -p statistic/discharge

mv tmp.txt statistic/discharge/$expname.txt
echo "Saving ..."
echo statistic/discharge/$expname.txt