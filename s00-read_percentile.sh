#! /bin/bash

### SET "mool PBS" @ IIS U-Tokyo
#PBS -q F20
#PBS -l select=1:ncpus=20:mem=100gb
#PBS -l place=scatter
#PBS -j oe
#PBS -m ea
#PBS -M menaka@rainbow.iis.u-tokyo.ac.jp
#PBS -V
#PBS -N read_percentile

#source ~/.bashrc
# import virtual environment
source ~/.bashrc
source ~/.bash_conda

source activate pydef

# OMP Settings
NCPUS=20
export OMP_NUM_THREADS=$NCPUS

# cd $PBS_O_WORKDIR
cd "/cluster/data7/menaka/Reveletal2022"

# Input Dir
indir="/cluster/data7/menaka/HydroDA/out"

# out dir
outdir="/cluster/data7/menaka/Reveletal2022/txt"


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

# N=`python src/calc_days.py $syear $smon $sday $eyear $emon $eday`

# expname="NOM_WSE_E2O_HWEB_001"
# expname="NOM_WSE_E2O_HWEB_002"
# expname="NOM_WSE_E2O_HWEB_003"
# expname="NOM_WSE_E2O_HWEB_004"
# expname="NOM_WSE_E2O_HWEB_005"
# expname="NOM_WSE_E2O_HWEB_006"
expname="NOM_WSE_E2O_HWEB_007"
# expname="ANO_WSE_E2O_HWEB_001"
# expname="ANO_WSE_E2O_HWEB_002"
# expname="ANO_WSE_E2O_HWEB_003"
# expname="ANO_WSE_E2O_HWEB_004"
# expname="ANO_WSE_E2O_HWEB_005"
# expname="ANO_WSE_E2O_HWEB_006"
# expname="DIR_WSE_E2O_HWEB_001"
# expname="DIR_WSE_E2O_HWEB_002"
# expname="DIR_WSE_E2O_HWEB_003"

ens_mem=49

rm -r $outdir/$expname/Q_percentile
mkdir -p $outdir/$expname/Q_percentile

python ./src/read_percentile.py $syear $eyear $CaMa_dir $mapname $expname $NCPUS &

wait

conda deactivate