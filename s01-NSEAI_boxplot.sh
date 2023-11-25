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
cd "/cluster/data7/menaka/Reveletal2021b"

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

# figname="fig9-NSEAI_boxplot"

# python src/NSEAI_boxplot.py $syear $eyear $CaMa_dir $mapname $figname $NCPUS &

# figname="fig10-deltacorr_boxplot"

# python src/deltacorr_boxplot.py $syear $eyear $CaMa_dir $mapname $figname $NCPUS &

figname="fig11-NSE_boxplot"

python src/NSE_boxplot.py $syear $eyear $CaMa_dir $mapname $figname $NCPUS &

# figname="fig12-corr_boxplot"

# python src/corr_boxplot.py $syear $eyear $CaMa_dir $mapname $figname $NCPUS &

# figname="fig13-KGE_boxplot"

# python src/KGE_boxplot.py $syear $eyear $CaMa_dir $mapname $figname $NCPUS &

# figname="fig14-ISS_boxplot"

# python src/ISS_boxplot.py $syear $eyear $CaMa_dir $mapname $figname $NCPUS &

wait

# conda deactivate