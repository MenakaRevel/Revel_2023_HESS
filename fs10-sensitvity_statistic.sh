#! /bin/bash

### SET "mool PBS" @ IIS U-Tokyo
#PBS -q E10
#PBS -l select=1:ncpus=10:mem=10gb
#PBS -l place=scatter
#PBS -j oe
#PBS -m ea
#PBS -M menaka@rainbow.iis.u-tokyo.ac.jp
#PBS -V
#PBS -N Figs10

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
eyear=2014
emon=12
eday=31

# figname="fig07-NSE_boxplot"
# figname="fig-KGE_boxplot"
# figname="fig-KGEAI_boxplot"
# figname="fig07-rISS_boxplot"
# figname="fig-DCORR_boxplot"
figname="fig10s-NSEAI_boxplot_dist"

#*** 0. experiment list
EXLIST="./Fig10s-experiment_list.nam"
rm -r $EXLIST
cat >> ${EXLIST} << EOF
Exp A1: ANO_WSE_E2O_HWEB_002
Exp A2: ANO_WSE_E2O_HWEB_001
Exp N1: NOM_WSE_E2O_HWEB_002
Exp N2: NOM_WSE_E2O_HWEB_005
Exp N3: NOM_WSE_E2O_HWEB_001
Exp N4: NOM_WSE_E2O_HWEB_004
EOF

# python src/NSE_boxplot.py $syear $eyear $CaMa_dir $mapname $figname $NCPUS &
# python src/rISS_boxplot.py $syear $eyear $CaMa_dir $mapname $figname $NCPUS &
python src/NSEAI_boxplot_distribution.py $syear $eyear $CaMa_dir $mapname $EXLIST $figname $NCPUS &
# python src/KGE_boxplot.py $syear $eyear $CaMa_dir $mapname $figname $NCPUS &
# python src/KGEAI_boxplot.py $syear $eyear $CaMa_dir $mapname $figname $NCPUS &
# python src/deltacorr_boxplot.py $syear $eyear $CaMa_dir $mapname $figname $NCPUS &

wait

rm -r $EXLIST

# conda deactivate