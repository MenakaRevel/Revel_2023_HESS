#! /bin/bash

### SET "mool PBS" @ IIS U-Tokyo
#PBS -q E10
#PBS -l select=1:ncpus=10:mem=10gb
#PBS -l place=scatter
#PBS -j oe
#PBS -m ea
#PBS -M menaka@rainbow.iis.u-tokyo.ac.jp
#PBS -V
#PBS -N Figure11

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
eyear=2009
emon=12
eday=31

figname="fig11-compare_runoff_bathymetry"

#*** 1. experiment list
EXLIST1="./Fig11-experiment_list_1.nam"
rm -r $EXLIST1
cat >> ${EXLIST1} << EOF
Exp D11: DIR_WSE_ECMWF_HWEB_011
Exp A11: ANO_WSE_ECMWF_HWEB_011
Exp N11: NOM_WSE_ECMWF_HWEB_011
EOF

#*** 2. experiment list
EXLIST2="./Fig11-experiment_list_2.nam"
rm -r $EXLIST2
cat >> ${EXLIST2} << EOF
Exp D12: DIR_WSE_ECMWF_HWEB_012
Exp A12: ANO_WSE_ECMWF_HWEB_012
Exp N12: NOM_WSE_ECMWF_HWEB_012
EOF

#*** 3. experiment list
EXLIST3="./Fig11-experiment_list_3.nam"
rm -r $EXLIST3
cat >> ${EXLIST3} << EOF
Exp D13: DIR_WSE_ECMWF_HWEB_013
Exp A13: ANO_WSE_ECMWF_HWEB_013
Exp N13: NOM_WSE_ECMWF_HWEB_013
EOF

#*** 4. experiment list
EXLIST4="./Fig11-experiment_list_4.nam"
rm -r $EXLIST4
cat >> ${EXLIST4} << EOF
Exp D14: DIR_WSE_ECMWF_HWEB_014
Exp A14: ANO_WSE_ECMWF_HWEB_014
Exp N14: NOM_WSE_ECMWF_HWEB_014
EOF

python src/boxplot_NSEAI_runoff_bathymetry.py $syear $eyear $CaMa_dir $mapname $EXLIST1 $EXLIST2 $EXLIST3 $EXLIST4 $figname $NCPUS &

wait

rm -r $EXLIST1
rm -r $EXLIST2
rm -r $EXLIST3
rm -r $EXLIST4

# conda deactivate