#! /bin/bash

### SET "mool PBS" @ IIS U-Tokyo
#PBS -q F10
#PBS -l select=1:ncpus=10:mem=100gb
#PBS -l place=scatter
#PBS -j oe
#PBS -m ea
#PBS -M menaka@rainbow.iis.u-tokyo.ac.jp
#PBS -V
#PBS -N multi_read_dis

# OMP Settings
NCPUS=10
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

N=`python src/calc_days.py $syear $smon $sday $eyear $emon $eday`

ens_mem=49

# for expname in "DIR_WSE_ECMWF_HWEB_011" "DIR_WSE_ECMWF_HWEB_012" "DIR_WSE_ECMWF_HWEB_013" "DIR_WSE_ECMWF_HWEB_014";
# for expname in "ANO_WSE_ECMWF_HWEB_012" "ANO_WSE_ECMWF_HWEB_013" "NOM_WSE_ECMWF_HWEB_012" "NOM_WSE_ECMWF_HWEB_013";
# for expname in "ANO_WSE_ECMWF_HWEB_011" "NOM_WSE_ECMWF_HWEB_011";
for expname in "DIR_WSE_ECMWF_HWEB_011" "ANO_WSE_ECMWF_HWEB_012"
do
    mkdir -p $outdir/$expname/outflow

    echo ./src/read_discharge $expname $mapname $syear $smon $sday $eyear $emon $eday $ens_mem $N $CaMa_dir $indir $outdir
    time ./src/read_discharge $expname $mapname $syear $smon $sday $eyear $emon $eday $ens_mem $N $CaMa_dir $indir $outdir &
        ## for parallel computation using multiple CPUs 
        NUM=`ps -U $USER | grep ./src/read_discharge | wc -l | awk '{print $1}'`
        while [ $NUM -gt $NCPUS ];
        do
            sleep 1
            NUM=`ps -U $USER | grep ./src/read_discharge | wc -l | awk '{print $1}'`
        done
        # enum=$(( $enum + 1 ))
done

wait