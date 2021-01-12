#!/bin/bash
#$ -cwd
#$ -j y
#$ -q hep.q
#$ -l h_rt=12:0:0 
#$ -o $JOB_ID_$TASK_ID.dat
#$ -t 1-3
 

source /home/hep/vc1117/.cms_setup.sh
export X509_USER_PROXY=/home/hep/vc1117/proxy

cd /home/hep/vc1117/LLP/CMSSW_10_2_18/src
eval `scramv1 runtime -sh`

cd /home/hep/vc1117/LLP/CMSSW_10_2_18/src/PhysicsTools/NanoAODTools/scripts/getEventYields

date
jobs=(2016 2017 2018)
year=${jobs[$SGE_TASK_ID-1]}
echo $year

./getEventYields /vols/cms/LLP/files_201117/analysis/$year $year

date
