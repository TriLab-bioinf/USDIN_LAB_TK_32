#!/usr/bin/bash 

str=$0
MYPATH="${str%/*}"

ln -s ${MYPATH}/adapters.fa
ln -s ${MYPATH}/aux-1_remove_chrM.sh
ln -s ${MYPATH}/aux1-exomedepth.R
ln -s ${MYPATH}/bruce_pipeline.sh
ln -s ${MYPATH}/check_packages.R
cp ${MYPATH}/config.cfg .
ln -s ${MYPATH}/run_exomeDepth.R


