#!/bin/bash

##runs a number of SCORE jobs one after another
##Markus Kuehbach, 2019/01/20, m.kuehbach at mpie.de


exenm="scorehybrid_rhcs_mmm2018_intel_febr2019_maws15"
ebsdnm="KAM_5deg_1NN.txt_MTEX_X2EastZInPlane.SEG2Deg.ang"
zero="0"
multiplier=111111


simid=64031
seed=1
targetseed=$(($multiplier * $seed))
uds=SCORE.Debug.MartinDiehl.64.SimID.6403.uds
mpirun -np 1 ./${exenm} ${simid} ${uds} ${ebsdnm} ${targetseed} 1>SCORE.SimID.${simid}.STDOUT.txt 2>SCORE.SimID.${simid}.STDERR.txt
