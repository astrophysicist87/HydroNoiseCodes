#! /usr/bin/env bash

outfilename="History/submit_jobs_record_`date +%F`.out"
outfile=`get_filename $outfilename`		# get_filename checks to see if outfilename already exists
										# and attaches an appropriate index to the end to ensure that it is a new file

for trajectory in 1 2 3
do
	for tauColor in 1.00 0.50 0.25 0.10 0.05 0.01
	do
		#pions
		outfile1=results_1sp/pion/w_delta_muB/trajectory_`echo $trajectory`_tauC_`echo $tauColor`.dat
		outfile2=results_HBT/pion/w_delta_muB/trajectory_`echo $trajectory`_tauC_`echo $tauColor`.dat
		qsub -v v0=`pwd`,v1=$trajectory,v2=1,v3=0,v4=$tauColor,v5=1,v6=$outfile1 run.pbs >> $outfile
		qsub -v v0=`pwd`,v1=$trajectory,v2=0,v3=1,v4=$tauColor,v5=1,v6=$outfile2 run.pbs >> $outfile
		#protons
		outfile3=results_1sp/proton/trajectory_`echo $trajectory`_tauC_`echo $tauColor`.dat
		outfile4=results_HBT/proton/trajectory_`echo $trajectory`_tauC_`echo $tauColor`.dat
		qsub -v v0=`pwd`,v1=$trajectory,v2=1,v3=0,v4=$tauColor,v5=2,v6=$outfile3 run.pbs >> $outfile
		qsub -v v0=`pwd`,v1=$trajectory,v2=0,v3=1,v4=$tauColor,v5=2,v6=$outfile4 run.pbs >> $outfile
	done
done
