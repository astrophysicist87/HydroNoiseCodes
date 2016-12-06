#! /usr/bin/env bash

#do1sp=1	# 1 - true
#doHBT=0	# 0 - false
#
#for chosen_rc in 0.001 0.0025 0.005 0.01 0.025 0.05 0.1 0.25 0.5 1.0
#do
#	for chosen_trajectory in 1 2 3
#	do
#		for chosen_tau_D in 0.01 0.02 0.03 0.04 0.05
#		do
#			echo $chosen_rc $chosen_trajectory $chosen_tau_D
#			./HNC $chosen_trajectory $do1sp $doHBT $chosen_rc $chosen_tau_D >> all_results_1sp.out
#		done
#	done
#done

for chosen_trajectory in 1 2 3
do
	for chosen_tau_D in 0.01 0.025 0.05 0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.5 0.55 0.6 0.65 0.7 0.75 0.8 0.85 0.9 0.95 1.0
	do
		echo $chosen_trajectory $chosen_tau_D
		./HNC $chosen_trajectory 1 0 $chosen_tau_D >> all_results_1sp_pointB.out
		./HNC $chosen_trajectory 0 1 $chosen_tau_D >> all_results_HBT_pointB.out
	done
done
