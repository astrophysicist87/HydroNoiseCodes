#! /usr/bin/env bash

rm all_results.out

do1sp=1	# 1 - true
doHBT=0	# 0 - false

for chosen_rc in 0.001 0.0025 0.005 0.01 0.025 0.05 0.1 0.25 0.5 1.0
do
	for chosen_trajectory in 1 2 3
	do
		for chosen_tau_D in 0.01 0.02 0.03 0.04 0.05
		do
			echo $chosen_rc $chosen_trajectory $chosen_tau_D
			./HNC $chosen_trajectory $do1sp $doHBT $chosen_rc $chosen_tau_D >> all_results.out
		done
	done
done
