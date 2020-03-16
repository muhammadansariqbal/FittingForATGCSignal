#!/bin/bash
# This is a script to run limit plots

for year in 16 17 18 Run2
do
	python plot1D_limit.py --POI cwww --pval 3.6 --year "$year"
	python plot1D_limit.py --POI ccw --pval 4.5 --year "$year"
	python plot1D_limit.py --POI cb --pval 20 --year "$year"

	python plot1D_limit.py --POI lZ --pval 0.014 --year "$year"
	python plot1D_limit.py --POI dg1z --pval 0.018 --year "$year"
	python plot1D_limit.py --POI dkz --pval 0.02 --year "$year"
done

for year in 16 17 18 Run2
do
	python plot2D_limit.py --POI cwww,ccw --pval 3.6,4.5 --year "$year"
	python plot2D_limit.py --POI cwww,cb --pval 3.6,20 --year "$year"
	python plot2D_limit.py --POI ccw,cb --pval 4.5,20 --year "$year"

	python plot2D_limit.py --POI lZ,dg1z --pval 0.014,0.018 --year "$year"
	python plot2D_limit.py --POI lZ,dkz --pval 0.014,0.02 --year "$year"
	python plot2D_limit.py --POI dg1z,dkz --pval 0.018,0.02 --year "$year"
done

