#!/bin/bash
# This is a script to run combine for different mWV cutsoffs

cd ..

for cutoff in 2400 2700 3000 3300 3600 3900 4200
do
	echo "====================================================================="
	echo "CUTOFF mWV: "$cutoff GeV
	echo "====================================================================="

	combine "cutoffLimits/workspace_simfit_cutoff""$cutoff"".root" -M MultiDimFit --floatOtherPOIs=0 --algo=grid --expectSignal=1 --points=100 --redefineSignalPOIs cwww -P cwww --freezeNuisances ccw,cb --setPhysicsModelParameters cwww=0,ccw=0,cb=0 --setPhysicsModelParameterRange cwww=-3.6,3.6 --minimizerStrategy=2 --cminPreScan -n "_cwww_3.6_cutoff""$cutoff"
	combine "cutoffLimits/workspace_simfit_cutoff""$cutoff"".root" -M MultiDimFit --floatOtherPOIs=0 --algo=grid --expectSignal=1 --points=100 --redefineSignalPOIs ccw -P ccw --freezeNuisances cwww,cb --setPhysicsModelParameters cwww=0,ccw=0,cb=0 --setPhysicsModelParameterRange ccw=-4.5,4.5 --minimizerStrategy=2 --cminPreScan -n "_ccw_4.5_cutoff""$cutoff"
	combine "cutoffLimits/workspace_simfit_cutoff""$cutoff"".root" -M MultiDimFit --floatOtherPOIs=0 --algo=grid --expectSignal=1 --points=100 --redefineSignalPOIs cb -P cb --freezeNuisances cwww,ccw --setPhysicsModelParameters cwww=0,ccw=0,cb=0 --setPhysicsModelParameterRange cb=-20,20 --minimizerStrategy=2 --cminPreScan -n "_cb_20_cutoff""$cutoff"

	mv "higgsCombine_cwww_3.6_cutoff""$cutoff"".MultiDimFit.mH120.root" cutoffLimits
	mv "higgsCombine_ccw_4.5_cutoff""$cutoff"".MultiDimFit.mH120.root" cutoffLimits
	mv "higgsCombine_cb_20_cutoff""$cutoff"".MultiDimFit.mH120.root" cutoffLimits

	combine "cutoffLimits/workspace_simfit_LEP_cutoff""$cutoff"".root" -M MultiDimFit --floatOtherPOIs=0 --algo=grid --expectSignal=1 --points=100 --redefineSignalPOIs lZ -P lZ --freezeNuisances dg1z,dkz --setPhysicsModelParameters lZ=0,dg1z=0,dkz=0 --setPhysicsModelParameterRange lZ=-0.014,0.014 --minimizerStrategy=2 --cminPreScan -n "_lZ_0.014_cutoff""$cutoff"
	combine "cutoffLimits/workspace_simfit_LEP_cutoff""$cutoff"".root" -M MultiDimFit --floatOtherPOIs=0 --algo=grid --expectSignal=1 --points=100 --redefineSignalPOIs dg1z -P dg1z --freezeNuisances lZ,dkz --setPhysicsModelParameters lZ=0,dg1z=0,dkz=0 --setPhysicsModelParameterRange dg1z=-0.018,0.018 --minimizerStrategy=2 --cminPreScan -n "_dg1z_0.018_cutoff""$cutoff"
	combine "cutoffLimits/workspace_simfit_LEP_cutoff""$cutoff"".root" -M MultiDimFit --floatOtherPOIs=0 --algo=grid --expectSignal=1 --points=100 --redefineSignalPOIs dkz -P dkz --freezeNuisances lZ,dg1z --setPhysicsModelParameters lZ=0,dg1z=0,dkz=0 --setPhysicsModelParameterRange dkz=-0.02,0.02 --minimizerStrategy=2 --cminPreScan -n "_dkz_0.02_cutoff""$cutoff"
	
	mv "higgsCombine_lZ_0.014_cutoff""$cutoff"".MultiDimFit.mH120.root" cutoffLimits
	mv "higgsCombine_dg1z_0.018_cutoff""$cutoff"".MultiDimFit.mH120.root" cutoffLimits
	mv "higgsCombine_dkz_0.02_cutoff""$cutoff"".MultiDimFit.mH120.root" cutoffLimits

done

