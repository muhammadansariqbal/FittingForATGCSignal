#!/bin/bash
# This is a script to run combine for different workspaces (Useful for bias testing)

for binWid in 95 85 75 65 55 45 35 25 15 5 1
do
	echo "====================================================================="
	echo "BIN WIDTH: "$binWid
	echo "====================================================================="
	combine "workspace_simfit_binWidth""$binWid"".root" -M MultiDimFit --floatOtherPOIs=0 --algo=grid --expectSignal=1 --points=4 --redefineSignalPOIs cwww,ccw -P cwww -P ccw --freezeNuisances cb --setPhysicsModelParameters cwww=0,ccw=0,cb=0 --setPhysicsModelParameterRange cwww=-3.6,3.6:ccw=-4.5,4.5 --minimizerStrategy=2 --cminPreScan -t -1 -n "_cwww_3.6_ccw_4.5_binWidth""$binWid"
	combine "workspace_simfit_binWidth""$binWid"".root" -M MultiDimFit --floatOtherPOIs=0 --algo=grid --expectSignal=1 --points=4 --redefineSignalPOIs cwww,cb -P cwww -P cb --freezeNuisances ccw --setPhysicsModelParameters cwww=0,ccw=0,cb=0 --setPhysicsModelParameterRange cwww=-3.6,3.6:cb=-20,20 --minimizerStrategy=2 --cminPreScan -t -1 -n "_cwww_3.6_cb_20_binWidth""$binWid"
	combine "workspace_simfit_binWidth""$binWid"".root" -M MultiDimFit --floatOtherPOIs=0 --algo=grid --expectSignal=1 --points=4 --redefineSignalPOIs ccw,cb -P ccw -P cb --freezeNuisances cwww --setPhysicsModelParameters cwww=0,ccw=0,cb=0 --setPhysicsModelParameterRange ccw=-4.5,4.5:cb=-20,20 --minimizerStrategy=2 --cminPreScan -t -1 -n "_ccw_4.5_cb_20_binWidth""$binWid"

done
