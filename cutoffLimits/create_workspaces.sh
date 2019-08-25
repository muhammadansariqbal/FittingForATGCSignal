#!/bin/bash
# This is a script to create workspaces and cards with different mWV cutoffs.

cd ..
for cutoff in 2100 2400 2700 3000 3300 3600 3900 4200 4500
do
	echo "====================================================================="
	echo "CUTOFF mWV: "$cutoff GeV
	echo "====================================================================="

	cp "cutoffLimits/wwlvj_el_HPV_900_""$cutoff""_workspace.root" "Input/wwlvj_el_HPV_workspace.root"
	cp "cutoffLimits/wwlvj_mu_HPV_900_""$cutoff""_workspace.root" "Input/wwlvj_mu_HPV_workspace.root"

	python make_PDF_input_oneCat.py -n -c elmu --cutoff $cutoff
	text2workspace.py aC_WWWZ_simfit.txt -o workspace_simfit.root -P CombinedEWKAnalysis.CommonTools.ACModel:par1par2par3_TF3_shape_Model --PO channels=WWWZ_sig_el,WWWZ_sig_mu,WWWZ_sb_lo_el,WWWZ_sb_lo_mu,WWWZ_sb_hi_el,WWWZ_sb_hi_mu --PO poi=cwww,ccw,cb --PO range_cwww=-20,20 --PO range_ccw=-30,30 --PO range_cb=-75,75

	mv "workspace_simfit.root" "cutoffLimits/workspace_simfit_cutoff""$cutoff"".root"

done

