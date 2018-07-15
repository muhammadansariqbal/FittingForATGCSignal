#!/bin/bash
# This is a script to create workspaces and cards with different binnings.

for binWid in 100 95 90 85 80 75 70 65 60 55 50 45 40 35 30 25 20 15 10 5
do
	echo "====================================================================="
	echo "BIN WIDTH: "$binWid
	echo "====================================================================="
#	touch workspace_simfit.root
#	touch aC_WWWZ_simfit.txt
	python make_PDF_input_oneCat.py -n -c elmu --binWidth $binWid
	text2workspace.py aC_WWWZ_simfit.txt -o workspace_simfit.root -P CombinedEWKAnalysis.CommonTools.ACModel:par1par2par3_TF3_shape_Model --PO channels=WWWZ_sig_el,WWWZ_sig_mu,WWWZ_sb_lo_el,WWWZ_sb_lo_mu,WWWZ_sb_hi_el,WWWZ_sb_hi_mu --PO poi=cwww,ccw,cb --PO range_cwww=-20,20 --PO range_ccw=-30,30 --PO range_cb=-75,75
	mv workspace_simfit.root "workspace_simfit_binWidth""$binWid"".root"
#	mv aC_WWWZ_simfit.txt "aC_WWWZ_simfit_binWidth""$binWid"".txt"
done
