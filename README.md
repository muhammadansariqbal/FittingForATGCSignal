SIGNAL PARAMETRIZATION
======================

```
This script creates the complete signal model and creates a datacard and a workspace per channel containing everything needed for the limit extraction.

Setup Instructions
------------------

# Setup CMSSW in ~/private/FittingForATGC/Signal/
cmsrel CMSSW_7_1_5
cd CMSSW_7_1_5/src
cmsenv

# Clone the repository
git clone git@github.com:muhammadansariqbal/FittingForATGCSignal.git

# Clone the Higgs combined limit tool repo; Look at the repo for further details
git clone git@github.com:muhammadansariqbal/HiggsAnalysis.git

# Clone the CombinedEWKAnalysis repo; Look at the repo for further details
git clone git@github.com:muhammadansariqbal/CombinedEWKAnalysis.git

# Build
scram b -j 20

# Make a folder Input in and FittingForSignal and copy the required files e.g.
cd FittingForATGCSignal; mkdir Input;
cp /afs/cern.ch/work/m/maiqbal/private/aTGC/Samples_80X_Working/WW-aTGC_mu.root ./Input
cp ../../../../Background/CMSSW_5_3_32/src/FittingForATGCBackground/cards_mu_HPV_900_3500/wwlvj_mu_HPV_workspace.root ./Input
# etc.

# Run the main script
python make_PDF_input_oneCat.py -n -c mu -p --savep
-n: Read the input trees and create RooDataHists(-> faster access); Needed at the first run or when the input trees are changed.
-c {channel}: Only run for {channel} (mu or el)
-p: Make plots
--savep: Save the plots
-b: Run in batch mode
--noatgcint: Set aTGC-interference terms to zero
--printatgc: Print the coefficients of the signal model
--atgc: Using different parametrization (Lagrangian approach instead of EFT)

# The workspaces for the different channels can now be combined with
text2workspace.py aC_WWWZ_simfit.txt -o workspace_simfit.root -P CombinedEWKAnalysis.CommonTools.ACModel:par1par2par3_TF3_shape_Model --PO channels=WWWZ_sig_el,WWWZ_sig_mu,WWWZ_sb_lo_el,WWWZ_sb_lo_mu,WWWZ_sb_hi_el,WWWZ_sb_hi_mu --PO poi=cwww,ccw,cb --PO range_cwww=-20,20 --PO range_ccw=-30,30 --PO range_cb=-75,75
-o: Name of the created workspace
-P: Name of the used model
--PO channels= Names of the channels
--PO poi= Names of the paramters of interest
--PO range_= Set paramter range (does't work atm but has to be added to avoid error message)

# This creates the final workspace called workspace_simfit.root. To inspect this workspace in ROOT you have to load the combined limit .so 
.L ../../lib/slc6_amd64_gcc481/libHiggsAnalysisCombinedLimit.so.

------------------

Limit Calculation
=================

# We can now run combine with (e.g. for 1d limits on cwww)
combine workspace_simfit.root -M MultiDimFit --floatOtherPOIs=0 --algo=grid --expectSignal=1 --points=1000 --redefineSignalPOIs cwww -P cwww --freezeNuisances ccw,cb --setPhysicsModelParameters cwww=0,ccw=0,cb=0 --setPhysicsModelParameterRange cwww=-10,10 --minimizerStrategy=2 -n Example
--points: is the number of scanned parameter values
--redefineSignalPOIs cwww -P cwww: defines the paramter(s) of interest
--freezeNuisances ccw,cb: fixes the other paramters
--setPhysicsModelParameters cwww=0,ccw=0,cb=0: sets the initial parameter values
--setPhysicsModelParameterRange cwww=-10,10: sets the parameter range to be scanned
-n: is added to the output name

# Results are saved in higgsCombineTestExample.MultiDimFit.mH120.root. To get 68% and 95% C.L. limits run e.g.
python buil1DInterval.py -10 10 higgsCombineExample.MultiDimFit.mH120.root cwww

# To get 2-Dimensional limits run e.g.
combine workspace_simfit.root -M MultiDimFit --floatOtherPOIs=0 --algo=grid --expectSignal=1 --points=1000 --redefineSignalPOIs cwww,ccw -P cwww -P ccw --freezeNuisances cb --setPhysicsModelParameters cwww=0,ccw=0,cb=0 --setPhysicsModelParameterRange cwww=-10,10:ccw=-20,20 --minimizerStrategy=2 -n Example

# To get the exact fit results for any point we need to run (e.g. for cwww=10)
combine workspace_simfit.root -M MaxLikelihoodFit --expectSignal=1 --freezeNuisances ccw,cb --setPhysicsModelParameters cwww=10,ccw=0,cb=0 --minimizerStrategy 2 --redefineSignalPOIs cwww --saveNormalizations --saveWithUncertainties --skipBOnlyFit -n Example
# The output is saved in mlfitExample.root containing a RooFitResult fit_s with all final parameter values as well as a RooArgSet norm_fit_s with the final normalizations.

# To plot the results we can use
python check_combine_results.py -n Example -c mu -P cwww:10
# which plots the mj and mlvj spectrum in all three regions while setting cwww to 10.

# If we want to get the background-only fit we have to freeze all aTGC-parameters and set a different POI, e.g.
combine workspace_simfit.root -M MaxLikelihoodFit --expectSignal=1 --freezeNuisances cwww,ccw,cb --setPhysicsModelParameters cwww=0,ccw=0,cb=0 --minimizerStrategy 2 --redefineSignalPOIs normvar_WJets_el --saveNormalizations --saveWithUncertainties --skipBOnlyFit -n BkgOnly

