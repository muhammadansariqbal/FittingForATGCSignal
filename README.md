SIGNAL PARAMETRIZATION
======================

```
This script creates the complete signal model and creates a datacard and a workspace per channel containing everything needed for the limit extraction.

Setup Instructions
==================

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
# -n: Read the input trees and create RooDataHists(-> faster access); Needed at the first run or when the input trees are changed.
# -c {channel}: Only run for {channel} (mu or el)
# -p: Make plots
# --savep: Save the plots
# -b: Run in batch mode
# --noatgcint: Set aTGC-interference terms to zero
# --printatgc: Print the coefficients of the signal model
# --atgc: Using different parametrization (Lagrangian approach instead of EFT)
# --binWidth: Use a different bin width than the standard, useful for Asimov data generation
# --cutoff: Specify mWV upper limit in GeV

# The workspaces for the different channels can now be combined with
text2workspace.py aC_WWWZ_simfit.txt -o workspace_simfit.root -P CombinedEWKAnalysis.CommonTools.ACModel:par1par2par3_TF3_shape_Model --PO channels=WWWZ_sig_el,WWWZ_sig_mu,WWWZ_sb_lo_el,WWWZ_sb_lo_mu,WWWZ_sb_hi_el,WWWZ_sb_hi_mu --PO poi=cwww,ccw,cb --PO range_cwww=-20,20 --PO range_ccw=-30,30 --PO range_cb=-75,75
# -o: Name of the created workspace
# -P: Name of the used model
# --PO channels= Names of the channels
# --PO poi= Names of the paramters of interest
# --PO range_= Set paramter range (does't work atm but has to be added to avoid error message)
# For vertex parametrization
text2workspace.py aC_WWWZ_simfit.txt -o workspace_simfit.root -P CombinedEWKAnalysis.CommonTools.ACModel:par1par2par3_TF3_shape_Model --PO channels=WWWZ_sig_el,WWWZ_sig_mu,WWWZ_sb_lo_el,WWWZ_sb_lo_mu,WWWZ_sb_hi_el,WWWZ_sb_hi_mu --PO poi=lZ,dg1z,dkz --PO range_lZ=-0.1,0.1 --PO range_dg1z=-0.1,0.1 --PO range_dkz=-0.1,0.1

# This creates the final workspace called workspace_simfit.root. To inspect this workspace in ROOT you have to load the combined limit .so 
.L ../../lib/slc6_amd64_gcc481/libHiggsAnalysisCombinedLimit.so.


Limit Calculation
=================
# We use the Higgs combine tool for doing the final data fits and calculating limits on aTGC parameters. Some of the arguments used by the combine tool are:
# -M: Likelihood method, we use MultiDimFit or MaxLikelihoodFit
# -t -1 or --expectSignal=1: -t -1 for Asimov data set generation
# --saveToys: The Asimov data set is saved in the toys directory as RooDataSet
# -n: This is added to the output name
# --points: The number of scanned parameter values
# --redefineSignalPOIs cwww -P cwww: Defines the paramter(s) of interest
# --freezeNuisances ccw,cb: Fixes the other paramters
# --setPhysicsModelParameters cwww=0,ccw=0,cb=0: Sets the initial parameter values
# --setPhysicsModelParameterRange cwww=-3.6,3.6: Sets the parameter range to be scanned
# --cminPreScan so that combine locates the correct global minimum

Fits for Single Point
---------------------
# To get the exact fit results for any point (e.g. cwww=3.6) we need to run
combine workspace_simfit.root -M MaxLikelihoodFit --expectSignal=1 --freezeNuisances ccw,cb --setPhysicsModelParameters cwww=3.6,ccw=0,cb=0 --minimizerStrategy 2 --cminPreScan --redefineSignalPOIs cwww --saveNormalizations --saveWithUncertainties --skipBOnlyFit -n _cwww_3.6
# The output is saved in mlfit_cwww_3.6.root containing a RooFitResult fit_s with all final parameter values as well as a RooArgSet norm_fit_s with the final normalizations.

# To get the results for all parameters zero
combine workspace_simfit.root -M MaxLikelihoodFit --expectSignal=1 --freezeNuisances ccw,cb --setPhysicsModelParameters cwww=0,ccw=0,cb=0 --minimizerStrategy 2 --cminPreScan --redefineSignalPOIs cwww --saveNormalizations --saveWithUncertainties --skipBOnlyFit -n AllZero

# We can also freeze all aTGC parameters and set a different POI
combine workspace_simfit.root -M MaxLikelihoodFit --expectSignal=1 --freezeNuisances cwww,ccw,cb --setPhysicsModelParameters cwww=0,ccw=0,cb=0 --minimizerStrategy 2 --cminPreScan --redefineSignalPOIs normvar_WJets_el --saveNormalizations --saveWithUncertainties --skipBOnlyFit -n BkgOnly

Asimov Data Set Generation
--------------------------
combine workspace_simfit.root -M MaxLikelihoodFit -t -1 --saveToys --freezeNuisances ccw,cb --setPhysicsModelParameters cwww=0,ccw=0,cb=0 --minimizerStrategy 2 --cminPreScan --redefineSignalPOIs cwww --saveNormalizations --saveWithUncertainties --skipBOnlyFit -n Asimov
# Add --toysFrequentist to generate Asimov data after getting optimal nuisance parameters values from a fit to data. Useful to get projected signal values.

Postfit Plots
-------------
python check_combine_result.py -n BkgOnly -c mu -P cwww:3.6
# -n: Addendum to the file name which contains the fit results
# -P: Set parameter value (default is -P cwww:0)
# -a: If we want to plot the Asimov data -a is the full file name which was saved from --saveToys (e.g. higgsCombineAsimov.MaxLikelihoodFit.mH120.123456.root)
# -r: Region

# For the new split plots, use the following for example
python check_combine_result_mJ.py -n AllZero -c el -P cwww:0 -r sig
python check_combine_result_mWV.py -n AllZero -c el -P cwww:0 -r sig


1-D Limits
----------
combine workspace_simfit.root -M MultiDimFit --floatOtherPOIs=0 --algo=grid --expectSignal=1 --points=1000 --redefineSignalPOIs cwww -P cwww --freezeNuisances ccw,cb --setPhysicsModelParameters cwww=0,ccw=0,cb=0 --setPhysicsModelParameterRange cwww=-3.6,3.6 --minimizerStrategy=2 --cminPreScan -n _cwww_3.6
combine workspace_simfit.root -M MultiDimFit --floatOtherPOIs=0 --algo=grid --expectSignal=1 --points=1000 --redefineSignalPOIs ccw -P ccw --freezeNuisances cwww,cb --setPhysicsModelParameters cwww=0,ccw=0,cb=0 --setPhysicsModelParameterRange ccw=-4.5,4.5 --minimizerStrategy=2 --cminPreScan -n _ccw_4.5
combine workspace_simfit.root -M MultiDimFit --floatOtherPOIs=0 --algo=grid --expectSignal=1 --points=1000 --redefineSignalPOIs cb -P cb --freezeNuisances cwww,ccw --setPhysicsModelParameters cwww=0,ccw=0,cb=0 --setPhysicsModelParameterRange cb=-20,20 --minimizerStrategy=2 --cminPreScan -n _cb_20

# For vertex parametrization
combine workspace_simfit.root -M MultiDimFit --floatOtherPOIs=0 --algo=grid --expectSignal=1 --points=1000 --redefineSignalPOIs lZ -P lZ --freezeNuisances dg1z,dkz --setPhysicsModelParameters lZ=0,dg1z=0,dkz=0 --setPhysicsModelParameterRange lZ=-0.014,0.014 --minimizerStrategy=2 --cminPreScan -n _lZ_0.014
combine workspace_simfit.root -M MultiDimFit --floatOtherPOIs=0 --algo=grid --expectSignal=1 --points=1000 --redefineSignalPOIs dg1z -P dg1z --freezeNuisances lZ,dkz --setPhysicsModelParameters lZ=0,dg1z=0,dkz=0 --setPhysicsModelParameterRange dg1z=-0.018,0.018 --minimizerStrategy=2 --cminPreScan -n _dg1z_0.018
combine workspace_simfit.root -M MultiDimFit --floatOtherPOIs=0 --algo=grid --expectSignal=1 --points=1000 --redefineSignalPOIs dkz -P dkz --freezeNuisances lZ,dg1z --setPhysicsModelParameters lZ=0,dg1z=0,dkz=0 --setPhysicsModelParameterRange dkz=-0.02,0.02 --minimizerStrategy=2 --cminPreScan -n _dkz_0.02

2-D Limits
----------
combine workspace_simfit.root -M MultiDimFit --floatOtherPOIs=0 --algo=grid --expectSignal=1 --points=1000 --redefineSignalPOIs cwww,ccw -P cwww -P ccw --freezeNuisances cb --setPhysicsModelParameters cwww=0,ccw=0,cb=0 --setPhysicsModelParameterRange cwww=-3.6,3.6:ccw=-4.5,4.5 --minimizerStrategy=2 --cminPreScan -n _cwww_3.6_ccw_4.5
combine workspace_simfit.root -M MultiDimFit --floatOtherPOIs=0 --algo=grid --expectSignal=1 --points=1000 --redefineSignalPOIs cwww,cb -P cwww -P cb --freezeNuisances ccw --setPhysicsModelParameters cwww=0,ccw=0,cb=0 --setPhysicsModelParameterRange cwww=-3.6,3.6:cb=-20,20 --minimizerStrategy=2 --cminPreScan -n _cwww_3.6_cb_20
combine workspace_simfit.root -M MultiDimFit --floatOtherPOIs=0 --algo=grid --expectSignal=1 --points=1000 --redefineSignalPOIs ccw,cb -P ccw -P cb --freezeNuisances cwww --setPhysicsModelParameters cwww=0,ccw=0,cb=0 --setPhysicsModelParameterRange ccw=-4.5,4.5:cb=-20,20 --minimizerStrategy=2 --cminPreScan -n _ccw_4.5_cb_20

# For vertex parametrization
combine workspace_simfit.root -M MultiDimFit --floatOtherPOIs=0 --algo=grid --expectSignal=1 --points=1000 --redefineSignalPOIs lZ,dg1z -P lZ -P dg1z --freezeNuisances dkz --setPhysicsModelParameters lZ=0,dg1z=0,dkz=0 --setPhysicsModelParameterRange lZ=-0.014,0.014:dg1z=-0.018,0.018 --minimizerStrategy=2 --cminPreScan -n _lZ_0.014_dg1z_0.018
combine workspace_simfit.root -M MultiDimFit --floatOtherPOIs=0 --algo=grid --expectSignal=1 --points=1000 --redefineSignalPOIs lZ,dkz -P lZ -P dkz --freezeNuisances dg1z --setPhysicsModelParameters lZ=0,dg1z=0,dkz=0 --setPhysicsModelParameterRange lZ=-0.014,0.014:dkz=-0.02,0.02 --minimizerStrategy=2 --cminPreScan -n _lZ_0.014_dkz_0.02
combine workspace_simfit.root -M MultiDimFit --floatOtherPOIs=0 --algo=grid --expectSignal=1 --points=1000 --redefineSignalPOIs dg1z,dkz -P dg1z -P dkz --freezeNuisances dkz --setPhysicsModelParameters lZ=0,dg1z=0,dkz=0 --setPhysicsModelParameterRange dg1z=-0.018,0.018:dkz=-0.02,0.02 --minimizerStrategy=2 --cminPreScan -n _dg1z_0.018_dkz_0.02

Get 68% and 95% Confidence Intervals
------------------------------------
python build1DInterval.py -3.6 3.6 higgsCombine_cwww_3.6.MultiDimFit.mH120.root cwww
python build1DInterval.py -4.5 4.5 higgsCombine_ccw_4.5.MultiDimFit.mH120.root ccw
python build1DInterval.py -20 20 higgsCombine_cb_20.MultiDimFit.mH120.root cb

python build1DInterval.py -0.014 0.014 higgsCombine_lZ_0.014.MultiDimFit.mH120.root lZ
python build1DInterval.py -0.018 0.018 higgsCombine_dg1z_0.018.MultiDimFit.mH120.root dg1z
python build1DInterval.py -0.02 0.02 higgsCombine_dkz_0.02.MultiDimFit.mH120.root dkz

Using screen sessions for long limit calculations
-------------------------------------------------

hostname
pagsh.krb -c 'kinit && screen'
screen -ls
screen -r

