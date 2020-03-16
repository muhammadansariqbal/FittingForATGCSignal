from ROOT import  *
from array import array
from optparse import OptionParser
import ROOT
import math as math
import os
import CMS_lumi, tdrstyle
import numpy as np

gSystem.Load("%s/lib/slc6_amd64_gcc481/libHiggsAnalysisCombinedLimit.so"%os.environ["CMSSW_BASE"])

gSystem.Load("PDFs/PdfDiagonalizer_cc.so")
gSystem.Load("PDFs/Util_cxx.so")
gSystem.Load("PDFs/hyperg_2F1_c.so")
gSystem.Load("PDFs/HWWLVJRooPdfs_cxx.so")
from ROOT import gROOT,draw_error_band
from ROOT import RooPoly3Pdf, RooChiSqPdf, RooErfExpPdf, RooErfPowExpPdf, RooErfPowPdf, RooErfPow2Pdf, RooExpNPdf, RooAlpha4ExpNPdf, RooExpTailPdf, RooAlpha4ExpTailPdf, Roo2ExpPdf, draw_error_band_extendPdf

gStyle.SetOptStat(0)
gStyle.SetOptTitle(0)
parser	= OptionParser()
parser.add_option('--POI',dest='POI',help='parameter of interest')
parser.add_option('--pval',dest='pval',help='value of parameter')
parser.add_option('-c',action='store_true',dest='close',default=False)
parser.add_option('--binWidths',dest='binWidths',default='0',help='plot best fit points with different bin widths for bias testing')
parser.add_option('--year',dest='year',default='16')
(options,args) = parser.parse_args()

POI		= []
pval		= []
binWidths	= []
for i in options.POI.split(','):
	POI.append(i)
for j in options.pval.split(','):
	pval.append(j)
for k in options.binWidths.split(','):
	binWidths.append(k)
year		= options.year
par_latex	= {'cwww' : 'c_{WWW} / #Lambda^{2} (TeV^{-2})', 'ccw' : 'c_{W} / #Lambda^{2} (TeV^{-2})', 'cb' : 'c_{B} / #Lambda^{2} (TeV^{-2})', 'lZ' : '#lambda_{Z}', 'dg1z' : '#Deltag_{1}^{Z}', 'dkz' : '#Delta#kappa_{Z}'}
par_noUnits	= {'cwww' : 'c_{WWW} / #Lambda^{2}', 'ccw' : 'c_{W} / #Lambda^{2}', 'cb' : 'c_{B} / #Lambda^{2}', 'lZ' : '#lambda_{Z}', 'dg1z' : '#Deltag_{1}^{Z}', 'dkz' : '#Delta#kappa_{Z}'}

def plots():
	pathExpected    = './ResultsExpected'+year+'/'
	pathObserved    = './ResultsObserved'+year+'/'

	# Make the expected TGraph	
	wsNameExp	= 'higgsCombine_%s_%s_%s_%s.MultiDimFit.mH120.root'%(POI[0],pval[0],POI[1],pval[1])
	print 'Reading expected '+wsNameExp
	fileInATGCExp	= TFile.Open(pathExpected+wsNameExp)
	treeExp		= fileInATGCExp.Get('limit')
	NEntriesExp	= treeExp.GetEntries()

	treeExp.GetEntry(1)
	par1	= POI[0]
	par2	= POI[1]
	xExp	= []
	yExp	= []
	zExp	= []

	for i in range(NEntriesExp-1):
		if i%1000==0:
                        print i
		treeExp.GetEntry(i+1)
		if 2*treeExp.deltaNLL < 599:
			xExp.append(treeExp.GetLeaf(par1).GetValue())
			yExp.append(treeExp.GetLeaf(par2).GetValue())
			zExp.append(2*treeExp.deltaNLL)

	graphExp	= TGraph2D(len(xExp),array('d',xExp),array('d',yExp),array('d',zExp))

	# Make the observed TGraph      
        wsNameObs       = 'higgsCombine_%s_%s_%s_%s.MultiDimFit.mH120.root'%(POI[0],pval[0],POI[1],pval[1])
        print 'Reading observed '+wsNameObs
        fileInATGCObs   = TFile.Open(pathObserved+wsNameObs)
        treeObs         = fileInATGCObs.Get('limit')
        NEntriesObs     = treeObs.GetEntries()

        treeObs.GetEntry(1)
        xObs    = []
        yObs    = []
        zObs    = []

        for i in range(NEntriesObs-1):
                if i%1000==0:
                        print i
                treeObs.GetEntry(i+1)
                if 2*treeObs.deltaNLL < 599:
                        xObs.append(treeObs.GetLeaf(par1).GetValue())
                        yObs.append(treeObs.GetLeaf(par2).GetValue())
                        zObs.append(2*treeObs.deltaNLL)

        graphObs        = TGraph2D(len(xObs),array('d',xObs),array('d',yObs),array('d',zObs))

	#Best fit point
	bestFitXBin	= ROOT.Long(0)
	bestFitYBin	= ROOT.Long(0)
	minDNLL		= ROOT.Long(0)
	graphObs.GetHistogram().GetMinimumBin(bestFitXBin,bestFitYBin,minDNLL)
	#bestFitX	= graphExp.GetHistogram().GetXaxis().GetBinCenter(bestFitXBin)
	#bestFitY	= graphExp.GetHistogram().GetYaxis().GetBinCenter(bestFitYBin)

	# Get best fit directly from tree instead of Delaunay histogram
	treeObs.GetEntry(0)
	bestFitX	= treeObs.GetLeaf(par1).GetValue()
	bestFitY        = treeObs.GetLeaf(par2).GetValue()

	c1              = TCanvas('c1','c1',800,750)
	c1.SetLeftMargin(0.13)
	c1.SetRightMargin(0.04)
	c1.SetBottomMargin(0.13)
	c1.SetTopMargin(0.1)

	# Draw temporary contours and extract the TGraphs from them
	contourLevels	= np.array([2.3, 5.99, 9.21])
	graphExp.GetHistogram().SetContour(3,contourLevels)
	graphExp.Draw("CONT LIST")
	c1.Update()
        c1.SetGrid()

	contourExp99    = TGraph(TObjArray(ROOT.gROOT.GetListOfSpecials().FindObject("contours"))[2].First())
	contourExp95    = TGraph(TObjArray(ROOT.gROOT.GetListOfSpecials().FindObject("contours"))[1].First())
	contourExp68    = TGraph(TObjArray(ROOT.gROOT.GetListOfSpecials().FindObject("contours"))[0].First())

	# Do the same for observed
	contourLevels   = np.array([5.99])
        graphObs.GetHistogram().SetContour(1,contourLevels)
        graphObs.Draw("CONT SAME LIST")
	c1.Update()
        c1.SetGrid()
	
	contourObs95    = TGraph(TObjArray(ROOT.gROOT.GetListOfSpecials().FindObject("contours"))[0].First())

	# 99% Expected =========================================================================================

	contourExp99.SetLineStyle(10)
	contourExp99.SetLineColor(kRed)
        contourExp99.SetLineWidth(2)

	contourExp99.GetXaxis().SetTitle(par_latex[par1])
        contourExp99.GetYaxis().SetTitle(par_latex[par2])
	contourExp99.GetXaxis().SetTitleSize(0.068)
	contourExp99.GetYaxis().SetTitleSize(0.068)
	contourExp99.GetXaxis().SetTitleOffset(0.8)
	contourExp99.GetYaxis().SetTitleOffset(0.85)
	contourExp99.GetXaxis().SetLabelSize(0.06)
	contourExp99.GetYaxis().SetLabelSize(0.06)
	contourExp99.GetYaxis().SetLabelOffset(0.015)
	contourExp99.GetXaxis().SetNdivisions(505)
        contourExp99.GetYaxis().SetNdivisions(505)
	if (par2=='dg1z' or par2=='dkz'):
		contourExp99.GetYaxis().SetTitle("")
	contourExp99.Draw('AC')
	c1.Update()

	contourExp99.GetYaxis().SetRangeUser(c1.GetFrame().GetY1(),1.5*c1.GetFrame().GetY2())
	contourExp99.Draw('AC')

	# 95% Expected =========================================================================================

        contourExp95.SetLineStyle(9)
        contourExp95.SetLineColor(kGreen+2)
        contourExp95.SetLineWidth(2)

        contourExp95.Draw('C SAME')

	# 68% Expected =========================================================================================

        contourExp68.SetLineStyle(2)
        contourExp68.SetLineColor(kBlue)
        contourExp68.SetLineWidth(2)

        contourExp68.Draw('C SAME')
	
	# 95% Observed =========================================================================================

        contourObs95.SetLineStyle(1)
        contourObs95.SetLineColor(kBlack)
        contourObs95.SetLineWidth(4)

        contourObs95.Draw('C SAME')

	# Points ===============================================================================================

	SMPoint	= TGraph(1)
	SMPoint.SetPoint(1,0,0)
	SMPoint.SetMarkerStyle(21)
	SMPoint.SetMarkerSize(2)
	SMPoint.Draw('P SAME')

	if float(binWidths[0])==0:
		bestFitPoint = TGraph(1)
        	bestFitPoint.SetPoint(0,bestFitX,bestFitY)
        	bestFitPoint.SetMarkerStyle(34)
		bestFitPoint.SetMarkerSize(3)
        	bestFitPoint.Draw('P SAME')
	else:
		bestFitPoints	= []
		colors	= []
		for binWid in binWidths:
			wsNameExp	= 'higgsCombine_%s_%s_%s_%s_binWidth%s.MultiDimFit.mH120.root'%(POI[0],pval[0],POI[1],pval[1],binWid)
			print 'Reading '+wsNameExp
			fileInBias	= TFile.Open(path+wsNameExp)
			treeInBias	= fileInBias.Get('limit')
			treeInBias.GetEntry(0)
			bestFitPoint = TGraph(1)
                	bestFitPoint.SetPoint(0,treeInBias.GetLeaf(par1).GetValue(),treeInBias.GetLeaf(par2).GetValue())
                	bestFitPoint.SetMarkerStyle(3)
			colorGray=(100-float(binWid))/150.0
			grayRGB=TColor(2000+int(binWid),colorGray,colorGray,colorGray)
			bestFitPoint.SetMarkerColor(2000+int(binWid))
                	bestFitPoint.Draw('P SAME')
			bestFitPoints.append(bestFitPoint)
			colors.append(grayRGB)
			fileInBias.Close()

	# ======================================================================================================

	if (par2=='dg1z' or par2=='dkz'):
		cover = TPaveText(0.005,0.78,0.12,0.92,"b1NDC")
		cover.SetFillColor(kWhite)
		cover.Draw("SAME")
		axisTitle = TLatex(0.01,0.85, par_latex[par2])
		axisTitle.SetNDC()
		axisTitle.SetTextAlign(23)
		axisTitle.SetTextAngle(90)
		axisTitle.SetTextSize(0.068)
		axisTitle.Draw("SAME")

	if year=='17':
		CMS_lumi.lumi_13TeV = "41.5 fb^{-1}"
	elif year=='18':
		CMS_lumi.lumi_13TeV = "59.9 fb^{-1}"
	elif year=='Run2':
		CMS_lumi.lumi_13TeV = "137 fb^{-1}"
	CMS_lumi.cmsTextSize	= 0.7
        CMS_lumi.relPosY        = -0.08
	CMS_lumi.lumiTextSize	= 0.525
	CMS_lumi.lumiTextOffset	= 0.1
	CMS_lumi.CMS_lumi(c1,4,11)

	leg=TLegend(0.14,0.73,0.95,0.89)
	leg.SetFillColor(kWhite)
        leg.SetBorderSize(0)
        leg.SetLineColor(0)
        leg.SetLineWidth(0)
        leg.SetLineStyle(0)
	leg.SetTextFont(42)
	leg.SetTextSize(0.043)
	leg.SetNColumns(2)
	leg.AddEntry(contourExp68,"Expected 68% CL","L")
	leg.AddEntry(contourExp95,"Expected 95% CL","L")
	leg.AddEntry(contourExp99,"Expected 99% CL","L")
	leg.AddEntry(contourObs95,"Observed 95% CL","L")
	leg.AddEntry(SMPoint,"SM point","P")
	leg.AddEntry(bestFitPoint,"Best-fit","P")
	leg.Draw("SAME")

	c1.Update()
	c1.SaveAs("limit2D_%s_%s_%s.pdf"%(year,par1,par2))

plots()
