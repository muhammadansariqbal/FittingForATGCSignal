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
(options,args) = parser.parse_args()

POI		= []
pval		= []
for i in options.POI.split(','):
	POI.append(i)
for j in options.pval.split(','):
	pval.append(j)
par_latex	= {'cwww' : 'c_{WWW} / #Lambda^{2} (TeV^{-2})', 'ccw' : 'c_{W} / #Lambda^{2} (TeV^{-2})', 'cb' : 'c_{B} / #Lambda^{2} (TeV^{-2})', 'lZ' : '#lambda_{Z}', 'dg1z' : '#Deltag_{1}^{Z}', 'dkz' : '#Delta#kappa_{Z}'}
par_noUnits	= {'cwww' : 'c_{WWW} / #Lambda^{2}', 'ccw' : 'c_{W} / #Lambda^{2}', 'cb' : 'c_{B} / #Lambda^{2}', 'lZ' : '#lambda_{Z}', 'dg1z' : '#Deltag_{1}^{Z}', 'dkz' : '#Delta#kappa_{Z}'}

def plots():
	path		= './'
	
	wsNameExp	= 'higgsCombine_%s_%s_%s_%s.MultiDimFit.mH120.root'%(POI[0],pval[0],POI[1],pval[1])
	print 'Reading '+wsNameExp
	fileInATGCExp	= TFile.Open(path+wsNameExp)
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
	bestFitXBin	= ROOT.Long(0)
	bestFitYBin	= ROOT.Long(0)
	minDNLL		= ROOT.Long(0)
	graphExp.GetHistogram().GetMinimumBin(bestFitXBin,bestFitYBin,minDNLL)
	#bestFitX	= graphExp.GetHistogram().GetXaxis().GetBinCenter(bestFitXBin)
	#bestFitY	= graphExp.GetHistogram().GetYaxis().GetBinCenter(bestFitYBin)

	# Get best fit directly from tree instead of Delaunay histogram
	treeExp.GetEntry(0)
	bestFitX	= treeExp.GetLeaf(par1).GetValue()
	bestFitY        = treeExp.GetLeaf(par2).GetValue()

	c1              = TCanvas('c1','c1',800,750)
	contourLevels	= np.array([2.3, 5.99, 9.21])
	graphExp.GetHistogram().SetContour(3,contourLevels)
	graphExp.Draw("CONT LIST")
	c1.Update()
	c1.SetGrid()

	# 99% Expected =========================================================================================

	contourExp99    = TObjArray(ROOT.gROOT.GetListOfSpecials().FindObject("contours"))[2].First()

	contourExp99.SetLineStyle(9)
	contourExp99.SetLineColor(kRed)
        contourExp99.SetLineWidth(1)

	contourExp99.GetXaxis().SetTitle(par_latex[par1])
        contourExp99.GetYaxis().SetTitle(par_latex[par2])
	contourExp99.GetXaxis().SetTitleSize(0.05)
	contourExp99.GetYaxis().SetTitleSize(0.05)
	contourExp99.GetXaxis().SetTitleOffset(0.75)
	contourExp99.GetYaxis().SetTitleOffset(0.9)
	contourExp99.GetXaxis().SetNdivisions(505)
        contourExp99.GetYaxis().SetNdivisions(505)
	contourExp99.Draw('AC')
	c1.Update()

	contourExp99.GetYaxis().SetRangeUser(c1.GetFrame().GetY1(),1.4*c1.GetFrame().GetY2())
	contourExp99.Draw('AC')

	# 95% Expected =========================================================================================

        contourExp95    = TObjArray(ROOT.gROOT.GetListOfSpecials().FindObject("contours"))[1].First()

        contourExp95.SetLineStyle(9)
        contourExp95.SetLineColor(kGreen+2)
        contourExp95.SetLineWidth(1)

        contourExp95.Draw('C SAME')

	# 68% Expected =========================================================================================

        contourExp68    = TObjArray(ROOT.gROOT.GetListOfSpecials().FindObject("contours"))[0].First()

        contourExp68.SetLineStyle(9)
        contourExp68.SetLineColor(kBlue)
        contourExp68.SetLineWidth(1)

        contourExp68.Draw('C SAME')
	
	# Points ===============================================================================================

	SMPoint	= TGraph(1)
	SMPoint.SetPoint(1,0,0)
	SMPoint.SetMarkerStyle(21)
	SMPoint.Draw('P SAME')

	bestFitPoint = TGraph(1)
        bestFitPoint.SetPoint(1,bestFitX,bestFitY)
        bestFitPoint.SetMarkerStyle(3)
        bestFitPoint.Draw('P SAME')

	# ======================================================================================================

	CMS_lumi.lumiTextOffset	= 0.1
	CMS_lumi.relPosY	= -0.05
	CMS_lumi.relExtraDY	= 0.24
	CMS_lumi.CMS_lumi(c1,4,11)

	leg=TLegend(0.125,0.73,0.875,0.89)
	leg.SetFillColor(kWhite)
        leg.SetBorderSize(0)
        leg.SetLineColor(0)
        leg.SetLineWidth(0)
        leg.SetLineStyle(0)
	leg.SetTextFont(42)
	leg.SetNColumns(2)
	leg.AddEntry(contourExp68,"Expected 68% C.L.","L")
	leg.AddEntry(contourExp95,"Expected 95% C.L.","L")
	leg.AddEntry(contourExp99,"Expected 99% C.L.","L")
	leg.AddEntry(SMPoint,"SM point","P")
	leg.AddEntry(bestFitPoint,"Asimov best fit","P")
	leg.Draw("SAME")

	c1.Update()
	c1.SaveAs("limit2D_%s_%s.pdf"%(par1,par2))
	raw_input('<>')

plots()
