from ROOT import  *
from array import array
from optparse import OptionParser
import math as math
import os
import CMS_lumi, tdrstyle

gSystem.Load("%s/lib/slc6_amd64_gcc481/libHiggsAnalysisCombinedLimit.so"%os.environ["CMSSW_BASE"])

gSystem.Load("PDFs/PdfDiagonalizer_cc.so")
gSystem.Load("PDFs/Util_cxx.so")
gSystem.Load("PDFs/hyperg_2F1_c.so")
gSystem.Load("PDFs/HWWLVJRooPdfs_cxx.so")
from ROOT import draw_error_band
from ROOT import RooPoly3Pdf, RooChiSqPdf, RooErfExpPdf, RooErfPowExpPdf, RooErfPowPdf, RooErfPow2Pdf, RooExpNPdf, RooAlpha4ExpNPdf, RooExpTailPdf, RooAlpha4ExpTailPdf, Roo2ExpPdf, draw_error_band_extendPdf

gStyle.SetOptStat(0)
gStyle.SetOptTitle(0)
parser	= OptionParser()
parser.add_option('--POI',dest='POI',help='parameter of interest')
parser.add_option('--pval',dest='pval',help='value of parameter')
parser.add_option('-c',action='store_true',dest='close',default=False)
parser.add_option('--year',dest='year',default='16')
(options,args) = parser.parse_args()

POI		= options.POI
pval		= options.pval
year		= options.year
par_latex	= {'cwww' : 'c_{WWW} / #Lambda^{2} (TeV^{-2})', 'ccw' : 'c_{W} / #Lambda^{2} (TeV^{-2})', 'cb' : 'c_{B} / #Lambda^{2} (TeV^{-2})', 'lZ' : '#lambda_{Z}', 'dg1z' : '#Deltag_{1}^{Z}', 'dkz' : '#Delta#kappa_{Z}'}
par_noUnits	= {'cwww' : 'c_{WWW} / #Lambda^{2}', 'ccw' : 'c_{W} / #Lambda^{2}', 'cb' : 'c_{B} / #Lambda^{2}', 'lZ' : '#lambda_{Z}', 'dg1z' : '#Deltag_{1}^{Z}', 'dkz' : '#Delta#kappa_{Z}'}

def plots():
	pathExpected	= './ResultsExpected'+year+'/'
	pathObserved	= './ResultsObserved'+year+'/'
	
	# Make the expected TGraph
	wsNameExp	= 'higgsCombine_%s_%s.MultiDimFit.mH120.root'%(POI,pval)
	print 'Reading expected '+wsNameExp
	fileInATGCExp	= TFile.Open(pathExpected+wsNameExp)
	treeExp		= fileInATGCExp.Get('limit')
	NEntriesExp	= treeExp.GetEntries()

	treeExp.GetEntry(1)
	par	= POI
	xExp	= []
	yExp	= []

	for i in range(NEntriesExp-1):
		if i%1000==0:
                        print i
		treeExp.GetEntry(i+1)
		if 2*treeExp.deltaNLL < 8:
			xExp.append(treeExp.GetLeaf(par).GetValue())
			yExp.append(2*treeExp.deltaNLL)

	graphExp	= TGraph(len(xExp),array('d',xExp),array('d',yExp))
	c1              = TCanvas()
	#c1.SetRightMargin(0.2)
	if yExp[0]<yExp[len(xExp)-1]:
                yMax=yExp[0]
        else:
                yMax=yExp[len(xExp)-1]
        graphExp.GetXaxis().SetTitle(par_latex[par])
	graphExp.GetYaxis().SetTitle('2#DeltaNLL')
        graphExp.GetYaxis().SetTitleSize(0.05)
        graphExp.GetYaxis().SetTitleOffset(0.75)
        graphExp.GetYaxis().SetRangeUser(0,yMax)
        graphExp.GetXaxis().SetTitleSize(0.05)
        graphExp.GetXaxis().SetTitleOffset(0.75)
        graphExp.SetLineStyle(1)
        graphExp.SetLineWidth(2)
        graphExp.SetLineColor(kGreen+2)

	# Make the observed TGraph
	wsNameObs       = 'higgsCombine_%s_%s.MultiDimFit.mH120.root'%(POI,pval)
        print 'Reading observed '+wsNameExp
        fileInATGCObs   = TFile.Open(pathObserved+wsNameExp)
        treeObs         = fileInATGCObs.Get('limit')
        NEntriesObs     = treeObs.GetEntries()

        treeObs.GetEntry(1)
        xObs    = []
        yObs    = []

        for i in range(NEntriesObs-1):
                if i%1000==0:
                        print i
                treeObs.GetEntry(i+1)
                if 2*treeObs.deltaNLL < 8:
                        xObs.append(treeObs.GetLeaf(par).GetValue())
                        yObs.append(2*treeObs.deltaNLL)

        graphObs        = TGraph(len(xObs),array('d',xObs),array('d',yObs))
	graphObs.SetLineStyle(1)
        graphObs.SetLineWidth(2)
        graphObs.SetLineColor(kBlack)

	# 95% Expected =========================================================================================

        for i in range((NEntriesExp-1)/2):
                j=i+1
                treeExp.GetEntry(j)
                if 2*treeExp.deltaNLL>3.84 and i<NEntriesExp/2-1:
                        continue
                else:
                        limLoExp95 = treeExp.GetLeaf(par).GetValue()
                        break
        for i in range((NEntriesExp-1)/2):
                j=i+ (NEntriesExp-1)/2 +1
                treeExp.GetEntry(j)
                if 2*treeExp.deltaNLL<3.84 and i<NEntriesExp/2-1:
                        continue
                else:
                        limHiExp95 = treeExp.GetLeaf(par).GetValue()
                        break

        line4           = TF1('line4','3.84', graphExp.GetXaxis().GetXmin(),graphExp.GetXaxis().GetXmax())
        line4.SetLineStyle(7)
        line4.SetLineColor(kBlack)
        line4.SetLineWidth(1)

        lineLoExp95 = TLine(limLoExp95,0,limLoExp95,yMax)
        lineLoExp95.SetLineStyle(7)
        #lineLoExp95.SetLineWidth(4)
        lineLoExp95.SetLineColor(kGreen+2)

        lineHiExp95 = TLine(limHiExp95,0,limHiExp95,yMax)
        lineHiExp95.SetLineStyle(7)
        #lineHiExp95.SetLineWidth
        lineHiExp95.SetLineColor(kGreen+2)

	# 95% Observed =========================================================================================

        for i in range((NEntriesObs-1)/2):
                j=i+1
                treeObs.GetEntry(j)
                if 2*treeObs.deltaNLL>3.84 and i<NEntriesObs/2-1:
                        continue
                else:
                        limLoObs95 = treeObs.GetLeaf(par).GetValue()
                        break
        for i in range((NEntriesObs-1)/2):
                j=i+ (NEntriesObs-1)/2 +1
                treeObs.GetEntry(j)
                if 2*treeObs.deltaNLL<3.84 and i<NEntriesObs/2-1:
                        continue
                else:
                        limHiObs95 = treeObs.GetLeaf(par).GetValue()
                        break

        lineLoObs95 = TLine(limLoObs95,0,limLoObs95,yMax)
        lineLoObs95.SetLineStyle(7)
        #lineLoObs95.SetLineWidth(4)
        lineLoObs95.SetLineColor(kBlack)

        lineHiObs95 = TLine(limHiObs95,0,limHiObs95,yMax)
        lineHiObs95.SetLineStyle(7)
        #lineHiObs95.SetLineWidth
        lineHiObs95.SetLineColor(kBlack)

	# 68% Expected =========================================================================================

	for i in range((NEntriesExp-1)/2):
                j=i+1
                treeExp.GetEntry(j)
                if 2*treeExp.deltaNLL>1.00 and i<NEntriesExp/2-1:
                        continue
                else:
                        limLoExp68 = treeExp.GetLeaf(par).GetValue()
                        break
        for i in range((NEntriesExp-1)/2):
                j=i+ (NEntriesExp-1)/2 +1
                treeExp.GetEntry(j)
                if 2*treeExp.deltaNLL<1.00 and i<NEntriesExp/2-1:
                        continue
                else:
                        limHiExp68 = treeExp.GetLeaf(par).GetValue()
                        break

        line1           = TF1('line1','1.00', graphExp.GetXaxis().GetXmin(),graphExp.GetXaxis().GetXmax())
        line1.SetLineStyle(7)
        line1.SetLineColor(kBlue)
        line1.SetLineWidth(1)

        lineLoExp68 = TLine(limLoExp68,0,limLoExp68,yMax)
        lineLoExp68.SetLineStyle(7)
        #lineLoExp68.SetLineWidth(4)
        lineLoExp68.SetLineColor(kBlue)

        lineHiExp68 = TLine(limHiExp68,0,limHiExp68,yMax)
        lineHiExp68.SetLineStyle(7)
        #lineHiExp68.SetLineWidth
        lineHiExp68.SetLineColor(kBlue)

	# 99% Expected =========================================================================================

        for i in range((NEntriesExp-1)/2):
                j=i+1
                treeExp.GetEntry(j)
                if 2*treeExp.deltaNLL>6.63 and i<NEntriesExp/2-1:
                        continue
                else:
                        limLoExp99 = treeExp.GetLeaf(par).GetValue()
                        break
        for i in range((NEntriesExp-1)/2):
                j=i+ (NEntriesExp-1)/2 +1
                treeExp.GetEntry(j)
                if 2*treeExp.deltaNLL<6.63 and i<NEntriesExp/2-1:
                        continue
                else:
                        limHiExp99 = treeExp.GetLeaf(par).GetValue()
                        break

        line6           = TF1('line6','6.63', graphExp.GetXaxis().GetXmin(),graphExp.GetXaxis().GetXmax())
        line6.SetLineStyle(7)
        line6.SetLineColor(kRed)
        line6.SetLineWidth(1)

        lineLoExp99 = TLine(limLoExp99,0,limLoExp99,yMax)
        lineLoExp99.SetLineStyle(7)
        #lineLoExp99.SetLineWidth(4)
        lineLoExp99.SetLineColor(kRed)

        lineHiExp99 = TLine(limHiExp99,0,limHiExp99,yMax)
        lineHiExp99.SetLineStyle(7)
        #lineHiExp99.SetLineWidth
        lineHiExp99.SetLineColor(kRed)

	# Shaded area between 68% and 99% ======================================================================

	boxLow	= TBox(limLoExp99,0,limLoExp68,yMax)
	boxLow.SetFillStyle(3002)
	boxLow.SetFillColor(kGray)
	boxLow.SetLineColor(kGray)

	boxHigh	= TBox(limHiExp99,0,limHiExp68,yMax)
        boxHigh.SetFillStyle(3002)
        boxHigh.SetFillColor(kGray)
	boxHigh.SetLineColor(kGray)

	# Draw everything in the order needed ==================================================================

	graphExp.Draw('AC')

	boxLow.Draw('SAME')
	boxHigh.Draw('SAME')

	graphExp.Draw('C SAME')
	graphObs.Draw('C SAME')

	line4.Draw('SAME')
        lineLoExp95.Draw('SAME')
        lineHiExp95.Draw('SAME')

	lineLoObs95.Draw('SAME')
        lineHiObs95.Draw('SAME')

	#line1.Draw('SAME')
        lineLoExp68.Draw('SAME')
        lineHiExp68.Draw('SAME')

	#line6.Draw('SAME')
        lineLoExp99.Draw('SAME')
        lineHiExp99.Draw('SAME')

	# ======================================================================================================

	if year=='17':
		CMS_lumi.lumi_13TeV = "41.5 fb^{-1}"
	elif year=='18':
		CMS_lumi.lumi_13TeV = "59.9 fb^{-1}"
	elif year=='Run2':
		CMS_lumi.lumi_13TeV = "137 fb^{-1}"
	CMS_lumi.lumiTextOffset=0.1
	CMS_lumi.relPosY    = -0.05
	CMS_lumi.writeExtraText=True
	CMS_lumi.extraText = "Private work"
	CMS_lumi.relExtraDX = 0.08
	CMS_lumi.relExtraDY = 0.25
	CMS_lumi.CMS_lumi(c1,4,11)

	leg=TLegend(0.35,0.625,0.65,0.85)
	leg.SetFillColor(0)
        leg.SetFillStyle(0)
        leg.SetBorderSize(0)
        leg.SetLineColor(0)
        leg.SetLineWidth(0)
        leg.SetLineStyle(0)
	leg.SetTextFont(42)
	leg.SetHeader(par_noUnits[par])
	leg.AddEntry(graphExp,"Expected","L")
	leg.AddEntry(graphObs,"Observed","L")
	leg.AddEntry(boxHigh,"Expected #pm 1#sigma","F")
	leg.Draw("SAME")

	c1.Update()
	c1.SaveAs("limit1D_%s_%s.pdf"%(year,POI))

	errorExp = float(xExp[1]) - float(xExp[0])
	errorObs = float(xObs[1]) - float(xObs[0])

	print '95%% expected C.L. limit on %s: [%s,%s] +- %s'%(par,round(limLoExp95,2),round(limHiExp95,2),round(errorExp,2))
	print '95%% observed C.L. limit on %s: [%s,%s] +- %s'%(par,round(limLoObs95,2),round(limHiObs95,2),round(errorObs,2))

plots()
