from ROOT import  *
from array import array
from optparse import OptionParser
import math as math
import os

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
(options,args) = parser.parse_args()

POI	= options.POI
pval	= options.pval

def plots():
	path		= './'
	wsname		= 'higgsCombine_%s_%s.MultiDimFit.mH120.root'%(POI,pval)
	print wsname
	fileInATGC	= TFile.Open(path+wsname)
	tree		= fileInATGC.Get('limit')
	NEntries	= tree.GetEntries()

	tree.GetEntry(1)
	if 'cwww' in wsname:
		par		= 'cwww'
		binlo		= tree.cwww
		tree.GetEntry(NEntries-1)
		binhi		= tree.cwww
	elif 'ccw' in wsname:
		par		= 'ccw'
		binlo		= tree.ccw
		tree.GetEntry(NEntries-1)
		binhi		= tree.ccw
	elif 'cb' in wsname:
		par		= 'cb'
		binlo		= tree.cb
		tree.GetEntry(NEntries-1)
		binhi		= tree.cb

	hist		= TH1F('hist','hist',NEntries, binlo-(0.5*(binhi-binlo)/(NEntries-1)), binhi+(0.5*(binhi-binlo)/(NEntries-1)))
	
	for i in range(NEntries-1):
		if i%1000==0:
			print i
		tree.GetEntry(i+1)
		if 2*tree.deltaNLL < 150:
			hist.SetBinContent(i+1,2*tree.deltaNLL)
		else:
			hist.SetBinContent(i+1,-1000)

	line4		= TF1('line4','3.84', binlo-(0.5*(binhi-binlo)/(NEntries-1)), binhi+(0.5*(binhi-binlo)/(NEntries-1)))
	line4.SetLineStyle(kDashed)
	line4.SetLineWidth(4)
	c1		= TCanvas('c1','c1',1)
	#hist.SetMarkerStyle(2)
	#hist.SetMarkerSize(0.3)
	hist.SetMarkerStyle(8)
	hist.SetMarkerSize(0.5)
	#hist.SetMarkerStyle(1)
	if par == 'cwww':
		hist.GetXaxis().SetTitle('c_{WWW} / \Lambda ^2 (1/ TeV ^2)')
	if par == 'ccw':
		hist.GetXaxis().SetTitle('c_{W} / \Lambda ^2 (1/ TeV ^2)')
	if par == 'cb':
		hist.GetXaxis().SetTitle('c_{B} / \Lambda ^2 (1/ TeV ^2)')

	for i in range((NEntries-1)/2):
		j=i+1
		tree.GetEntry(j)
		if 2*tree.deltaNLL>3.84:
			continue
		else:
			limlo = tree.GetLeaf(par).GetValue()
			break
	for i in range((NEntries-1)/2):
		j=i+ (NEntries-1)/2 +1
		tree.GetEntry(j)
		if 2*tree.deltaNLL<3.84:
			continue
		else:
			limhi = tree.GetLeaf(par).GetValue()
			break
	
	linelo = TLine(limlo,-0.1,limlo,hist.GetMaximum())
	linelo.SetLineStyle(kDashed)
	linelo.SetLineWidth(4)
	linelo.SetLineColor(kRed)
	linehi = TLine(limhi,-0.1,limhi,hist.GetMaximum())
	linehi.SetLineStyle(kDashed)
	linehi.SetLineColor(kRed)
	linehi.SetLineWidth(4)
	hist.GetYaxis().SetTitle('2*deltaNLL')
	#hist.GetYaxis().SetRangeUser(0,4.2)
	#hist.GetXaxis().SetRangeUser(10,14.5)
	hist.SetMinimum(-0.1)
	hist.Draw('p')
	line4.Draw('SAME')
	linelo.Draw('SAME')
	linehi.Draw('SAME')
	c1.Update()

	error = (float(binhi)-float(binlo))/float(NEntries-1)

	print '95%% C.L. limit on %s: [%s,%s] +- %s'%(par,round(limlo,2),round(limhi,2),round(error,2))

	#tree.GetEntry(751)
	#print tree.cwww
	#print tree.deltaNLL
	#tree.GetEntry(750)
	#print tree.cwww
	#print tree.deltaNLL
	#tree.GetEntry(752)
	#print tree.cwww
	#print tree.deltaNLL

	raw_input('<>')

plots()
