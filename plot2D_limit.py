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
parser.add_option('--POI',dest='POI',help='parameters of interest')
parser.add_option('--pval',dest='pval',help='values of parameters')
(options,args) = parser.parse_args()

POI	= []
pval	= []
for i in options.POI.split(','):
	POI.append(i)
for j in options.pval.split(','):
	pval.append(j)

def plot2d():
	#raise RuntimeError("fix this")
	path		= './'
	wsname		= 'higgsCombine_%s_%s_%s_%s.MultiDimFit.mH120.root'%(POI[0],pval[0],POI[1],pval[1])
	fileInATGC	= TFile.Open(path+wsname)
	tree		= fileInATGC.Get('limit')
	
	NEntries	= tree.GetEntries()
	nbins		= int(math.sqrt(NEntries))
	binslo		= [-20,-30]
	binshi		= [20,30]

	if 'cwww' in POI:
		#tree.GetEntry(1)
		#binslo[POI.index('cwww')]	= int(tree.cwww%-100)-1
		#tree.GetEntry(NEntries-1)
		#binshi[POI.index('cwww')]	= int(tree.cwww%100)+1
		if POI.index('cwww') == 0:
			labelx = 'c_{WWW}'
		else:
			labely = 'c_{WWW}'
	if 'ccw' in wsname:
		#tree.GetEntry(1)
		#binslo[POI.index('ccw')]	= int(tree.ccw%-100)-1
		#tree.GetEntry(NEntries-1)
		#binshi[POI.index('ccw')]	= int(tree.ccw%100)+1
		if POI.index('ccw') == 0:
			labelx = 'c_{W}'
		else:
			labely = 'c_{W}'
	if 'cb' in wsname:
		#tree.GetEntry(1)
		#binslo[POI.index('cb')]		= int(tree.cb%-100)-1
		#tree.GetEntry(NEntries-1)
		#binshi[POI.index('cb')]		= int(tree.cb%100)+1
		if POI.index('cb') == 0:
			labelx = 'c_{B}'
		else:
			labely = 'c_{B}'

	hist		= TH2F('hist','hist',nbins,binslo[0],binshi[0],nbins,binslo[1],binshi[1])
	#hist		= TH2F('hist','hist',315,-18,18,315,-30,30)

	for i in range(1,NEntries):
		if i%2500==0:
			print i
		tree.GetEntry(i)
		if 2*limit.deltaNLL<599:
			hist.Fill(limit.GetLeaf(POI[0]).GetValue(),limit.GetLeaf(POI[1]).GetValue(),2*limit.deltaNLL)
		else:
			hist.Fill(limit.GetLeaf(POI[0]).GetValue(),limit.GetLeaf(POI[1]).GetValue(),-1000)

	c1 = TCanvas('c1','c1',1)
	hist.SetMinimum(-0.1)
	#hist.SetMaximum(3.84)
	hist.GetXaxis().SetTitle(labelx)
	hist.GetYaxis().SetTitle(labely)

	hist.GetZaxis().SetTitle('2*deltaNLL')
	hist.Draw('prof cont1')
	c1.Update()
	raw_input('</>')

plot2d()
