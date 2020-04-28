from ROOT import  *
from array import array
from optparse import OptionParser
import math as math
import os

def limits(cutoff='4500'):

	mWVCutoff.append(float(cutoff))

	fileIn   = TFile.Open('higgsCombine_%s_%s_cutoff%s.MultiDimFit.mH120.root'%(POI,pval,cutoff))
        treeObs      = fileIn.Get('limit')
        nEntries     = treeObs.GetEntries()

        treeObs.GetEntry(1)
        xObs    = []
        yObs    = []

        for i in range(nEntries-1):
                treeObs.GetEntry(i+1)
                if 2*treeObs.deltaNLL < 8:
                        xObs.append(treeObs.GetLeaf(POI).GetValue())
                        yObs.append(2*treeObs.deltaNLL)
        graphObs        = TGraph(len(xObs),array('d',xObs),array('d',yObs))

	#graphObs.GetHistogram().GetMinimumBin(minBin,minDNLL)
	treeObs.GetEntry(0)
	c.append(treeObs.GetLeaf(POI).GetValue())
	
	for i in range(1000):
                if (graphObs.Eval(float(pval)*i/1000.) > 3.84):
                        limHi.append(float(pval)*i/1000.)
                        break

	for i in range(1000):
		if (graphObs.Eval(-float(pval)*i/1000.) > 3.84):
			limLo.append(-float(pval)*i/1000.)
			break

parser	= OptionParser()
parser.add_option('--POI',dest='POI',help='parameter of interest')
parser.add_option('--pval',dest='pval',help='value of parameter')
(options,args) = parser.parse_args()

POI		= options.POI
pval		= options.pval

mWVCutoff = []
c = []
limLo = []
limHi = []

for cutoff in ['2400','2700','3000','3300','3600','3900','4200','4500']:
	limits(cutoff)

print "Central values:"
print c
print "Lower limits:"
print limLo
print "Upper limits:"
print limHi
