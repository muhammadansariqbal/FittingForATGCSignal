import sys

import ROOT
from ROOT import *
import math as math

inputFileName = sys.argv[3]
minmax = [float(sys.argv[1]),float(sys.argv[2])]
parm = sys.argv[4]
file = TFile.Open(inputFileName,'READ')
limit = file.Get('limit')

nEntries = limit.GetEntries()-1

histo = TH1F('LLHscan','',nEntries,minmax[0],minmax[1])
limit.Draw("%s >> LLHscan"%parm,"deltaNLL",'goff')

minBin = histo.GetMinimumBin()
minBinCenter = histo.GetBinCenter(minBin)
minVal = histo.GetBinContent(minBin)#should be 0

print 'Found profile-likelihood minimum value %.3f at %.3f'%(minVal,
                                                             minBinCenter)
#skip best fit value which is entry 0
bounds68 = []
bounds95 = []
lastBelowErr68 = False
lastAboveErr68 = True
lastBelowErr95 = False
lastAboveErr95 = True
for i in xrange(nEntries):
    limit.GetEntry(i+1)
    NLL = limit.GetLeaf("deltaNLL").GetValue()
    if NLL - minVal <= 0.5:
        if lastAboveErr68:
            bounds68.append([])
            bounds68[-1].append(limit.GetLeaf(parm).GetValue())
        lastBelowErr68 = True
        lastAboveErr68 = False
    else:
        if lastBelowErr68:
            bounds68[-1].append(limit.GetLeaf(parm).GetValue())
        lastBelowErr68 = False
        lastAboveErr68 = True
    if NLL - minVal <= 1.92:
        if lastAboveErr95:
            bounds95.append([])
            bounds95[-1].append(limit.GetLeaf(parm).GetValue())
        lastBelowErr95 = True
        lastAboveErr95 = False
    else:
        if lastBelowErr95:
            bounds95[-1].append(limit.GetLeaf(parm).GetValue())
        lastBelowErr95 = False
        lastAboveErr95 = True
file.Close()

boundstext68 = ['[%.3g,%.3g]'%(b[0],b[1]) for b in bounds68]
boundstext95 = ['[%.3g,%.3g]'%(b[0],b[1]) for b in bounds95]
print '68 CL Limit on %s:'%parm,'U'.join(boundstext68)
print '95 CL Limit on %s:'%parm,'U'.join(boundstext95)

mass_Z		= 0.0911876
mass_W		= 0.080385
GF		= 11.663787
g2		= (8*mass_W**2*GF)/math.sqrt(2)
tan_theta	= math.tan(math.acos(mass_W/mass_Z))
if parm=='ccw':
	print 'limits on delta_g1z'
	limlo	= bounds95[0][0] * (mass_Z**2)/(2)
	limhi	= bounds95[0][1] * (mass_Z**2)/(2)
	print '[' + str(limlo) + ',' + str(limhi) + ']' 
elif parm=='cwww':
	print 'limits on lambda_Z :'
	limlo	= bounds95[0][0] * (3*g2*mass_W**2)/2
	limhi	= bounds95[0][1] * (3*g2*mass_W**2)/2 
	print '[' + str(limlo) + ',' + str(limhi) + ']'
elif parm=='cb':
	print 'limits on delta_kappaZ (wrong):'
	limlo	= bounds95[0][0] * (mass_W**2*tan_theta**2)/2
	limhi	= bounds95[0][1] * (mass_W**2*tan_theta**2)/2
	print '[' + str(limlo) + ',' + str(limhi) + ']'

