from ROOT import  *
from array import array
import math as math
import os
import CMS_lumi, tdrstyle

gStyle.SetOptStat(0)
gStyle.SetOptTitle(0)

def plot(par='cwww'):

	mWVCutoff = [2400, 2700, 3000, 3300, 3600, 3900, 4200, 4500]

	cwwwLimitsLo = [-2.96, -2.47, -2.23, -2.20, -1.95, -1.78, -1.78, -1.58]
	cwwwLimitsHi = [2.99, 2.48, 2.25, 2.21, 1.956, 1.77, 1.79, 1.59]
	cwLimitsLo = [-3.79, -3.05, -2.69, -2.81, -2.46, -2.25, -2.27, -2.00]
	cwLimitsHi = [4.46, 3.80, 3.50, 3.43, 3.11, 2.93, 2.87, 2.65]
	cbLimitsLo = [-14.14, -13.66, -12.74, -13.04, -11.4, -9.9, -9.96, -8.78]
	cbLimitsHi = [13.04, 12.8, 12.12, 12.42, 10.94, 9.54, 9.6, 8.74]

	lzLimitsLo = [-0.0122, -0.0102, -0.0092, -0.0091, -0.0080, -0.0074, -0.0073, -0.0065]
	lzLimitsHi = [0.0123, 0.0103, 0.0093, 0.0092, 0.0081, 0.0074, 0.0074, 0.0066]
	dg1zLimitsLo = [-0.0110, -0.0097, -0.0087, -0.0091, -0.0078, -0.0069, -0.0069, -0.0061]
	dg1zLimitsHi = [0.0123, 0.0110, 0.0101, 0.0101, 0.009, 0.0082, 0.0081, 0.0074]
	dkzLimitsLo = [-0.0121, -0.0118, -0.0112, -0.0115, -0.0101, -0.0088, -0.0088, -0.0079]
	dkzLimitsHi = [0.0131, 0.0126, 0.0118, 0.0121, 0.0106, 0.0092, 0.0092, 0.0082]

	limitsLo = {'cwww':cwwwLimitsLo, 'cw':cwLimitsLo, 'cb':cbLimitsLo, 'lz':lzLimitsLo, 'dg1z':dg1zLimitsLo, 'dkz':dkzLimitsLo}
	limitsHi = {'cwww':cwwwLimitsHi, 'cw':cwLimitsHi, 'cb':cbLimitsHi, 'lz':lzLimitsHi, 'dg1z':dg1zLimitsHi, 'dkz':dkzLimitsHi}
	limLo = limitsLo[par]
	limHi = limitsHi[par]

	c = []
	xerr = []
	cerr = []

	for i in range(len(mWVCutoff)):
		c.append((limLo[i]+limHi[i])/2)
		xerr.append(0)
		cerr.append((limHi[i]-limLo[i])/2)

	# Add one point outside to make smooth curve behave properly
	mWVCutoff.append(4800)
	c.append(c[len(c)-1])
	xerr.append(xerr[len(xerr)-1])
	cerr.append(cerr[len(cerr)-1])

	graph = TGraphErrors(len(mWVCutoff), array('d',mWVCutoff), array('d',c), array('d',xerr), array('d',cerr))

	graph.GetXaxis().SetTitle("m_{WV} cutoff (GeV)")
	graph.GetXaxis().SetTitleSize(0.06)
	graph.GetXaxis().SetTitleOffset(0.75)
	graph.GetXaxis().SetRangeUser(2400,4500)
	#graph.GetYaxis().SetRangeUser(-10,12)
	label = {'cwww':'c_{WWW}/#Lambda^{2} (TeV^{-2})', 'cw':'c_{W}/#Lambda^{2} (TeV^{-2})', 'cb':'c_{B}/#Lambda^{2} (TeV^{-2})', 'lz':'#lambda_{Z}', 'dg1z':'#Deltag_{1}^{Z}', 'dkz':'#Delta#kappa_{Z}'}
	graph.GetYaxis().SetTitle(label[par])
	graph.GetYaxis().SetTitleSize(0.06)
	graph.GetYaxis().SetTitleOffset(0.7)

	graph.SetFillColor(kYellow-3)
	graph.SetFillStyle(3004)
	graph.SetLineColor(kMagenta)
	graph.SetLineWidth(2)

	c1 = TCanvas("can","can",800,640)
	pad = TPad("pad","pad",0.,0.,1.,1.)
	c1.cd()
	pad.Draw()
	pad.cd()

	graph.Draw("A4 SAME")
	medianLine = TLine(2400,0.,4500.,0.); medianLine.SetLineWidth(2); medianLine.SetLineColor(kBlack); medianLine.SetLineStyle(9); medianLine.Draw("SAME")
	graph.Draw("P SAME")

	# Lumi text
        CMS_lumi.lumiTextSize = 0.0
        CMS_lumi.writeExtraText = True
        CMS_lumi.extraText = "Private work"
        CMS_lumi.cmsTextSize = 0.6
	CMS_lumi.extraOverCmsTextSize = 0.57
        CMS_lumi.relPosX    = 0.01
        CMS_lumi.relPosY    = -0.070
        CMS_lumi.relExtraDX = 0.13
        CMS_lumi.relExtraDY = 0.40
        CMS_lumi.CMS_lumi(pad,4,11)
        CMS_lumi.cmsTextSize=0.0
        CMS_lumi.writeExtraText = False
        CMS_lumi.lumiTextSize = 0.52
        CMS_lumi.lumiTextOffset = 0.1
        CMS_lumi.CMS_lumi(pad,4,11)

	leg = TLegend(0.62,0.76,0.94,0.9)
	leg.SetFillColor(0)
	leg.SetFillStyle(0)
	leg.SetBorderSize(0)
	leg.SetLineColor(0)
	leg.SetLineWidth(0)
	leg.SetLineStyle(0)
	leg.SetTextFont(42)
	leg.AddEntry(graph,"95% CL interval","E")
	leg.AddEntry(medianLine,"SM value","L")
	leg.Draw()
	
	#c1.Update()
	c1.SaveAs("cutoffLimits_%s.pdf"%par)
	pad.Delete()
	c1.Delete()
	
for par in ['cwww','cw','cb','lz','dg1z','dkz']:
	plot(par)

