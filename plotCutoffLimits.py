from ROOT import  *
from array import array
import math as math
import os
import CMS_lumi, tdrstyle

gSystem.Load("%s/lib/slc6_amd64_gcc481/libHiggsAnalysisCombinedLimit.so"%os.environ["CMSSW_BASE"])

gStyle.SetOptStat(0)
gStyle.SetOptTitle(0)

def plot():

	mWVCutoff = [2400, 2700, 3000, 3300, 3600, 3900, 4200, 4500]
	cwwwLimitsLo = [-2.92, -2.41, -2.2, -2.2, -1.91, -1.76, -1.76, -1.58]
	cwwwLimitsHi = [2.99, 2.48, 2.27, 2.27, 1.98, 1.84, 1.84, 1.59]
	ccwLimitsLo = [-1000, -3.02, -2.65, -2.74, -2.38, -2.2, -2.2, -2.00]
	ccwLimitsHi = [1000, 3.83, 3.56, 3.46, 3.19, 3.02, 2.92, 2.65]
	cbOrigLimitsLo = [-13.8, -13.4, -12.6, -13, -11, -9.8, -9.8, -8.78]
	cbOrigLimitsHi = [13.4, 13, 12.2, 12.6, 11, 9.8, 9.8, 8.74]
	cbLimitsLo = []
	cbLimitsHi = []

	for i in range(len(mWVCutoff)):
		cbLimitsLo.append(cbOrigLimitsLo[i]/2.0)
		cbLimitsHi.append(cbOrigLimitsHi[i]/2.0)

	graphcwwwLimitsLo = TGraph(len(mWVCutoff), array('d',mWVCutoff), array('d',cwwwLimitsLo))
	graphcwwwLimitsHi = TGraph(len(mWVCutoff), array('d',mWVCutoff), array('d',cwwwLimitsHi))
	graphccwLimitsLo = TGraph(len(mWVCutoff), array('d',mWVCutoff), array('d',ccwLimitsLo))
        graphccwLimitsHi = TGraph(len(mWVCutoff), array('d',mWVCutoff), array('d',ccwLimitsHi))
	graphcbLimitsLo = TGraph(len(mWVCutoff), array('d',mWVCutoff), array('d',cbLimitsLo))
        graphcbLimitsHi = TGraph(len(mWVCutoff), array('d',mWVCutoff), array('d',cbLimitsHi))

	graphcwwwLimitsLo.GetXaxis().SetTitle("m_{WV} cutoff (GeV)")
	graphcwwwLimitsLo.GetXaxis().SetTitleSize(0.06)
	graphcwwwLimitsLo.GetXaxis().SetTitleOffset(0.75)
	graphcwwwLimitsLo.GetYaxis().SetRangeUser(-10,12)
	graphcwwwLimitsLo.GetYaxis().SetTitle("95% CL limits (TeV^{-2})")
	graphcwwwLimitsLo.GetYaxis().SetTitleSize(0.06)
	graphcwwwLimitsLo.GetYaxis().SetTitleOffset(0.7)

	graphcwwwLimitsLo.SetMarkerStyle(22)
	graphcwwwLimitsLo.SetMarkerColor(kBlue)
	graphcwwwLimitsLo.SetMarkerSize(1.5)

	graphcwwwLimitsHi.SetMarkerStyle(22)
        graphcwwwLimitsHi.SetMarkerColor(kBlue)
        graphcwwwLimitsHi.SetMarkerSize(1.5)

	graphccwLimitsLo.SetMarkerStyle(20)
        graphccwLimitsLo.SetMarkerColor(kRed)
        graphccwLimitsLo.SetMarkerSize(1.5)

        graphccwLimitsHi.SetMarkerStyle(20)
        graphccwLimitsHi.SetMarkerColor(kRed)
        graphccwLimitsHi.SetMarkerSize(1.5)

	graphcbLimitsLo.SetMarkerStyle(21)
        graphcbLimitsLo.SetMarkerColor(kGreen+2)
        graphcbLimitsLo.SetMarkerSize(1.5)

        graphcbLimitsHi.SetMarkerStyle(21)
        graphcbLimitsHi.SetMarkerColor(kGreen+2)
        graphcbLimitsHi.SetMarkerSize(1.5)

	c1 = TCanvas()

	graphcwwwLimitsLo.Draw('AP')
	graphcbLimitsLo.Draw('SAME P')
        graphcbLimitsHi.Draw('SAME P')
	graphccwLimitsLo.Draw('SAME P')
        graphccwLimitsHi.Draw('SAME P')
	graphcwwwLimitsLo.Draw('SAME P')
	graphcwwwLimitsHi.Draw('SAME P')

	CMS_lumi.lumiTextOffset=0.1
	CMS_lumi.relPosY    = -0.05
	CMS_lumi.relExtraDY = 0.24
	CMS_lumi.CMS_lumi(c1,4,11)

	leg=TLegend(0.12,0.8,0.92,0.9)
	leg.SetTextSize(0.05)
	leg.SetNColumns(3)
	leg.SetFillColor(0)
        leg.SetFillStyle(0)
        leg.SetBorderSize(0)
        leg.SetLineColor(0)
        leg.SetLineWidth(0)
        leg.SetLineStyle(0)
	leg.SetTextFont(42)
	leg.AddEntry(graphcwwwLimitsLo,"c_{WWW} / #Lambda^{2}","P")
	leg.AddEntry(graphccwLimitsLo,"c_{W} / #Lambda^{2}","P")
	leg.AddEntry(graphcbLimitsLo,"0.5 #times c_{B} / #Lambda^{2}","P")
	leg.Draw("SAME")

	c1.Update()
	c1.SaveAs("cutoffLimits.pdf")

	raw_input('Done')

plot()

