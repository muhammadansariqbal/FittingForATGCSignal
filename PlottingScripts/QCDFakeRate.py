from ROOT import *
from array import array
from optparse import OptionParser
import math as math
import os
import CMS_lumi, tdrstyle

gStyle.SetOptStat(0)
gStyle.SetOptTitle(0)

def plot(ch="el"):

	fileMC = TFile.Open("/afs/cern.ch/work/m/maiqbal/private/aTGC/QCDEstimate/MC_%s.root"%ch)
	fileIso = TFile.Open("/afs/cern.ch/work/m/maiqbal/private/aTGC/QCDEstimate/tree_data_%s.root"%ch)
	fileAntiIso = TFile.Open("/afs/cern.ch/work/m/maiqbal/private/aTGC/QCDEstimate/data_%s_antiIso.root"%ch)
	treeMC = fileMC.Get("BasicTree")
	treeIso = fileIso.Get("BasicTree")
	treeAntiIso = fileAntiIso.Get("treeDumper/BasicTree")

	# These things just have to be kept in memory so that ROOT does not make them disappear in the next loop
	pads = []
	paveTexts = []
	legends = []

	canvas = TCanvas(ch,ch,800,640)
	pad = TPad("pad","pad",0.,0.,1.,1.)
	pads.append(pad)

	canvas.cd()
	pad.Draw()
	pad.cd()

	if ch == "el":
		METCut = 110
	else:
		METCut = 40

	# Draw the histograms
	treeIso.Draw("jet_mass_softdrop_PUPPI >> histIso(30,0,45)", "(100) * (jet_mass_softdrop_PUPPI>0 && jet_mass_softdrop_PUPPI<160 && (jet_mass_softdrop_PUPPI<40. || jet_mass_softdrop_PUPPI>150.) && MWW_SD>900)", "HIST")
        histIso = gDirectory.Get("histIso")

	treeMC.Draw("jet_mass_softdrop_PUPPI >> histMC(30,0,45)", "(100*totEventWeight) * (jet_mass_softdrop_PUPPI>0 && jet_mass_softdrop_PUPPI<160 && (jet_mass_softdrop_PUPPI<40. || jet_mass_softdrop_PUPPI>150.) && MWW_SD>900)", "HIST SAME")
        histMC = gDirectory.Get("histMC")

	treeAntiIso.Draw("jet_mass_softdrop_PUPPI >> histAntiIso(30,0,45)", "(1) * (jet_mass_softdrop_PUPPI>0 && jet_mass_softdrop_PUPPI<160 && (jet_mass_softdrop_PUPPI<40. || jet_mass_softdrop_PUPPI>150.) && nbtag==0)", "HIST")
	histAntiIso = gDirectory.Get("histAntiIso")

	# Subtract MC from isolated
	histIso.Add(histMC,-1)
	histIso.Draw("HIST SAME")

	# Set properties
	histAntiIso.SetLineColor(kAzure+7)
	histAntiIso.SetLineWidth(2)
	histAntiIso.GetXaxis().SetTitle("m_{SD} (GeV)")
	histAntiIso.GetYaxis().SetTitle("Events / 1.5 GeV")
	histAntiIso.GetXaxis().SetTitleSize(0.06)
	histAntiIso.GetXaxis().SetTitleOffset(0.68)
        histAntiIso.GetYaxis().SetTitleSize(0.06)
        histAntiIso.GetYaxis().SetTitleOffset(0.7)
	
	histIso.SetLineColor(kGray+1)
	histIso.SetLineWidth(2)

	# Lumi text
	CMS_lumi.lumiTextSize = 0.0
	CMS_lumi.writeExtraText = False
	CMS_lumi.extraText = "Simulation"
	CMS_lumi.cmsTextSize = 0.6
	CMS_lumi.relPosX    = 0.128
	CMS_lumi.relPosY    = -0.070
	CMS_lumi.relExtraDX = 0.16
	CMS_lumi.relExtraDY = 0.25
	CMS_lumi.CMS_lumi(pad,4,11)
	CMS_lumi.cmsTextSize=0.0
	CMS_lumi.writeExtraText = False
	CMS_lumi.lumiTextSize = 0.52
	CMS_lumi.lumiTextOffset = 0.16
	CMS_lumi.CMS_lumi(pad,4,11)

	# Channel text
	pt = TPaveText(0.19,0.8,0.39,0.975, "blNDC")
	pt.SetFillStyle(0)
	pt.SetBorderSize(0)
	pt.SetTextAlign(13)
	pt.SetTextSize(0.0575)
	if ch=="el":
	    pt.AddText("Electron channel")
	else:
	    pt.AddText("Muon channel")
	pt.Draw("SAME")
	paveTexts.append(pt)

	# Legend
	leg = TLegend(0.45,0.55,0.95,0.65)
	leg.SetFillColor(0)
	leg.SetFillStyle(0)
	leg.SetBorderSize(0)
	leg.SetLineColor(0)
	leg.SetLineWidth(0)
	leg.SetLineStyle(0)
	leg.SetTextFont(42)
	leg.AddEntry(histAntiIso,"QCD anti-isolated side-sideband","L")
	leg.AddEntry(histIso,"QCD isolated side-sideband #times 100","L")
	leg.Draw()
	legends.append(leg)
	
	canvas.Update()
	canvas.SaveAs('QCDFakeRate_%s.pdf'%(ch))
	#raw_input("Drawn")
	pad.Delete()
	canvas.Delete()

for channel in ["el","mu"]:
	plot(channel)
