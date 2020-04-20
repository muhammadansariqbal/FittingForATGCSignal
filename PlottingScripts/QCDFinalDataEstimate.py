from ROOT import *
from array import array
from optparse import OptionParser
import math as math
import os
import CMS_lumi, tdrstyle

gStyle.SetOptStat(0)
gStyle.SetOptTitle(0)
TGaxis.SetMaxDigits(3)

def plot(ch="el"):

	fileMC = TFile.Open("/afs/cern.ch/user/m/maiqbal/private/FittingForATGC/Background/CMSSW_5_3_32/src/FittingForATGCBackground/InputTrees/%s/tree_WJets_%s.root"%(ch,ch))
	fileAntiIso = TFile.Open("/afs/cern.ch/work/m/maiqbal/private/aTGC/QCDEstimate/data_%s_antiIso.root"%ch)
	treeMC = fileMC.Get("BasicTree")
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
		fakeRate = 0.0298 # percent
	else:
		METCut = 40
		fakeRate = 0.0166 # percent

	# Draw the histograms
	treeAntiIso.Draw("jet_mass_softdrop_PUPPI >> histAntiIso(37,40,150)", "(1) * (jet_mass_softdrop_PUPPI>40. && jet_mass_softdrop_PUPPI<150. && nbtag==0)", "HIST")
	histAntiIso = gDirectory.Get("histAntiIso")

	treeMC.Draw("jet_mass_softdrop_PUPPI >> histMC(37,40,150)", "(100 * totEventWeight) * (jet_pt>200. && jet_tau21_PUPPI<0.55 && jet_mass_softdrop_PUPPI<150. && jet_mass_softdrop_PUPPI>40. && W_pt>200. && fabs(deltaR_LeptonWJet)>TMath::Pi()/2 && fabs(deltaPhi_WJetMet)>2. && fabs(deltaPhi_WJetWlep)>2. && nbtag==0 && pfMET>%s && MWW_SD>900 && MWW_SD<4500)"%(METCut), "HIST SAME")
	histMC = gDirectory.Get("histMC")

	treeAntiIso.Draw("jet_mass_softdrop_PUPPI >> histIso(37,40,150)", "(%s) * (jet_mass_softdrop_PUPPI>40. && jet_mass_softdrop_PUPPI<150. && nbtag==0)"%fakeRate, "HIST SAME")
	histIso = gDirectory.Get("histIso")

	# Set properties
	histAntiIso.SetLineColor(kAzure+10)
	histAntiIso.SetLineWidth(2)
	histAntiIso.GetXaxis().SetTitle("m_{SD} (GeV)")
	histAntiIso.GetYaxis().SetTitle("Events / 3 GeV")
	histAntiIso.GetXaxis().SetTitleSize(0.06)
	histAntiIso.GetXaxis().SetTitleOffset(0.68)
        histAntiIso.GetYaxis().SetTitleSize(0.06)
        histAntiIso.GetYaxis().SetTitleOffset(0.7)

	histMC.SetLineColor(kGreen+2)
	histMC.SetLineWidth(2)
	
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
	leg = TLegend(0.45,0.65,0.95,0.75)
	leg.SetFillColor(0)
	leg.SetFillStyle(0)
	leg.SetBorderSize(0)
	leg.SetLineColor(0)
	leg.SetLineWidth(0)
	leg.SetLineStyle(0)
	leg.SetTextFont(42)
	leg.AddEntry(histAntiIso,"QCD anti-isolated signal+sideband","L")
	leg.AddEntry(histIso,"Final QCD isolated estimate #times 100","L")
	leg.AddEntry(histMC,"W+jets #times 100","L")
	leg.Draw()
	legends.append(leg)
	
	canvas.Update()
	canvas.SaveAs('QCDFinalDataEstimate_%s.pdf'%(ch))
	#raw_input("Drawn")
	pad.Delete()
	canvas.Delete()

for channel in ["el","mu"]:
	plot(channel)
