from ROOT import *
from array import array
from optparse import OptionParser
import math as math
import os
import CMS_lumi, tdrstyle

gStyle.SetOptStat(0)
gStyle.SetOptTitle(0)

def plot(ch="el",par=61,legStr="c_{WWW}/#Lambda^{2}=-3.6 TeV^{-2}"):

	fileInWW = TFile.Open("Input/WW-aTGC_%s.root"%(ch))
	fileInWZ = TFile.Open("Input/WZ-aTGC_%s.root"%(ch))
	treeInWW = fileInWW.Get("BasicTree")
	treeInWZ = fileInWZ.Get("BasicTree")

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

	# Set channel specific things
	if ch == "el":
		METCut = 110
		if par == 59 or par == 63:
			yRange = 60
		else:
			yRange = 40
	else:
		METCut = 40
		if par == 59 or par == 63:
			yRange = 80
		else:
			yRange = 60

	# Draw the four histograms
	treeInWW.Draw("jet_mass_softdrop_PUPPI >> histSMWW(27,65,105)", "(totEventWeight*aTGCWeights[61]) * (jet_pt>200. && jet_tau21_PUPPI<0.55 && jet_mass_softdrop_PUPPI<105. && jet_mass_softdrop_PUPPI>65. && W_pt>200. && fabs(deltaR_LeptonWJet)>TMath::Pi()/2 && fabs(deltaPhi_WJetMet)>2. && fabs(deltaPhi_WJetWlep)>2. && nbtag==0 && pfMET>%s && MWW>900 && MWW<4500)"%(METCut), "HIST")
	histSMWW = gDirectory.Get("histSMWW")
	histSMWW.GetYaxis().SetRangeUser(0,yRange)

	treeInWZ.Draw("jet_mass_softdrop_PUPPI >> histSMWZ(27,65,105)", "(totEventWeight*aTGCWeights[61]) * (jet_pt>200. && jet_tau21_PUPPI<0.55 && jet_mass_softdrop_PUPPI<105. && jet_mass_softdrop_PUPPI>65. && W_pt>200. && fabs(deltaR_LeptonWJet)>TMath::Pi()/2 && fabs(deltaPhi_WJetMet)>2. && fabs(deltaPhi_WJetWlep)>2. && nbtag==0 && pfMET>%s && MWW>900 && MWW<4500)"%(METCut), "SAME HIST")
	histSMWZ = gDirectory.Get("histSMWZ")

	treeInWW.Draw("jet_mass_softdrop_PUPPI >> histWW(27,65,105)", "(totEventWeight*aTGCWeights[%s]) * (jet_pt>200. && jet_tau21_PUPPI<0.55 && jet_mass_softdrop_PUPPI<105. && jet_mass_softdrop_PUPPI>65. && W_pt>200. && fabs(deltaR_LeptonWJet)>TMath::Pi()/2 && fabs(deltaPhi_WJetMet)>2. && fabs(deltaPhi_WJetWlep)>2. && nbtag==0 && pfMET>%s && MWW>900 && MWW<4500)"%(par,METCut), "SAME HIST")
	histWW = gDirectory.Get("histWW")

	treeInWZ.Draw("jet_mass_softdrop_PUPPI >> histWZ(27,65,105)", "(totEventWeight*aTGCWeights[%s]) * (jet_pt>200. && jet_tau21_PUPPI<0.55 && jet_mass_softdrop_PUPPI<105. && jet_mass_softdrop_PUPPI>65. && W_pt>200. && fabs(deltaR_LeptonWJet)>TMath::Pi()/2 && fabs(deltaPhi_WJetMet)>2. && fabs(deltaPhi_WJetWlep)>2. && nbtag==0 && pfMET>%s && MWW>900 && MWW<4500)"%(par,METCut), "SAME HIST")
	histWZ = gDirectory.Get("histWZ")

	# Set properties
	histSMWW.SetLineColor(kBlack)
	histSMWW.SetFillColor(kRed)
	histSMWW.GetXaxis().SetTitle("m_{SD} (GeV)")
	histSMWW.GetXaxis().SetTitleSize(0.06)
	histSMWW.GetXaxis().SetTitleOffset(0.68)
	histSMWW.GetYaxis().SetTitle("Events / 1.5 GeV")
        histSMWW.GetYaxis().SetTitleSize(0.06)
        histSMWW.GetYaxis().SetTitleOffset(0.7)
	
	histSMWZ.SetLineColor(kBlack)
	histSMWZ.SetFillColor(kRed+2)

	histWW.SetLineColor(kBlue)
	histWW.SetLineWidth(2)

	histWZ.SetLineColor(kBlue+2)
	histWZ.SetLineWidth(2)

	# Lumi text
	CMS_lumi.lumiTextSize = 0.0
	#CMS_lumi.writeExtraText = True
	CMS_lumi.cmsTextSize = 0.65625
	CMS_lumi.relPosY    = -0.07875
	CMS_lumi.relExtraDX = 0.175
	CMS_lumi.relExtraDY = 0.21
	CMS_lumi.CMS_lumi(pad,4,11)
	CMS_lumi.cmsTextSize=0.0
	CMS_lumi.writeExtraText = False
	CMS_lumi.lumiTextSize = 0.56875
	CMS_lumi.lumiTextOffset = 0.175
	CMS_lumi.CMS_lumi(pad,4,11)

	# Channel text
	pt = TPaveText(0.125,0.8,0.325,0.975, "blNDC")
	pt.SetFillStyle(0)
	pt.SetBorderSize(0)
	pt.SetTextAlign(13)
	pt.SetTextSize(0.05775)
	if ch=="el":
	    pt.AddText("Electron channel")
	else:
	    pt.AddText("Muon channel")
	pt.Draw("SAME")
	paveTexts.append(pt)

	# Legend
	leg = TLegend(0.55,0.65,0.89,0.9)
	leg.SetFillColor(0)
	leg.SetFillStyle(0)
	leg.SetBorderSize(0)
	leg.SetLineColor(0)
	leg.SetLineWidth(0)
	leg.SetLineStyle(0)
	leg.SetTextFont(42)
	leg.AddEntry(histSMWW,"SM WW","F")
	leg.AddEntry(histSMWZ,"SM WZ","F")
	leg.AddEntry(histWW,"WW %s"%(legStr),"L")
	leg.AddEntry(histWZ,"WZ %s"%(legStr),"L")
	leg.Draw()
	legends.append(leg)
	
	canvas.Update()
	canvas.SaveAs('signalHist_mJ_%s_%s.pdf'%(ch,par))
	#raw_input("Drawn")
	pad.Delete()
	canvas.Delete()

par = [11, 111, 51, 71, 59, 63]
legStr = ["c_{WWW}/#Lambda^{2}=3.6 TeV^{-2}",
	"c_{WWW}/#Lambda^{2}=-3.6 TeV^{-2}",
	"c_{W}/#Lambda^{2}=4.5 TeV^{-2}",
	"c_{W}/#Lambda^{2}=-4.5 TeV^{-2}",
	"c_{B}/#Lambda^{2}=20 TeV^{-2}",
	"c_{B}/#Lambda^{2}=-20 TeV^{-2}"
]

for channel in ["el","mu"]:
	for ind in range(6):
		plot(channel,par[ind],legStr[ind])
