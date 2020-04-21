from ROOT import *
from array import array
from optparse import OptionParser
import math as math
import os
import CMS_lumi, tdrstyle

gStyle.SetOptStat(0)
gStyle.SetOptTitle(0)

def plot(ch="el",var="msd"):

	fileInWJets = TFile.Open("/afs/cern.ch/user/m/maiqbal/private/FittingForATGC/Background/CMSSW_5_3_32/src/FittingForATGCBackground/InputTrees/%s/tree_WJets_%s.root"%(ch,ch))
	fileInQCD = TFile.Open("/afs/cern.ch/work/m/maiqbal/private/aTGC/FittingInputTrees/%s/tree_QCD_%s.root"%(ch,ch))
	treeInWJets = fileInWJets.Get("BasicTree")
	treeInQCD = fileInQCD.Get("BasicTree")

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

	yRange = 260
	if ch == "el":
		METCut = 110
	else:
		METCut = 40
		if var == "msd":
			yRange = 480

	# Draw the histograms
	if var == "msd":
		treeInWJets.Draw("jet_mass_softdrop_PUPPI >> histWJets(43,30,160)", "(totEventWeight) * (jet_pt>200. && jet_tau21_PUPPI<0.55 && jet_mass_softdrop_PUPPI<150. && jet_mass_softdrop_PUPPI>40. && W_pt>200. && fabs(deltaR_LeptonWJet)>TMath::Pi()/2 && fabs(deltaPhi_WJetMet)>2. && fabs(deltaPhi_WJetWlep)>2. && nbtag==0 && pfMET>%s && MWW_SD>900 && MWW_SD<4500)"%(METCut), "HIST")
		histWJets = gDirectory.Get("histWJets")
		histWJets.GetYaxis().SetRangeUser(0,yRange)
	else:
		treeInWJets.Draw("MWW_SD >> histWJets(36,900,4500)", "(totEventWeight) * (jet_pt>200. && jet_tau21_PUPPI<0.55 && jet_mass_softdrop_PUPPI<150. && jet_mass_softdrop_PUPPI>40. && W_pt>200. && fabs(deltaR_LeptonWJet)>TMath::Pi()/2 && fabs(deltaPhi_WJetMet)>2. && fabs(deltaPhi_WJetWlep)>2. && nbtag==0 && pfMET>%s && MWW_SD>900 && MWW_SD<4500)"%(METCut), "HIST")
		pad.SetLogy()
		histWJets = gDirectory.Get("histWJets")

	if var == "msd":
		treeInQCD.Draw("jet_mass_softdrop_PUPPI >> histQCD(43,30,160)", "(totEventWeight) * (jet_pt>200. && jet_tau21_PUPPI<0.55 && jet_mass_softdrop_PUPPI<150. && jet_mass_softdrop_PUPPI>40. && nbtag==0 && MWW_SD>900 && MWW_SD<4500)", "SAME HIST")
	else:
		treeInQCD.Draw("MWW_SD >> histQCD(36,900,4500)", "(totEventWeight) * (jet_pt>200. && jet_tau21_PUPPI<0.55 && jet_mass_softdrop_PUPPI<150. && jet_mass_softdrop_PUPPI>40. && nbtag==0 && MWW_SD>900 && MWW_SD<4500)", "SAME HIST")
	histQCD = gDirectory.Get("histQCD")

	# Set properties
	histWJets.SetLineColor(kGreen+2)
	histWJets.SetLineWidth(2)
	if var == "msd":
		histWJets.GetXaxis().SetTitle("m_{SD} (GeV)")
		histWJets.GetYaxis().SetTitle("Events / 3 GeV")
	else:
		histWJets.GetXaxis().SetTitle("m_{WV} (GeV)")
		histWJets.GetYaxis().SetTitle("Events / 100 GeV")
	histWJets.GetXaxis().SetTitleSize(0.06)
	histWJets.GetXaxis().SetTitleOffset(0.68)
        histWJets.GetYaxis().SetTitleSize(0.06)
        histWJets.GetYaxis().SetTitleOffset(0.7)
	
	histQCD.SetLineColor(kGray+1)
	histQCD.SetLineWidth(2)

	# Lumi text
        CMS_lumi.lumiTextSize = 0.0
        CMS_lumi.writeExtraText = True
        CMS_lumi.extraText = "Simulation, private work"
        CMS_lumi.cmsTextSize = 0.6
        CMS_lumi.relPosX    = 0.01
        CMS_lumi.relPosY    = -0.070
        CMS_lumi.relExtraDX = 0.13
        CMS_lumi.relExtraDY = 0.22
        CMS_lumi.CMS_lumi(pad,4,11)
        CMS_lumi.cmsTextSize=0.0
        CMS_lumi.writeExtraText = False
        CMS_lumi.lumiTextSize = 0.52
        CMS_lumi.lumiTextOffset = 0.1
        CMS_lumi.CMS_lumi(pad,4,11)

	# Channel text
	pt = TPaveText(0.125,0.8,0.325,0.975, "blNDC")
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
	leg = TLegend(0.70,0.73,0.94,0.85)
	leg.SetFillColor(0)
	leg.SetFillStyle(0)
	leg.SetBorderSize(0)
	leg.SetLineColor(0)
	leg.SetLineWidth(0)
	leg.SetLineStyle(0)
	leg.SetTextFont(42)
	leg.AddEntry(histWJets,"W+jets","L")
	leg.AddEntry(histQCD,"QCD","L")
	leg.Draw()
	legends.append(leg)
	
	canvas.Update()
	canvas.SaveAs('QCDMC_%s_%s.pdf'%(ch,var))
	#raw_input("Drawn")
	pad.Delete()
	canvas.Delete()

for channel in ["el","mu"]:
	for var in ["msd","mwv"]:
		plot(channel,var)
