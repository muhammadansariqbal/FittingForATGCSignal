from ROOT import  *
from array import array
from optparse import OptionParser
from ConfigParser import SafeConfigParser
import math as math
import random
import os
import CMS_lumi, tdrstyle

gSystem.Load('%s/lib/slc6_amd64_gcc481/libHiggsAnalysisCombinedLimit.so'%os.environ['CMSSW_BASE'])
from ROOT import RooPoly3Pdf, RooChiSqPdf, RooErfExpPdf, RooErfPowExpPdf, RooErfPowPdf, RooErfPow2Pdf, RooExpNPdf, RooAlpha4ExpNPdf, RooExpTailPdf, RooAlpha4ExpTailPdf, Roo2ExpPdf


parser        = OptionParser()
parser.add_option('-n', '--newtrees', action='store_true', dest='newtrees', default=False, help='recreate aTGC histograms')
parser.add_option('-p', '--plots', action='store_true', dest='Make_plots', default=False, help='make plots')
parser.add_option('--savep', action='store_true', dest='savep', default=False, help='save plots')
parser.add_option('-b', action='store_true', dest='batch', default=False, help='batch mode')
parser.add_option('-c', '--ch', dest='chan', default='elmu', help='channel, el, mu or elmu')
parser.add_option('--noatgcint', action='store_true', dest='noatgcint', default=False, help='set atgc-interference coefficients to zero')
parser.add_option('--printatgc', action='store_true', default=False, help='print atgc-interference contribution')
parser.add_option('--atgc', action='store_true', dest='atgc', default=False, help='use anomalous coupling parametrization instead of EFT')
parser.add_option('--binWidth', dest='binWidth', default='100.', help='Use different MWV binnings, required for Asimov datasets')
parser.add_option('--cutoff', dest='mWVCutoff', default='4500.', help='mWV range upper limit in GeV')

(options,args) = parser.parse_args()


gStyle.SetOptStat(0)
gStyle.SetOptTitle(0)
if options.batch:
        gROOT.SetBatch(True)

if not os.path.isdir('Output'):
        os.system('mkdir Output')


class Prepare_workspace_4limit:
    
        def __init__(self,ch,mlvj_lo,mlvj_hi,mlvj_binWid):
            
            self.POI                    = ['cwww','ccw','cb']
            self.PAR_TITLES             = {'cwww' : '#frac{c_{WWW}}{#Lambda^{2}}', 'ccw' : '#frac{c_{W}}{#Lambda^{2}}', 'cb' : '#frac{c_{B}}{#Lambda^{2}}'}#latex titles 
            self.PAR_MAX                = {'cwww' : 3.6, 'ccw' : 4.5, 'cb' : 20}#atgc points
            self.ch                     = ch

            self.binlo                  = mlvj_lo                #lower bound on invariant mass
            self.binhi                  = mlvj_hi                #upper bound

            self.channel                = "WV_"+self.ch
            self.nbins                  = int((self.binhi-self.binlo)/mlvj_binWid)
            
            self.WS                     = RooWorkspace("w")        #final workspace
            self.wtmp                   = RooWorkspace('wtmp')
            
            self.fitresults             = []
            ##nuisance parameter to change all slope parameters by certain percentage (bigger for cb in WZ-cateogry)
            self.eps                    = RooRealVar('slope_nuis','slope_nuis',1,0,2)
            self.eps.setConstant(kTRUE)
            self.eps4cbWZ               = RooFormulaVar('rel_slope_nuis4cbWZ','rel_slope_nuis4cbWZ','1+3*(@0-1)',RooArgList(self.eps))
            
            ##read workspace containing background pdfs
            fileInWs                    = TFile.Open('Input/wwlvj_%s_HPV_workspace.root'%(self.ch))
            w                           = fileInWs.Get('workspace4limit_')
            self.rrv_mass_lvj           = w.var('rrv_mass_lvj')
            self.rrv_mass_lvj.SetTitle('M_{WV}')
            self.rrv_mass_lvj.setRange(self.binlo,self.binhi)
            self.rrv_mass_j             = w.var('rrv_mass_j')

            #names of bkg contributions and regions
            self.bkgs        = ['WJets','TTbar','WW','WZ','STop']
            self.regions     = ['sig','sb_lo','sb_hi']


        #read trees containing aTGC WW and WZ events and fill them into histograms
        def Read_ATGCtree(self,ch='el'):

            #cut on MET has to be applied
            if self.ch=='el':
                METCUT        = 80.
            elif self.ch=='mu':
                METCUT        = 40.
            else:
                raise RuntimeError('no such channel %s'%self.ch)
            
            print '######### Making histograms for aTGC working points #########'
            hists4scale        = {}
            
            for WV in ['WW','WZ']:
                #create 3 histograms for each aTGC parameter (positive, negative and positive-negative working point)
                for para in self.POI:
                    hists4scale['c_pos_%s_hist_%s'%(WV,para)] = TH1F('c_pos_%s_hist_%s'%(WV,para),'c_pos_%s_hist_%s'%(WV,para),self.nbins,self.binlo,self.binhi);
                    hists4scale['c_neg_%s_hist_%s'%(WV,para)] = TH1F('c_neg_%s_hist_%s'%(WV,para),'c_neg_%s_hist_%s'%(WV,para),self.nbins,self.binlo,self.binhi);
                    hists4scale['c_dif_%s_hist_%s'%(WV,para)] = TH1F('c_dif_%s_hist_%s'%(WV,para),'c_dif_%s_hist_%s'%(WV,para),self.nbins,self.binlo,self.binhi);
                    hists4scale['c_pos_%s_hist_%s'%(WV,para)].Sumw2(kTRUE)
                    hists4scale['c_neg_%s_hist_%s'%(WV,para)].Sumw2(kTRUE)
                    hists4scale['c_dif_%s_hist_%s'%(WV,para)].Sumw2(kTRUE)
                #add histograms for SM and all aTGC parameters unequal to zero
                hists4scale['c_SM_%s_hist'%WV]                  = TH1F('c_SM_%s_hist'%WV,'c_SM_%s_hist'%WV,self.nbins,self.binlo,self.binhi);                
                hists4scale['c_%s_histall3'%WV]                 = TH1F('c_%s_histall3'%WV,'c_%s_histall3'%WV,self.nbins,self.binlo,self.binhi);
                hists4scale['c_SM_%s_hist'%WV].Sumw2(kTRUE)
                hists4scale['c_%s_histall3'%WV].Sumw2(kTRUE)

		# Add histograms for two aTGC parameters positive
		hists4scale['c_cwww_ccw_%s_hist'%WV]=TH1F('c_cwww_ccw_%s_hist'%WV,'c_cwww_ccw_%s_hist'%WV,self.nbins,self.binlo,self.binhi);
		hists4scale['c_ccw_cb_%s_hist'%WV]=TH1F('c_ccw_cb_%s_hist'%WV,'c_ccw_cb_%s_hist'%WV,self.nbins,self.binlo,self.binhi);
		hists4scale['c_cwww_ccw_%s_hist'%WV].Sumw2(kTRUE)
		hists4scale['c_ccw_cb_%s_hist'%WV].Sumw2(kTRUE)

		# Add histograms for aTGC-aTGC interference terms
		hists4scale['c_int_cwww_ccw_%s_hist'%WV]=TH1F('c_int_cwww_ccw_%s_hist'%WV,'c_int_cwww_ccw_%s_hist'%WV,self.nbins,self.binlo,self.binhi);
		hists4scale['c_int_ccw_cb_%s_hist'%WV]=TH1F('c_int_ccw_cb_%s_hist'%WV,'c_int_ccw_cb_%s_hist'%WV,self.nbins,self.binlo,self.binhi);
		hists4scale['c_int_cwww_ccw_%s_hist'%WV].Sumw2(kTRUE)
		hists4scale['c_int_ccw_cb_%s_hist'%WV].Sumw2(kTRUE)

                print 'reading %s-aTGC_%s.root'%(WV,self.ch)
                fileInATGC        = TFile.Open('Input/%s-aTGC_%s.root'%(WV,self.ch))
                treeInATGC        = fileInATGC.Get('BasicTree')

                lumi_tmp         = 35922.

                for i in range(treeInATGC.GetEntries()):
                    if i%10000==0:
                            print str(i) + '/' + str(treeInATGC.GetEntries())
                    treeInATGC.GetEntry(i)
                    MWW                = treeInATGC.MWW_SD
                    #apply cuts
                    #using whole mj-range (sideband and signal region)
                    if treeInATGC.jet_pt>200. and treeInATGC.jet_tau21_PUPPI<0.55 and treeInATGC.W_pt>200. and treeInATGC.deltaR_LeptonWJet>math.pi/2. and treeInATGC.jet_mass_softdrop_PUPPI>40 and treeInATGC.jet_mass_softdrop_PUPPI<150 and abs(treeInATGC.deltaPhi_WJetMet)>2. and abs(treeInATGC.deltaPhi_WJetWlep)>2. and treeInATGC.nbtag==0 and treeInATGC.pfMET>METCUT and MWW>self.binlo:
			weight_part = treeInATGC.totEventWeight
			aTGC        = treeInATGC.aTGCWeights                #contains weights for different workingpoints
			#all3hists4scale['c_%s_histall3'%WV].Fill(MWW,aTGC[123] * weight_part)
			hists4scale['c_%s_histall3'%WV].Fill(MWW,aTGC[123] * weight_part)
			#SM
			hists4scale['c_SM_%s_hist'%WV].Fill(MWW,aTGC[61] * weight_part)
			#cwww
			hists4scale['c_pos_%s_hist_cwww'%WV].Fill(MWW,aTGC[11] * weight_part)
			hists4scale['c_neg_%s_hist_cwww'%WV].Fill(MWW,aTGC[111] * weight_part)
			#ccw
			hists4scale['c_pos_%s_hist_ccw'%WV].Fill(MWW,aTGC[51] * weight_part)
			hists4scale['c_neg_%s_hist_ccw'%WV].Fill(MWW,aTGC[71] * weight_part)
			#cb
			hists4scale['c_pos_%s_hist_cb'%WV].Fill(MWW,aTGC[59] * weight_part)
			hists4scale['c_neg_%s_hist_cb'%WV].Fill(MWW,aTGC[63] * weight_part)
			#ccw-SM interference
			hists4scale['c_dif_%s_hist_ccw'%WV].Fill(MWW,(aTGC[51]-aTGC[71]) * weight_part)
			#cb-SM interference
			hists4scale['c_dif_%s_hist_cb'%WV].Fill(MWW,(aTGC[59]-aTGC[63]) * weight_part)
			#cwww+ccw
			hists4scale['c_cwww_ccw_%s_hist'%WV].Fill(MWW,aTGC[1] * weight_part)
			#ccw+cb
			hists4scale['c_ccw_cb_%s_hist'%WV].Fill(MWW,aTGC[49] * weight_part)
			#cwww-ccw interference
			hists4scale['c_int_cwww_ccw_%s_hist'%WV].Fill(MWW, ((aTGC[1]-aTGC[101])-(aTGC[11]-aTGC[111])) * weight_part)
			#ccw-cb interference
			hists4scale['c_int_ccw_cb_%s_hist'%WV].Fill(MWW, ((aTGC[49]-aTGC[69])-(aTGC[51]-aTGC[71])) * weight_part)

		# Fit exponential to the aTGC-aTGC interference histograms (This is to avoid doing this via gen-level files)
		hists4scale['c_int_cwww_ccw_%s_hist'%WV].Fit("expo")
		a5_val=hists4scale['c_int_cwww_ccw_%s_hist'%WV].GetFunction("expo").GetParameter(1)
		hists4scale['c_int_ccw_cb_%s_hist'%WV].Fit("expo")
		a7_val=hists4scale['c_int_ccw_cb_%s_hist'%WV].GetFunction("expo").GetParameter(1)
		
		# Write the slopes to workspace
		a5=RooRealVar('a5_%s'%WV,'a5_%s'%WV,-0.0001,-0.01,0.01)
		a7=RooRealVar('a7_%s'%WV,'a7_%s'%WV,-0.001,-0.01,0.01)
		a5.setVal(a5_val)
		a7.setVal(a7_val)
		a5.setConstant(true)
		a7.setConstant(true)
		self.Import_to_ws(self.wtmp, [a5,a7])
		
		# Write normalizations to workspace
		# Not used now (histograms are retrieved later, RooDataHists are created and sumEntries is used on them); although can be done via this too, the results are identical
		N3645=RooRealVar('N_cwww_ccw_36_45_%s'%WV,'N_cwww_ccw_36_45_%s'%WV,hists4scale['c_cwww_ccw_%s_hist'%WV].Integral())
		N4520=RooRealVar('N_ccw_cb_45_20_%s'%WV,'N_ccw_cb_45_20_%s'%WV,hists4scale['c_ccw_cb_%s_hist'%WV].Integral())
		N36=RooRealVar('N_cwww_36_%s'%WV,'N_cwww_36_%s'%WV,hists4scale['c_pos_%s_hist_cwww'%WV].Integral())
		N36_=RooRealVar('N_cwww__36_%s'%WV,'N_cwww__36_%s'%WV,hists4scale['c_neg_%s_hist_cwww'%WV].Integral())
		N45=RooRealVar('N_ccw_45_%s'%WV,'N_ccw_45_%s'%WV,hists4scale['c_pos_%s_hist_ccw'%WV].Integral())
                N45_=RooRealVar('N_ccw__45_%s'%WV,'N_ccw__45_%s'%WV,hists4scale['c_neg_%s_hist_ccw'%WV].Integral())
		N20=RooRealVar('N_cb_20_%s'%WV,'N_cb_20_%s'%WV,hists4scale['c_pos_%s_hist_cb'%WV].Integral())
                N20_=RooRealVar('N_cb__20_%s'%WV,'N_cb__20_%s'%WV,hists4scale['c_neg_%s_hist_cb'%WV].Integral())
		self.Import_to_ws(self.wtmp, [N3645,N4520,N36,N36_,N45,N45_,N20,N20_])

            #write histograms to file
            fileOut        = TFile.Open('Output/hists4scale_%s_WV_aTGC-%s_%s.root'%(self.ch,self.binlo,self.binhi),'recreate')
            for key in hists4scale:
                hists4scale[key].Write()
            print '--------> Written to file ' + fileOut.GetName()
            fileOut.Close()

        def Make_plots(self,rrv_x,cat,fitres):
                
            can     = []
            can2    = []
            plots   = []
            plots2  = []
            pads    = []
            
            channel = self.ch+'_'+cat

            for i in range(3):
                rrv_x.setRange(self.binlo,self.binhi)
                p       = rrv_x.frame(self.binlo,self.binhi)
                p2      = rrv_x.frame(self.binlo,self.binhi)
                c       = TCanvas(self.POI[i]+'-',self.POI[i]+'-',1)
                c.cd()
                pad1        = TPad('pad1_%s'%self.POI[i],'pad1_%s'%self.POI[i],0.,0.25,1.,1.)
                pad2        = TPad('pad2_%s'%self.POI[i],'pad2_%s'%self.POI[i],0.,0.02,1.,0.25)
                c2          = TCanvas(self.POI[i]+'+',self.POI[i]+'+',1)
                c2.cd()
                pad3        = TPad('pad3_%s'%self.POI[i],'pad3_%s'%self.POI[i],0.,0.25,1.,1.)
                pad4        = TPad('pad4_%s'%self.POI[i],'pad4_%s'%self.POI[i],0.,0.02,1.,0.25)
                p2pads      = [pad1,pad2,pad3,pad4]
                can.append(c)
                can2.append(c2)
                plots.append(p)
                plots2.append(p2)
                pads.append(p2pads)

            for i in range(3):

                can[i].cd()
                pads[i][0].Draw()
                pads[i][1].Draw()
                pads[i][0].SetLeftMargin(0.1)
                pads[i][1].SetLeftMargin(0.1)
                
                norm = self.wtmp.function('normfactor_3d_%s'%channel)

                for j in range(3):
                        self.wtmp.var(self.POI[j]).setVal(0)
                self.wtmp.data('SMdatahist_%s'%cat).plotOn(plots[i],RooFit.MarkerColor(kBlack),RooFit.LineColor(kBlack),RooFit.DataError(RooAbsData.SumW2),RooFit.DrawOption('E0'),RooFit.Name('SMdata'))
                normvalSM        = norm.getVal() * self.wtmp.data('SMdatahist_%s'%cat).sumEntries()

		#print "SM data for fitting in Matlab"
	        #histPrint=plots[i].getHist();
	        #for i in range(histPrint.GetN()):
		#	print histPrint.GetY()[i]
		#print "The corresponding sumw2 errors"
		#for i in range(histPrint.GetN()):
		#	print histPrint.GetErrorY(i)
		#raw_input('SM data printed.')

                self.wtmp.pdf('aTGC_model_%s'%channel).plotOn(plots[i],RooFit.LineColor(kBlack),RooFit.Normalization(normvalSM, RooAbsReal.NumEvent),RooFit.Name('SMmodel'))

                self.wtmp.data('neg_datahist_%s_%s'%(cat,self.POI[i])).plotOn(plots[i],RooFit.MarkerColor(kBlue),RooFit.LineColor(kBlue),RooFit.DataError(RooAbsData.SumW2),RooFit.DrawOption('E0'),RooFit.Name('atgcdata'))
                self.wtmp.var(self.POI[i]).setVal(-self.PAR_MAX[self.POI[i]])
                normvalneg = norm.getVal() * self.wtmp.data('SMdatahist_%s'%cat).sumEntries()
                self.wtmp.pdf('aTGC_model_%s'%channel).plotOn(plots[i],RooFit.LineColor(kBlue),RooFit.Normalization(normvalneg, RooAbsReal.NumEvent),RooFit.Name('atgcmodel'))
        
                pullhist = plots[i].pullHist('atgcdata','atgcmodel')
                
                plotmax        = 1000
                if self.ch == 'el':
                    plotmin = 1e-4
                    if cat == 'WZ':
                        plotmin = 3e-4
                elif self.ch == 'mu':
                    plotmin = 1e-4
                    if cat == 'WZ':
                        plotmin = 1e-4
                        plotmax = 500
                if cat == 'WV':
                    plotmin = 1e-2
                    plotmax = 1.5e2
                plots[i].GetYaxis().SetRangeUser(plotmin,plotmax)
                pads[i][0].cd()
                pads[i][0].SetLogy()
                plots[i].SetTitle('')
                plots[i].GetYaxis().SetTitle('arbitrary units')
                #plots[i].GetXaxis().SetRangeUser(900,options.mWVCutoff)

                plots[i].Draw()
                ndof        = (self.binhi-self.binlo)/100 - 4
                plots[i].Print()
                
                parlatex        = ['#frac{c_{WWW}}{#Lambda^{2}}','#frac{c_{W}}{#Lambda^{2}}','#frac{c_{B}}{#Lambda^{2}}']
                leg        = TLegend(0.11,0.2,0.4,0.6)
                leg.SetFillStyle(0)
                leg.SetBorderSize(0)
                leg.AddEntry(plots[i].findObject('SMdata'),'MC '+parlatex[i]+'=0 TeV^{-2}','le')
                leg.AddEntry(plots[i].findObject('SMmodel'),'signal model '+parlatex[i]+'=0 TeV^{-2}','l')
                leg.AddEntry(plots[i].findObject('atgcdata'),'MC '+parlatex[i]+'='+str(-self.PAR_MAX[self.POI[i]])+' TeV^{-2}','le')
                leg.AddEntry(plots[i].findObject('atgcmodel'),'signal model '+parlatex[i]+'='+str(-self.PAR_MAX[self.POI[i]])+' TeV^{-2}','l')
                leg.Draw()
                leg.Print()
                
                pads[i][1].cd()
                ratio_style = TH1D('ratio_style','ratio_style',int(self.binhi-self.binlo)/100,self.binlo,self.binhi)
                ratio_style.SetMarkerStyle(21)
                ratio_style.SetMaximum(3)
                ratio_style.SetMinimum(-3)
                ratio_style.GetYaxis().SetNdivisions(7)
                ratio_style.GetYaxis().SetTitle('#frac{MC-Fit}{Error}')
                ratio_style.GetYaxis().SetLabelSize(0.125)
                ratio_style.GetYaxis().SetTitleSize(0.2)
                ratio_style.GetYaxis().SetTitleOffset(0.2)
                ratio_style.Draw("")
                pullhist.SetLineColor(kBlue)
                pullhist.Draw("SAME E1")

                CMS_lumi.CMS_lumi(can[i], 4, 11)
                pt = TPaveText(0.6,0.82,0.90,0.90, "blNDC")
                pt.SetFillStyle(0)
                pt.SetBorderSize(0)
                pt.SetTextAlign(32)
                pt.SetTextSize(0.035)
                if (self.ch=='el'):
                    pt.AddText("Electron channel")
                elif (self.ch=='mu'):
                    pt.AddText("Muon channel")
                if (cat=='WW'):
                    pt.AddText("WW")
                elif (cat=='WZ'):
                    pt.AddText("WZ")
                pt.Draw("SAME")

                can[i].Update()
                if options.savep:
                    if not os.path.isdir('Plots'):
                        os.system('mkdir Plots')
                    can[i].SaveAs('Plots/%s_neg_%s.pdf'%(self.POI[i],channel))
                    can[i].SaveAs('Plots/%s_neg_%s.png'%(self.POI[i],channel))
                

                for j in range(3):
                        self.wtmp.var(self.POI[j]).setVal(0)
                self.wtmp.data('SMdatahist_%s'%cat).plotOn(plots2[i],RooFit.MarkerColor(kBlack),RooFit.LineColor(kBlack),RooFit.DataError(RooAbsData.SumW2),RooFit.DrawOption('E0'))
                self.wtmp.data('pos_datahist_%s_%s'%(cat,self.POI[i])).plotOn(plots2[i],RooFit.MarkerColor(kBlue),RooFit.LineColor(kBlue),RooFit.DataError(RooAbsData.SumW2),RooFit.DrawOption('E0'))

                self.wtmp.pdf('aTGC_model_%s'%channel).plotOn(plots2[i],RooFit.LineColor(kBlack),RooFit.Normalization(normvalSM, RooAbsReal.NumEvent))
                self.wtmp.var(self.POI[i]).setVal(self.PAR_MAX[self.POI[i]])
                normvalpos = norm.getVal() * self.wtmp.data('SMdatahist_%s'%cat).sumEntries()

                self.wtmp.pdf('aTGC_model_%s'%channel).plotOn(plots2[i],RooFit.LineColor(kBlue),RooFit.Normalization(normvalpos, RooAbsReal.NumEvent))
                
                self.wtmp.data('pos_datahist_%s_%s'%(cat,self.POI[i])).plotOn(plots2[i],RooFit.MarkerColor(kBlue),RooFit.LineColor(kBlue),RooFit.DataError(RooAbsData.SumW2),RooFit.DrawOption('E'))
                plots2[i].GetYaxis().SetRangeUser(plotmin,plotmax)
                plots2[i].GetYaxis().SetTitle('arbitrary units')
                can2[i].cd()
                pads[i][2].Draw()
                pads[i][3].Draw()
                pads[i][2].SetLeftMargin(0.1)
                pads[i][3].SetLeftMargin(0.1)
                plots2[i].SetTitle('')
                pads[i][2].SetLogy()
                pads[i][2].cd()
                plots2[i].Draw()
		leg2        = TLegend(0.11,0.2,0.4,0.6)
                leg2.SetFillStyle(0)
                leg2.SetBorderSize(0)
                leg2.AddEntry(plots[i].findObject('SMdata'),'MC '+parlatex[i]+'=0 TeV^{-2}','le')
                leg2.AddEntry(plots[i].findObject('SMmodel'),'signal model '+parlatex[i]+'=0 TeV^{-2}','l')
                leg2.AddEntry(plots[i].findObject('atgcdata'),'MC '+parlatex[i]+'=+'+str(self.PAR_MAX[self.POI[i]])+' TeV^{-2}','le')
                leg2.AddEntry(plots[i].findObject('atgcmodel'),'signal model '+parlatex[i]+'=+'+str(self.PAR_MAX[self.POI[i]])+' TeV^{-2}','l')
		leg2.Draw()
		leg2.Print()
                pullhist2 = plots2[i].pullHist('h_pos_datahist_%s_%s'%(cat,self.POI[i]),'aTGC_model_%s_Norm[rrv_mass_lvj]'%channel)
                pads[i][3].cd()
                ratio_style.Draw("")
                pullhist2.SetLineColor(kBlue)
                pullhist2.Draw("E1")

                CMS_lumi.CMS_lumi(can2[i], 4, 11)
                pt.Draw("SAME")

                can2[i].Update()
                if options.savep:
                    can2[i].SaveAs('Plots/%s_pos_%s.pdf'%(self.POI[i],channel))
                    can2[i].SaveAs('Plots/%s_pos_%s.png'%(self.POI[i],channel))
                    
            if not options.batch:
                raw_input('plots plotted')

            
        #function to import multiple items from a list into a workspace
        def Import_to_ws(self,workspace,items,recycle=0):
            for item in items:
                if recycle:
                    getattr(workspace,'import')(item,RooFit.RecycleConflictNodes())
                else:
                    getattr(workspace,'import')(item)


        def Make_signal_pdf(self,rrv_x,sample):
            
            channel        = self.ch+'_'+sample                #needed for variables that differ for WW and WZ
            
            cwww    = RooRealVar('cwww','cwww',0,-36,36);
            ccw     = RooRealVar('ccw','ccw',0,-45,45);
            cb      = RooRealVar('cb','cb',0,-200,200);
            cwww.setConstant(kTRUE);
            ccw.setConstant(kTRUE);
            cb.setConstant(kTRUE);
   
            #get SM and other histograms and make RooDataHists
            fileInHist      = TFile.Open('Output/hists4scale_%s_WV_aTGC-%s_%s.root'%(self.ch,self.binlo,self.binhi))
            rrv_x.setRange(self.binlo,self.binhi)
            SMdatahist      = RooDataHist('SMdatahist_%s'%sample,'SMdatahist_%s'%sample,RooArgList(rrv_x),fileInHist.Get('c_SM_%s_hist'%sample))
	    cwwwPosDataHist = RooDataHist('cwwwPosDataHist_%s'%sample,'cwwwPosDataHist_%s'%sample,RooArgList(rrv_x),fileInHist.Get('c_pos_%s_hist_cwww'%sample))
	    cwwwNegDataHist = RooDataHist('cwwwNegDataHist_%s'%sample,'cwwwNegDataHist_%s'%sample,RooArgList(rrv_x),fileInHist.Get('c_neg_%s_hist_cwww'%sample))
	    ccwPosDataHist  = RooDataHist('ccwPosDataHist_%s'%sample,'ccwPosDataHist_%s'%sample,RooArgList(rrv_x),fileInHist.Get('c_pos_%s_hist_ccw'%sample))
	    ccwNegDataHist  = RooDataHist('ccwNegDataHist_%s'%sample,'ccwNegDataHist_%s'%sample,RooArgList(rrv_x),fileInHist.Get('c_neg_%s_hist_ccw'%sample))
	    cbPosDataHist   = RooDataHist('cbPosDataHist_%s'%sample,'cbPosDataHist_%s'%sample,RooArgList(rrv_x),fileInHist.Get('c_pos_%s_hist_cb'%sample))
	    cbNegDataHist   = RooDataHist('cbNegDataHist_%s'%sample,'cbNegDataHist_%s'%sample,RooArgList(rrv_x),fileInHist.Get('c_neg_%s_hist_cb'%sample))
	    cwwwccwDataHist = RooDataHist('cwwwccwDataHist_%s'%sample,'cwwwccwDataHist_%s'%sample,RooArgList(rrv_x),fileInHist.Get('c_cwww_ccw_%s_hist'%sample))
	    ccwcbDataHist   = RooDataHist('ccwcbDataHist_%s'%sample,'ccwcbDataHist_%s'%sample,RooArgList(rrv_x),fileInHist.Get('c_ccw_cb_%s_hist'%sample))
            fileInHist.Close()

            #make SM pdf, simple exponential
            #SMdatahist.printMultiline(SMdatahist.defaultPrintStream(), 1, kTRUE);
            #raw_input('SM data printed.');
            a1_4fit         = RooRealVar('a_SM_4fit_%s'%channel,'a_SM_4fit_%s'%channel,-0.1,-2.,0.)
            a1              = RooFormulaVar('a_SM_%s'%channel,'a_SM_%s'%channel,'@0*@1',RooArgList(a1_4fit,self.eps))
            SMPdfMain       = RooExponential('SMPdfMain_%s'%channel,'SMPdfMain_%s'%channel,rrv_x,a1)
            aCorr           = RooRealVar('a_SM_Corr_%s'%channel,'a_SM_Corr_%s'%channel,-0.1,-2.,0.)
            SMPdfCorr       = RooExponential('SMPdfCorr_%s'%channel,'SMPdfCorr_%s'%channel,rrv_x,aCorr)
            SMCombCoeff     = RooRealVar('SMCombCoeff_%s'%channel,'SMCombCoeff_%s',0.5)
            SMPdf           = RooAddPdf('SMPdf_%s'%channel,'SMPdf_%s'%channel,SMPdfMain,SMPdfCorr,SMCombCoeff)
            #actual fit to determine SM shape parameter a1_4fit
            fitresSM        = SMPdf.fitTo(SMdatahist, RooFit.SumW2Error(kTRUE), RooFit.Save(kTRUE))
            self.fitresults.append(fitresSM)
            a1_4fit.setConstant(kTRUE)
            aCorr.setConstant(kTRUE)
            #coefficient for SM term and other terms in final signal function
            N_SM            = RooRealVar('N_SM_%s'%channel,'N_SM_%s'%channel,SMdatahist.sumEntries())
	    N_3645          = RooRealVar('N_3645_%s'%channel,'N_3645_%s'%channel,cwwwccwDataHist.sumEntries())
	    N_4520          = RooRealVar('N_4520_%s'%channel,'N_4520_%s'%channel,ccwcbDataHist.sumEntries())
            N_36            = RooRealVar('N_36_%s'%channel,'N_36_%s'%channel,cwwwPosDataHist.sumEntries())
            N__36           = RooRealVar('N__36_%s'%channel,'N__36_%s'%channel,cwwwNegDataHist.sumEntries())
            N_45            = RooRealVar('N_45_%s'%channel,'N_45%s'%channel,ccwPosDataHist.sumEntries())
            N__45           = RooRealVar('N__45%s'%channel,'N__45%s'%channel,ccwNegDataHist.sumEntries())
            N_20            = RooRealVar('N_20%s'%channel,'N_20%s'%channel,cbPosDataHist.sumEntries())
            N__20           = RooRealVar('N__20%s'%channel,'N__20%s'%channel,cbNegDataHist.sumEntries())

            self.Import_to_ws(self.wtmp,[cwww,ccw,cb,self.eps4cbWZ,SMdatahist,SMdatahist,N_SM])
            
            #define parameter ranges for error function, only needed for WZ
            if sample=='WZ':
                if self.ch=='el':
                    Erf_width_cwww      = RooRealVar('Erf_width_cwww_%s'%channel,'Erf_width_cwww_%s'%channel,1000.,500.,1500.)
                    Erf_width_ccw       = RooRealVar('Erf_width_ccw_%s'%channel,'Erf_width_ccw_%s'%channel,1500.,1000.,2000.)
                elif self.ch=='mu':
                    Erf_width_cwww      = RooRealVar('Erf_width_cwww_%s'%channel,'Erf_width_cwww_%s'%channel,1000.,500.,7500.)
                    Erf_width_ccw       = RooRealVar('Erf_width_ccw_%s'%channel,'Erf_width_ccw_%s'%channel,1500.,500.,2000.)
                Erf_offset_cwww     = RooRealVar('Erf_offset_cwww_%s'%channel,'Erf_offset_cwww_%s'%channel,1000.,500.,1500.)
                Erf_offset_ccw      = RooRealVar('Erf_offset_ccw_%s'%channel,'Erf_offset_ccw_%s'%channel,1500.,500.,2500.)
                Erf_offset_cwww.setConstant(kTRUE)
                Erf_width_cwww.setConstant(kTRUE)
                Erf_offset_ccw.setConstant(kTRUE)
                Erf_width_ccw.setConstant(kTRUE)
                self.Import_to_ws(self.wtmp,[Erf_width_cwww,Erf_offset_cwww,Erf_width_ccw,Erf_offset_ccw])
                
            for i in range(len(self.POI)):
                s_name          = self.POI[i] + '_' + channel #added to parameter names
                fileInHist      = TFile.Open('Output/hists4scale_%s_WV_aTGC-%s_%s.root'%(self.ch,self.binlo,self.binhi))
                rrv_x.setRange(self.binlo,self.binhi)                
                pos_datahist    = RooDataHist('pos_datahist_%s_%s'%(sample,self.POI[i]),'pos_datahist_%s_%s'%(sample,self.POI[i]),RooArgList(rrv_x),fileInHist.Get('c_pos_%s_hist_%s'%(sample,self.POI[i])))
                neg_datahist    = RooDataHist('neg_datahist_%s_%s'%(sample,self.POI[i]),'neg_datahist_%s_%s'%(sample,self.POI[i]),RooArgList(rrv_x),fileInHist.Get('c_neg_%s_hist_%s'%(sample,self.POI[i])))
                dif_datahist    = RooDataHist('dif_datahist_%s_%s'%(sample,self.POI[i]),'dif_datahist_%s_%s'%(sample,self.POI[i]),RooArgList(rrv_x),fileInHist.Get('c_dif_%s_hist_%s'%(sample,self.POI[i])))

                SMWW        = RooDataHist('SMWW_4scale','SMWW_4scale',RooArgList(rrv_x),fileInHist.Get('c_SM_WW_hist'))
                posWW       = RooDataHist('posWW_4scale_%s'%self.POI[i],'posWW_4scale_%s'%self.POI[i],RooArgList(rrv_x),fileInHist.Get('c_pos_WW_hist_%s'%self.POI[i]))
                negWW       = RooDataHist('negWW_4scale_%s'%self.POI[i],'negWW_4scale_%s'%self.POI[i],RooArgList(rrv_x),fileInHist.Get('c_neg_WW_hist_%s'%self.POI[i]))
                SMWZ        = RooDataHist('SMWZ_4scale','SMWZ_4scale',RooArgList(rrv_x),fileInHist.Get('c_SM_WZ_hist'))
                posWZ       = RooDataHist('posWZ_4scale_%s'%self.POI[i],'posWZ_4scale_%s'%self.POI[i],RooArgList(rrv_x),fileInHist.Get('c_pos_WZ_hist_%s'%self.POI[i]))
                negWZ       = RooDataHist('negWZ_4scale_%s'%self.POI[i],'negWZ_4scale_%s'%self.POI[i],RooArgList(rrv_x),fileInHist.Get('c_neg_WZ_hist_%s'%self.POI[i]))
                fileInHist.Close()
                
                #import datasets to wtmp and final workspace WS
                self.Import_to_ws(self.wtmp,[pos_datahist,neg_datahist,dif_datahist])
                self.Import_to_ws(self.WS,[pos_datahist,neg_datahist,dif_datahist])
                
                #get scaling parabel from yields
                #FIXME scaling to the sum of WW and WZ leads to over-estimating WW and under-estimating WZ
                #FIXME scaling to WW and WZ separately leads to a really high scaling factor for WZ
                hist4scale = TH1F('hist4scale_%s'%self.POI[i],'hist4scale_%s'%self.POI[i],3,-1.5*self.PAR_MAX[self.POI[i]],1.5*self.PAR_MAX[self.POI[i]])
                hist4scale.SetBinContent(2,1)
		if sample=='WW':
		    hist4scale.SetBinContent(1,(negWW.sumEntries())/(SMWW.sumEntries()))
                    hist4scale.SetBinContent(3,(posWW.sumEntries())/(SMWW.sumEntries()))
		else:
		    hist4scale.SetBinContent(1,(negWZ.sumEntries())/(SMWZ.sumEntries()))
                    hist4scale.SetBinContent(3,(posWZ.sumEntries())/(SMWZ.sumEntries()))
                #fit parabel
                hist4scale.Fit('pol2','0')
                fitfunc     = hist4scale.GetFunction('pol2')
                par1        = RooRealVar('par1_%s'%s_name,'par1_%s'%s_name,fitfunc.GetParameter(1));
                par1.setConstant(kTRUE);
                par2        = RooRealVar('par2_%s'%s_name,'par2_%s'%s_name,fitfunc.GetParameter(2));
                par2.setConstant(kTRUE);

                N_pos_tmp   = pos_datahist.sumEntries()
                N_neg_tmp   = neg_datahist.sumEntries()
                N_quad      = RooRealVar('N_quad_%s'%s_name,'N_quad_%s'%s_name, ((N_pos_tmp+N_neg_tmp)/2)-N_SM.getVal() )
		#N_quad      = RooRealVar('N_quad_%s'%s_name,'N_quad_%s'%s_name, 0 )
                
                #scaleshape is the relative change to SM
                scaleshape  = RooFormulaVar('scaleshape_%s'%s_name,'scaleshape_%s'%s_name, '(@0*@2+@1*@2**2)', RooArgList(par1,par2,self.wtmp.var(self.POI[i])))
                #FIXME only very few atgc events for cb in WZ sample, fit doesn't work yet -> different parametrization, starting values+ranges or leave out completely
                if sample=='WZ' and self.POI[i]=='cb':
                    N_lin       = RooRealVar('N_lin_%s'%s_name,'N_lin_%s'%s_name, 0)
                    a2_4fit     = RooRealVar('a_quad_4fit_%s'%s_name,'a_quad_4fit_%s'%s_name,-0.1,-2,0.)
                    a2          = RooFormulaVar('a_quad_nuis_%s'%s_name,'a_quad_nuis_%s'%s_name,'@0*@1',RooArgList(a2_4fit,self.eps4cbWZ))
                    a3_4fit     = RooRealVar('a_lin_4fit_%s'%s_name,'a_lin_4fit_%s'%s_name,-0.0001,-0.1,0.)
                    a3          = RooFormulaVar('a_lin_nuis_%s'%s_name,'a_lin_nuis_%s'%s_name,'@0*@1',RooArgList(a3_4fit,self.eps4cbWZ))
                    cPdf_quad   = RooExponential('Pdf_quad_%s'%s_name,'Pdf_quad_%s'%s_name,rrv_x,a2)
                else:
		    #N_lin       = RooRealVar('N_lin_%s'%s_name,'N_lin_%s'%s_name, 0 )
                    N_lin       = RooRealVar('N_lin_%s'%s_name,'N_lin_%s'%s_name, (N_pos_tmp-N_neg_tmp)/2 )
                    a2_4fit     = RooRealVar('a_quad_4fit_%s'%s_name,'a_quad_4fit_%s'%s_name,-0.001,-0.01,0.)
                    a2          = RooFormulaVar('a_quad_nuis_%s'%s_name,'a_quad_nuis_%s'%s_name,'@0*@1',RooArgList(a2_4fit,self.eps))
                    a3_4fit     = RooRealVar('a_lin_4fit_%s'%s_name,'a_lin_4fit_%s'%s_name,-0.001,-0.01,0.)
                    a3          = RooFormulaVar('a_lin_nuis_%s'%s_name,'a_lin_nuis_%s'%s_name,'@0*@1',RooArgList(a3_4fit,self.eps))
                    #simple exponential sufficient for WW events
                    if   sample == 'WW':
                        cPdf_quad       = RooExponential('Pdf_quad_%s'%s_name,'Pdf_quad_%s'%s_name,rrv_x,a2)
                    #additional error function to describe turn-on for WZ events
                    elif sample == 'WZ':
                        cPdf_quad       = RooErfExpPdf('Pdf_quad_%s'%s_name,'Pdf_quad_%s'%s_name,rrv_x,a2,self.wtmp.var('Erf_offset_%s'%s_name),self.wtmp.var('Erf_width_%s'%s_name))
                a2_4fit.setConstant(kTRUE)
                a3_4fit.setConstant(kTRUE)
                #PDF for SM interference
                cPdf_lin        = RooExponential('Pdf_lin_%s'%s_name,'Pdf_lin_%s'%s_name,rrv_x,a3)

                self.Import_to_ws(self.wtmp,[cPdf_quad,cPdf_lin],1)
                self.Import_to_ws(self.wtmp,[N_quad,N_lin,scaleshape])
                
            ###make model
            #list of all coefficients
            paralist    = RooArgList(N_SM)

            # Include aTGC-interference
            # Get parameter values of aTGC-interference from tmp workspace where they are saved in the start 
            a5_tmp      = RooRealVar('a_cwww_ccw_%s'%channel,'a_cwww_ccw_%s'%channel, self.wtmp.var('a5_%s'%sample).getVal())
            a7_tmp      = RooRealVar('a_ccw_cb_%s'%channel,'a_ccw_cb_%s'%channel, self.wtmp.var('a7_%s'%sample).getVal())
            a5_tmp.setConstant(kTRUE)
            a7_tmp.setConstant(kTRUE)
            #apply uncertainty parameter, bigger uncertainty for c_B in WZ
            a5          = RooFormulaVar('a_cwww_ccw_nuis_%s'%channel,'a_cwww_ccw_nuis_%s'%channel,'@0*@1',RooArgList(a5_tmp,self.eps))
            if sample=='WZ':
                a7           = RooFormulaVar('a_ccw_cb_nuis_%s'%channel,'a_ccw_cb_nuis_%s'%channel,'@0*@1',RooArgList(a7_tmp,self.eps4cbWZ))
            else:
                a7           = RooFormulaVar('a_ccw_cb_nuis_%s'%channel,'a_ccw_cb_nuis_%s'%channel,'@0*@1',RooArgList(a7_tmp,self.eps))
            
            Pdf_cwww_ccw    = RooExponential('Pdf_cwww_ccw_%s'%channel,'Pdf_cwww_ccw_%s'%channel,rrv_x,a5)
            Pdf_ccw_cb      = RooExponential('Pdf_ccw_cb_%s'%channel,'Pdf_ccw_cb_%s'%channel,rrv_x,a7)

            if options.noatgcint:
                cf = 0
            else:
                # cf was used originally for scaling MC to gen level interference terms, not needed anymore for that, hence set to 1
                cf = 1

            # Get other coefficients
            NSM         = N_SM.getVal()
            N3645       = N_3645.getVal()
            N4520       = N_4520.getVal()
            N36         = N_36.getVal()
            N36_        = N__36.getVal()
            N45         = N_45.getVal()
            N45_        = N__45.getVal()
            N20         = N_20.getVal()
            N20_        = N__20.getVal()

            ##define final coefficients, scaled by cf
            N_cwww_ccw      = RooRealVar('N_cwww_ccw_%s'%channel,'N_cwww_ccw_%s'%channel,\
                                            cf*((N3645+NSM)-(N36+N45)))
            N_ccw_cb        = RooRealVar('N_ccw_cb_%s'%channel,'N_ccw_cb_%s'%channel,\
                                            cf*((N4520+NSM)-(N45+N20)))

            paralist.add(RooArgList(self.wtmp.function('N_quad_%s_%s'%(self.POI[0],channel)),self.wtmp.var('cwww'),\
                                    self.wtmp.function('N_quad_%s_%s'%(self.POI[1],channel)),self.wtmp.function('N_lin_%s_%s'%(self.POI[1],channel)),self.wtmp.var('ccw'),\
                                    self.wtmp.function('N_quad_%s_%s'%(self.POI[2],channel)),self.wtmp.function('N_lin_%s_%s'%(self.POI[2],channel)),self.wtmp.var('cb')))
            paralist.add(RooArgList(N_cwww_ccw,N_ccw_cb))
            
            #parts of final signal model formula
            cwww_s      = '+@1*(@2/3.6)**2'
            ccw_s       = '+@3*(@5/4.5)**2+@4*(@5/4.5)'
            cb_s        = '+@6*(@8/20)**2+@7*(@8/20)'
            cwww_ccw_s  = '+@9*(@2/3.6)*(@5/4.5)'
            ccw_cb_s    = '+@10*(@5/4.5)*(@8/20)'
            Pdf_norm    = RooFormulaVar( 'Pdf_norm_%s'%channel, 'Pdf_norm_%s'%channel, '@0'+cwww_s+ccw_s+cb_s+cwww_ccw_s+ccw_cb_s, paralist)
            paralistN   = RooArgList()
            for i in range(11):
                paralistN.add(RooArgList(paralist.at(i)))
            paralistN.add(RooArgList(Pdf_norm))

            N1                = RooFormulaVar( 'N1_%s'%channel, 'N1_%s'%channel, '@0/@11', paralistN )
            N2                = RooFormulaVar( 'N2_%s'%channel, 'N2_%s'%channel, '(@1*(@2/3.6)**2)/@11', paralistN )
            #N3 ->no SM-interference for c_WWW
            N4                = RooFormulaVar( 'N4_%s'%channel, 'N4_%s'%channel, '(@3*(@5/4.5)**2)/@11', paralistN )
            N5                = RooFormulaVar( 'N5_%s'%channel, 'N5_%s'%channel, '(@4*(@5/4.5))/@11', paralistN )
            N6                = RooFormulaVar( 'N6_%s'%channel, 'N6_%s'%channel, '(@6*(@8/20)**2)/@11', paralistN )
            N7                = RooFormulaVar( 'N7_%s'%channel, 'N7_%s'%channel, '(@7*(@8/20))/@11', paralistN )
            N8                = RooFormulaVar( 'N8_%s'%channel, 'N8_%s'%channel, '(@9*(@2/3.6)*(@5/4.5))/@11', paralistN )
            #N9 ->no aTGC-interference for c_WWW/c_B #FIXME should be added for WZ
            N10                = RooFormulaVar( 'N10_%s'%channel,'N10_%s'%channel,'(@10*(@5/4.5)*(@8/20))/@11', paralistN )

            N_list        = RooArgList(N1,N2,N4,N5,N6,N7)
            N_list.add(RooArgList(N8,N10))
            Pdf_list        = RooArgList(SMPdf)
            Pdf_list.add(RooArgList(self.wtmp.pdf('Pdf_quad_cwww_%s'%channel),\
                                    self.wtmp.pdf('Pdf_quad_ccw_%s'%channel),self.wtmp.pdf('Pdf_lin_ccw_%s'%channel),\
                                    self.wtmp.pdf('Pdf_quad_cb_%s'%channel),self.wtmp.pdf('Pdf_lin_cb_%s'%channel)))
            Pdf_list.add(RooArgList(Pdf_cwww_ccw,Pdf_ccw_cb))
            model                = RooAddPdf('aTGC_model_%s'%channel,'aTGC_model_%s'%channel, Pdf_list, N_list)

            scale_list        = RooArgList(self.wtmp.function('scaleshape_cwww_%s'%channel), self.wtmp.function('scaleshape_ccw_%s'%channel), self.wtmp.function('scaleshape_cb_%s'%channel))
            normfactor_3d        = RooFormulaVar('normfactor_3d_%s'%channel,'normfactor_3d_%s'%channel,'1+@0+@1+@2',scale_list)
            self.wtmp.Print()

            #fit 3 pdfs for 3 atgc parameters
            for i in range(3):
                s_name        = self.POI[i] + '_' + channel
                for j in range(3):
                    self.wtmp.var(self.POI[j]).setVal(0)
                self.wtmp.var(self.POI[i]).setVal(self.PAR_MAX[self.POI[i]])

                #fit SM-interference first
                ##no SM-interference for cwww; not enough aTGC events for cb in WZ sample
                if not self.POI[i] == 'cwww' and not (sample=='WZ' and self.POI[i]=='cb'):
                    #set SM and quadratical terms to zero so only the linear term is fitted
                    N_SM_tmp = N_SM.getVal()
                    N_quad_tmp = self.wtmp.var('N_quad_%s'%s_name).getVal()
                    N_SM.setVal(0)
                    self.wtmp.var('N_quad_%s'%s_name).setVal(0)
                    
                    self.wtmp.var('a_lin_4fit_%s'%s_name).setConstant(kFALSE)
                    fitres1                = model.fitTo(self.wtmp.data('dif_datahist_%s_%s'%(sample,self.POI[i])),RooFit.Save(kTRUE), RooFit.SumW2Error(kTRUE), RooFit.Minimizer('Minuit2'))
                    self.wtmp.var('a_lin_4fit_%s'%s_name).setConstant(kTRUE)
                    self.fitresults.append(fitres1)
                    
                    N_SM.setVal(N_SM_tmp)
                    self.wtmp.var('N_quad_%s'%s_name).setVal(N_quad_tmp)

                #fit quadratical term
                self.wtmp.var('a_quad_4fit_%s'%s_name).setConstant(kFALSE)
                if sample=='WZ' and self.POI[i]!='cb':
                    self.wtmp.var('Erf_offset_%s'%s_name).setConstant(kFALSE)
                    self.wtmp.var('Erf_width_%s'%s_name).setConstant(kFALSE)
                fitres2         = model.fitTo(self.wtmp.data('pos_datahist_%s_%s'%(sample,self.POI[i])), RooFit.Save(kTRUE), RooFit.SumW2Error(kTRUE))
                fitres2         = model.fitTo(self.wtmp.data('pos_datahist_%s_%s'%(sample,self.POI[i])), RooFit.Save(kTRUE), RooFit.SumW2Error(kTRUE), RooFit.Minimizer('Minuit2'))
                self.wtmp.var('a_quad_4fit_%s'%s_name).setConstant(kTRUE)
                if sample=='WZ' and self.POI[i]!='cb':
                    self.wtmp.var('Erf_offset_%s'%s_name).setConstant(kTRUE)
                    self.wtmp.var('Erf_width_%s'%s_name).setConstant(kTRUE)
                self.fitresults.append(fitres2)
                
            for i in range(3):
                self.wtmp.var(self.POI[i]).setVal(0)

            if options.atgc:
                #go from EFT parametrization to Lagrangian approach, taken from SI-HEP-2011-17
                self.Transform2Lagrangian(N_list,Pdf_list,scale_list,channel)
            else:
                model.Print()
                self.Import_to_ws(self.wtmp,[normfactor_3d,model])
                self.Import_to_ws(self.WS,[normfactor_3d,model])
            
            #print coefficients to see contribution for all atgc-parameter positive
            for i in range(3):
                self.wtmp.var(self.POI[i]).setVal(self.PAR_MAX[self.POI[i]])
            for i in range(11):
                print paralist.at(i).GetName() + ' : ' + str(paralist.at(i).getVal())

            #print self.fitresults
            for i in range(8):
                print N_list.at(i).GetName() + ' : ' + str(N_list.at(i).getVal())


        def Write_datacard(self,w,region):
            ### make the card for this channel and plane ID
            codename    = 'WWWZ_' + region + '_' + self.ch
            bkgs_4card  = ['WJets','TTbar','STop']
            Nbkg_int    = len(bkgs_4card)
            uncert_map  = {}
            ##define uncertainties
            #                                              |------------el----------------------||------------mu----------------------|
            #                                                WW   WZ     WJets   TTbar   STop      WW     WZ     WJets   TTbar   STop
            uncert_map['lumi_13TeV']                    = [1.025,1.025  ,'-'    ,1.025  ,1.025  ,1.025  ,1.025  ,'-'    ,1.025  ,1.025]
            uncert_map['CMS_eff_vtag_tau21_sf_13TeV']   = [1.140,1.140  ,'-'    ,1.140  ,1.140  ,1.140  ,1.140  ,'-'    ,1.140  ,1.140]
            uncert_map['pdf_qqbar']                     = [1.019,1.024  ,'-'    ,1.028  ,1.002  ,1.018  ,1.025  ,'-'    ,1.027  ,1.003]
            uncert_map['QCD_scale_VV']                  = [1.058,1.048  ,'-'    ,'-'    ,'-'    ,1.060  ,1.043  ,'-'    ,'-'    ,'-'  ]
            uncert_map['QCD_scale_TTbar']               = ['-'  ,'-'    ,'-'    ,1.180  ,1.009  ,'-'    ,'-'    ,'-'    ,1.178  ,1.011]
            uncert_map['CMS_scale_j']                   = [1.043,1.024  ,'-'    ,1.044  ,1.049  ,1.038  ,1.025  ,'-'    ,1.035  ,1.030]
            uncert_map['CMS_res_j']                     = [1.018,1.027  ,'-'    ,1.018  ,1.034  ,1.006  ,1.029  ,'-'    ,1.008  ,1.009]
            uncert_map['CMS_scale_met']                 = [1.010,1.006  ,'-'    ,1.009  ,1.006  ,1.002  ,1.002  ,'-'    ,1.006  ,1.010]
            uncert_map['CMS_eff_b']                     = ['-'  ,1.001  ,'-'    ,1.010  ,1.008  ,'-'    ,1.001  ,'-'    ,1.010  ,1.008]
            uncert_map['CMS_eff_e']                     = [1.023,1.023  ,'-'    ,1.021  ,1.022  ,'-'    ,'-'    ,'-'    ,'-'    ,'-'  ]
            uncert_map['CMS_eff_m']                     = ['-'  ,'-'    ,'-'    ,'-'    ,'-'    ,1.026  ,1.024  ,'-'    ,1.018  ,1.020]
            uncert_map['CMS_scale_e']                   = [1.015,1.009  ,'-'    ,1.008  ,1.014  ,'-'     ,'-'   ,'-'    ,'-'    ,'-'  ]
            uncert_map['CMS_scale_m']                   = ['-'  ,'-'    ,'-'    ,'-'    ,'-'    ,1.017  ,1.012  ,'-'    ,1.007  ,1.011]
            uncert_map['CMS_res_e']                     = [1.001,1.002  ,'-'    ,1.003  ,1.012  ,'-'     ,'-'   ,'-'    ,'-'    ,'-'  ]
            uncert_map['CMS_res_m']                     = ['-'  ,'-'    ,'-'    ,'-'    ,'-'    ,1.001  ,1.003  ,'-'    ,'-'    ,1.003]
            bkgs                                        = ['WW','WZ','TTbar','STop']
            NlnN    = len(uncert_map)

            card = """\nimax 1  number of channels\njmax {Nbkg_int}  number of backgrounds\nkmax *  number of nuisance parameters (sources of systematical uncertainties)\n-------------""".format(Nbkg_int=Nbkg_int+1)
            for i in range(0,Nbkg_int):
                card += """\nshapes {bkg_name}\t\t\t\t {codename} {codename}_ws.root\t proc_{codename}:$PROCESS""".format(codename=codename,bkg_name=bkgs_4card[i])
            card += """\nshapes data_obs\t\t\t\t {codename} {codename}_ws.root\t proc_{codename}:$PROCESS""".format(codename=codename)    
            card += """\nshapes aTGC_WW_{region}_{channel}\t\t {codename} {codename}_ws.root\t proc_{codename}:$PROCESS""".format(codename=codename,region=region,channel=self.ch)
            card += """\nshapes aTGC_WZ_{region}_{channel}\t\t {codename} {codename}_ws.root\t proc_{codename}:$PROCESS""".format(codename=codename,region=region,channel=self.ch)
            card += """\n------------\nbin\t\t\t{codename}\nobservation {obs}\n------------\nbin\t\t\t{codename}\t\t""".format(codename=codename,obs=w.data('dataset_2d_%s_%s'%(region,self.ch)).sumEntries())
            for i in range(0,Nbkg_int+1):
                card += """\t\t{codename}""".format(codename=codename)
            card += """\nprocess\t\taTGC_WW_{region}_{channel}    """.format(region=region,channel=self.ch)
            card += """\t\taTGC_WZ_{region}_{channel}""".format(region=region,channel=self.ch)
            for i in range(0,Nbkg_int):
                card += """\t\t{bkg_name}""".format(bkg_name=bkgs_4card[i])
            card += """\nprocess\t\t-1\t\t\t\t\t\t\t0"""
            for i in range(0,Nbkg_int):
                card += """ \t\t\t\t{i}""".format(i=i+1)
            card += """\nrate\t\t1\t\t1"""
            for i in range(0,Nbkg_int):
                card += """ \t\t\t1"""
            card += """\n------------\n"""
            for uncert in uncert_map:
                card += """{uncert}\t\t\tlnN\t\t""".format(uncert=uncert)
                for i in range(5):
                    if self.ch=='mu':
                        i += 5
                    card += """{value}\t\t\t""".format(value=uncert_map[uncert][i])
                card += """\n"""

# Copy these in data card below as needed
# rrv_c_ChiSq_WJets0_{ch}  flatParam
# rrv_shift_ChiSq_WJets0_{ch}  flatParam
# rrv_b0_Poly3_WJets0_{ch}  flatParam
# rrv_b1_Poly3_WJets0_{ch}  flatParam
# rrv_b2_Poly3_WJets0_{ch}  flatParam
# rrv_b3_Poly3_WJets0_{ch}  flatParam
# Deco_WJets0_sim_{ch}_HPV_mlvj_13TeV_eig0 param 0.0 1.4
# Deco_WJets0_sim_{ch}_HPV_mlvj_13TeV_eig1 param 0.0 1.4
# Deco_WJets0_sim_{ch}_HPV_mlvj_13TeV_eig2 param 0.0 1.4
# Deco_WJets0_sim_{ch}_HPV_mlvj_13TeV_eig3 param 0.0 1.4
# Deco_TTbar_sb_{ch}_HPV_mlvj_13TeV_eig0 param 0.0 2.0
# Deco_TTbar_sb_{ch}_HPV_mlvj_13TeV_eig1 param 0.0 2.0
# Deco_TTbar_sig_{ch}_HPV_mlvj_13TeV_eig0 param 0.0 2.0
# Deco_TTbar_sig_{ch}_HPV_mlvj_13TeV_eig1 param 0.0 2.0

# Write in the data card            
            card += '''
normvar_WJets_{ch}  flatParam
rrv_c_ChiSq_WJets0_{ch}  flatParam
rrv_shift_ChiSq_WJets0_{ch}  flatParam
rrv_n_ExpN_WJets0_sb_{ch}  flatParam
rrv_c_ExpN_WJets0_sb_{ch}  flatParam
Deco_WJets0_sim_{ch}_HPV_mlvj_13TeV_eig0 param 0.0 1.815
Deco_WJets0_sim_{ch}_HPV_mlvj_13TeV_eig1 param 0.0 1.815
Deco_WJets0_sim_{ch}_HPV_mlvj_13TeV_eig2 param 0.0 1.815
Deco_WJets0_sim_{ch}_HPV_mlvj_13TeV_eig3 param 0.0 1.815
slope_nuis    param  1.0 0.05'''.format(ch=self.ch)

            print card
            cardfile = open('aC_%s.txt'%(codename),'w')
            cardfile.write(card)
            cardfile.close()


        def Transform2Lagrangian(self,N_list,Pdf_list,scale_list,channel):
            #go from EFT parametrization to Lagrangian approach, taken from SI-HEP-2011-17
            Z_mass          = 0.0911876
            W_mass          = 0.080385
            G_F             = 11.663787
            g_weak          = math.sqrt((8*G_F*W_mass**2)/(math.sqrt(2)))
            theta_W         = math.acos(W_mass/Z_mass)
            tan_theta_W     = math.tan(theta_W)
            sin_theta_W     = math.sin(theta_W)

            coeff_cb1       = RooRealVar('coeff_cb1','coeff_cb1',2/(tan_theta_W*tan_theta_W*Z_mass*Z_mass))
            coeff_cb2       = RooRealVar('coeff_cb2','coeff_cb2',2/(sin_theta_W*sin_theta_W*Z_mass*Z_mass))
            coeff_ccw       = RooRealVar('coeff_ccw','coeff_ccw',2/(Z_mass*Z_mass))
            coeff_cwww      = RooRealVar('coeff_cwww','coeff_cwww',2/(3*g_weak*g_weak*W_mass*W_mass))

            lZ              = RooRealVar('lZ','lZ',0,-10*self.PAR_MAX['cwww']*(3*g_weak*g_weak*W_mass*W_mass/2), 10*self.PAR_MAX['cwww']*(3*g_weak*g_weak*W_mass*W_mass/2))
            dg1z            = RooRealVar('dg1z','dg1z',0,-10*self.PAR_MAX['ccw']*(Z_mass*Z_mass/2), 10*self.PAR_MAX['ccw']*(Z_mass*Z_mass/2))
            dkz             = RooRealVar('dkz','dkz',0,(-10*self.PAR_MAX['ccw']-(10*self.PAR_MAX['cb']*tan_theta_W*tan_theta_W))*(W_mass*W_mass/2),(10*self.PAR_MAX['ccw']-(-10*self.PAR_MAX['cb']*tan_theta_W*tan_theta_W))*(W_mass*W_mass/2))
            dg1z.setConstant(kTRUE)
            lZ.setConstant(kTRUE)
            dkz.setConstant(kTRUE)

            cwww            = RooFormulaVar('cwww_atgc','cwww_atgc','@0*@1',RooArgList(lZ,coeff_cwww))
            ccw             = RooFormulaVar('ccw_atgc','ccw_atgc','@0*@1',RooArgList(dg1z,coeff_ccw))
            cb              = RooFormulaVar('cb_atgc','cb_atgc','@0*@1-@2*@3',RooArgList(dg1z,coeff_cb1,dkz,coeff_cb2))
            atgc_pars       = RooArgList(cwww,ccw,cb)
            N_list_atgc     = RooArgList()
            scale_list_atgc = RooArgList()

            for i in range(8):
                    customize_N        = RooCustomizer(N_list.at(i),'customize_N')
                    for j in range(3):
                            customize_N.replaceArg(self.wtmp.var(self.POI[j]),atgc_pars.at(j))
                    N_list_atgc.add(customize_N.build())
                    N_list_atgc.at(i).SetName(N_list.at(i).GetName())
            model_atgc        = RooAddPdf('aTGC_model_%s'%channel,'aTGC_model_%s'%channel, Pdf_list, N_list_atgc)
            for i in range(3):
                    customize_scale        = RooCustomizer(scale_list.at(i),'customize_scale')
                    for j in range(3):
                            customize_scale.replaceArg(self.wtmp.var(self.POI[j]),atgc_pars.at(j))
                    scale_list_atgc.add(customize_scale.build())
                    
            normfactor_3d        = RooFormulaVar('normfactor_3d_%s'%channel,'normfactor_3d_%s'%channel,'1+@0+@1+@2',scale_list_atgc)
    
            getattr(self.wtmp,'import')(model_atgc,RooFit.RecycleConflictNodes())
            getattr(self.WS,'import')(model_atgc,RooFit.RecycleConflictNodes())
            getattr(self.wtmp,'import')(normfactor_3d,RooFit.RecycleConflictNodes())
            getattr(self.WS,'import')(normfactor_3d,RooFit.RecycleConflictNodes())
            self.WS.Print()
            raw_input(channel)


        ########################
        ######MAIN CODE#########
        ########################
        def Make_input(self):

            #prepare variables, parameters and temporary workspace
            if options.newtrees:
                self.Read_ATGCtree(self.ch)
            
            #make and fit signal pdf for WW and WZ
            self.Make_signal_pdf(self.rrv_mass_lvj,'WW')
            self.Make_signal_pdf(self.rrv_mass_lvj,'WZ')

            #read, rename and write bkg pdfs and bkg rates
            fileInWs    = TFile.Open('Input/wwlvj_%s_HPV_workspace.root'%self.ch)
            w_bkg       = fileInWs.Get('workspace4limit_') 
            w_bkg.var('rrv_mass_j').setBins(self.nbins)

            path        ='%s/src/CombinedEWKAnalysis/CommonTools/data/anomalousCoupling'%os.environ["CMSSW_BASE"]

            for bkg in ['WJets','TTbar','STop','WW','WZ']:
                w_bkg.var('norm_%s_%s'%(bkg,self.ch)).Print()
                getattr(self.WS,'import')(w_bkg.var('norm_%s_%s'%(bkg,self.ch)))

            #import m_pruned and define ranges
            getattr(self.WS,'import')(w_bkg.var('rrv_mass_j'))
            self.WS.var('rrv_mass_j').setRange('sb_lo',40,65)
            self.WS.var('rrv_mass_j').setRange('sig',65,105)
            self.WS.var('rrv_mass_j').setRange('sb_hi',105,150)
            self.WS.var('rrv_mass_lvj').setRange(900,float(options.mWVCutoff))

            #bkg-pdfs have the format '[bkg-name]_mlvj_[region]_[ch]' or '[bkg-name]_mj_[region]_[ch]'

            #create a workspace for each component in each region
            for region in self.regions:
                self.WS2 = self.WS.Clone("w")        #temporary 
                set_mj        = RooArgSet(self.WS2.var('rrv_mass_j'))
                for bkg in self.bkgs:
                    #define global norm for whole mj spectrum
                    norm_var    = RooRealVar('normvar_%s_%s'%(bkg,self.ch),'normvar_%s_%s'%(bkg,self.ch),self.WS2.var("norm_%s_%s"%(bkg,self.ch)).getVal(),0,5e4)
                    norm_var.setConstant(kTRUE)
                    #define integral over region
                    reg_Int     = w_bkg.pdf('mj_%s_%s'%(bkg,self.ch)).createIntegral(set_mj,set_mj, region)
                    if bkg=='WJets':        #norm floating for WJets, integral depends on (floating) shape parameter
                        norm        = RooFormulaVar('%s_%s_%s_norm'%(bkg,region,self.ch),'%s_%s_%s_norm'%(bkg,region,self.ch),'@0*@1',RooArgList(reg_Int,norm_var))
                    else:#norm and integral fixed for rest
                        norm        = RooFormulaVar('%s_%s_%s_norm'%(bkg,region,self.ch),'%s_%s_%s_norm'%(bkg,region,self.ch),'%s*@0'%reg_Int.getVal(),RooArgList(norm_var))
                    if region == 'sig':
                        bkg_MWV     = w_bkg.pdf('%s_mlvj_sig_%s'%(bkg,self.ch))
                    else:#pdfs from the sb fit are fitted simultaneously in the lower and upper sb
                        bkg_MWV = w_bkg.pdf('%s_mlvj_sb_%s'%(bkg,self.ch)).clone('%s_mlvj_%s_%s'%(bkg,region,self.ch))
                    bkg_mj          = w_bkg.pdf('%s_mj_%s_%s'%(bkg,region,self.ch))
                    #make 2d pdf
                    bkg_2d_pdf      = RooProdPdf(bkg,bkg,RooArgList(bkg_MWV,bkg_mj))
                    bkg_MWV.Print();bkg_mj.Print();bkg_2d_pdf.Print();
                    norm.SetName(bkg_2d_pdf.GetName()+'_norm')#the normalization variable must have the corresponding pdf-name + _norm
                    self.Import_to_ws(self.WS2,[bkg_2d_pdf,norm],1)

                #signal function for WW and WZ in signal region and lower/upper sideband
                ##FIXME? signal shape is not explicitly evaluated in the sideband region since its contribution is assumed to be negligible there
                data_obs            = RooDataSet('data_obs','data_obs',w_bkg.data('dataset_2d_%s_%s'%(region,self.ch)),RooArgSet(self.WS2.var('rrv_mass_lvj'),self.WS2.var('mj_%s'%region)))
                getattr(self.WS2,'import')(data_obs)

                for VV in ['WW','WZ']:
                    pdf_atgc_mlvj_VV    = self.WS2.pdf('aTGC_model_%s_%s'%(self.ch,VV))
                    pdf_atgc_mj_VV      = w_bkg.pdf('%s_mj_%s_%s'%(VV,region,self.ch))
                    pdf_atgc_VV_2d      = RooProdPdf('aTGC_%s_%s_%s'%(VV,region,self.ch),'aTGC_%s_%s_%s'%(VV,region,self.ch),RooArgList(pdf_atgc_mlvj_VV,pdf_atgc_mj_VV))
                    #final normalization
                    norm_VV_reg         =self.WS2.function("%s_norm"%VV).Clone("%s_norm_%s_%s"%(VV,region,self.ch))
                    signal_norm_VV      = RooFormulaVar(pdf_atgc_VV_2d.GetName()+'_norm',pdf_atgc_VV_2d.GetName()+'_norm','@0*@1',RooArgList(self.WS2.function('normfactor_3d_%s_%s'%(self.ch,VV)),norm_VV_reg))
                    self.Import_to_ws(self.WS2,[pdf_atgc_VV_2d,signal_norm_VV],1)

                ##define which parameters are floating (also has to be done in the datacard)
                self.WS2.var("rrv_shift_ChiSq_WJets0_%s"%self.ch).setConstant(kFALSE)
                self.WS2.var("rrv_c_ChiSq_WJets0_%s"%self.ch).setConstant(kFALSE)
		#self.WS2.var("rrv_b0_Poly3_WJets0_%s"%self.ch).setConstant(kFALSE)
		#self.WS2.var("rrv_b1_Poly3_WJets0_%s"%self.ch).setConstant(kFALSE)
		#self.WS2.var("rrv_b2_Poly3_WJets0_%s"%self.ch).setConstant(kFALSE)
		#self.WS2.var("rrv_b3_Poly3_WJets0_%s"%self.ch).setConstant(kFALSE)
                self.WS2.var("normvar_WJets_%s"%self.ch).setConstant(kFALSE)
                if 'sb' in region:
                    self.WS2.var("rrv_c_ExpN_WJets0_sb_%s"%self.ch).setConstant(kFALSE)
                    self.WS2.var("rrv_n_ExpN_WJets0_sb_%s"%self.ch).setConstant(kFALSE)
                else:
                    self.WS2.var("Deco_WJets0_sim_%s_HPV_mlvj_13TeV_eig0"%self.ch).setConstant(kTRUE)
                    self.WS2.var("Deco_WJets0_sim_%s_HPV_mlvj_13TeV_eig1"%self.ch).setConstant(kTRUE)
                    self.WS2.var("Deco_WJets0_sim_%s_HPV_mlvj_13TeV_eig2"%self.ch).setConstant(kTRUE)
                    self.WS2.var("Deco_WJets0_sim_%s_HPV_mlvj_13TeV_eig3"%self.ch).setConstant(kTRUE)
                if region == 'sig':
                   self.WS2.var("Deco_TTbar_sig_%s_HPV_mlvj_13TeV_eig0"%self.ch).setConstant(kTRUE)
                   self.WS2.var("Deco_TTbar_sig_%s_HPV_mlvj_13TeV_eig1"%self.ch).setConstant(kTRUE)
                else:
                   self.WS2.var("Deco_TTbar_sb_%s_HPV_mlvj_13TeV_eig0"%self.ch).setConstant(kTRUE)
                   self.WS2.var("Deco_TTbar_sb_%s_HPV_mlvj_13TeV_eig1"%self.ch).setConstant(kTRUE)

                output        = TFile('WWWZ_%s_%s_ws.root'%(region,self.ch),'recreate')
                self.WS2.SetName('proc_WWWZ_%s_%s'%(region,self.ch))
                self.WS2.Write();
                output.Close()
                print 'Write to file ' + output.GetName()


            ##create the datacards for all regions
            self.Write_datacard(w_bkg,"sb_lo")
            self.Write_datacard(w_bkg,"sig")
            self.Write_datacard(w_bkg,"sb_hi")

            #make some plots
            if options.Make_plots:
                self.Make_plots(self.rrv_mass_lvj,'WW',self.fitresults)
                self.Make_plots(self.rrv_mass_lvj,'WZ',self.fitresults)
            
            for i in range(len(self.fitresults)):
                    self.fitresults[i].Print()

            if options.printatgc:
                self.wtmp.var('cwww').setVal(3.6)
                self.wtmp.var('ccw').setVal(4.5)
                self.wtmp.var('cb').setVal(0)
                print 'cwww and ccw positive:'
                for i in range(8):
                    print N_list.at(i).GetName() + ' : ' + str(N_list.at(i).getVal())
                self.wtmp.var('cwww').setVal(3.6)
                self.wtmp.var('ccw').setVal(0)
                self.wtmp.var('cb').setVal(20)
                print 'cwww and cb positive:'
                for i in range(8):
                    print N_list.at(i).GetName() + ' : ' + str(N_list.at(i).getVal())
                self.wtmp.var('cwww').setVal(0)
                self.wtmp.var('ccw').setVal(4.5)
                self.wtmp.var('cb').setVal(20)
                print 'ccw and cb positive:'
                for i in range(8):
                    print N_list.at(i).GetName() + ' : ' + str(N_list.at(i).getVal())

                #actual yields
                for i in range(3):
                    for j in range(3):
                        self.wtmp.var(self.POI[j]).setVal(0)
                    self.wtmp.var(self.POI[i]).setVal(self.PAR_MAX[self.POI[i]])
                    print channel + ' ' + self.POI[i] + ' : ' + str(w.var('rate_VV').getVal()*normfactor_3d.getVal())

                raw_input(channel)


###run code###

if __name__ == '__main__':
    if options.chan=='elmu':
        makeWS_el        = Prepare_workspace_4limit('el',900,float(options.mWVCutoff),float(options.binWidth))
        makeWS_el.Make_input()
        makeWS_mu        = Prepare_workspace_4limit('mu',900,float(options.mWVCutoff),float(options.binWidth))
        makeWS_mu.Make_input()
    else:
        makeWS        = Prepare_workspace_4limit(options.chan,900,float(options.mWVCutoff),float(options.binWidth))
        makeWS.Make_input()
    #combine the created datacards
    output_card_name = 'aC_WWWZ_simfit'
    cmd = 'combineCards.py aC_WWWZ_sig_el.txt aC_WWWZ_sig_mu.txt aC_WWWZ_sb_lo_el.txt aC_WWWZ_sb_lo_mu.txt aC_WWWZ_sb_hi_el.txt aC_WWWZ_sb_hi_mu.txt > %s.txt'%output_card_name
    #cmd = 'combineCards.py aC_WWWZ_sig_mu.txt aC_WWWZ_sb_lo_mu.txt aC_WWWZ_sb_hi_mu.txt > %s.txt'%output_card_name
    print cmd
    os.system(cmd)
    print 'generated Card : %s.txt'%output_card_name
