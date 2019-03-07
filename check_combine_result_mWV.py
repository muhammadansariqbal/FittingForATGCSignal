from ROOT import  *
from array import array
from optparse import OptionParser
import math as math
import os
import CMS_lumi, tdrstyle

gSystem.Load("%s/lib/slc6_amd64_gcc481/libHiggsAnalysisCombinedLimit.so"%os.environ["CMSSW_BASE"])

gSystem.Load("PDFs/PdfDiagonalizer_cc.so")
gSystem.Load("PDFs/Util_cxx.so")
gSystem.Load("PDFs/hyperg_2F1_c.so")
gSystem.Load("PDFs/HWWLVJRooPdfs_cxx.so")
from ROOT import draw_error_band
from ROOT import RooPoly3Pdf, RooChiSqPdf, RooErfExpPdf, RooErfPowExpPdf, RooErfPowPdf, RooErfPow2Pdf, RooExpNPdf, RooAlpha4ExpNPdf, RooExpTailPdf, RooAlpha4ExpTailPdf, Roo2ExpPdf, draw_error_band_extendPdf


ch_nums = {"sig_el":"ch1","sb_lo_el":"ch3","sb_hi_el":"ch5","sig_mu":"ch2","sb_lo_mu":"ch4","sb_hi_mu":"ch6"}
mj_bins = {"sb_lo":5,"sig":8,"sb_hi":9}

parser        = OptionParser()
parser.add_option('-c', '--ch', dest='ch', default='el', help='channel, el or mu')
parser.add_option('-r', '--reg', dest='reg', default='sig')
parser.add_option('-P','--POI', dest='poi', default='cwww:0')
parser.add_option('-n', dest='name', default='')
parser.add_option('-a', '--asimov', dest='asimovDataFile', default='none')
(options,args) = parser.parse_args()
gStyle.SetOptStat(0)
gStyle.SetOptTitle(0)


def plot(w,fitres,normset,ch,region):
    
    rrv_x       = w.var("rrv_mass_lvj")

    p           = rrv_x.frame(RooFit.Bins(rrv_x.getBins()),RooFit.Title(region+"_"+ch))
    if ch=='el':
        p.GetYaxis().SetRangeUser(1e-4,100)
    elif ch=='mu':
        p.GetYaxis().SetRangeUser(1e-4,130)

    ch_num      = ch_nums[region+'_'+ch]
    bkgs        = ['WJets','TTbar','STop']
    if options.asimovDataFile=='none':
	data = w.data("data_obs")
    else:
	fileAsimov = TFile.Open(options.asimovDataFile)
	data = fileAsimov.Get("toys/toy_asimov")

    colors      = {'WJets':kGreen+1, 'TTbar':kOrange, 'STop':kBlue, 'WW':kRed, 'WZ':kRed+2, 'atgc':kMagenta}


    bkg_pdfs    = {}
    bkg_norms    = {}
    bkg_uncs    = {}
    for bkg in bkgs:
        bkg_pdfs[bkg]       = RooExtendPdf(bkg,bkg,w.pdf("shapeBkg_%s_%s"%(bkg,ch_num)),w.function("shapeBkg_%s_%s__norm"%(bkg,ch_num)))
        #bkg_norms[bkg]      = w.function("shapeBkg_%s_%s__norm"%(bkg,ch_num))
        bkg_norms[bkg]      = RooRealVar("norm_%s"%bkg,"norm_%s"%bkg,normset.getRealValue("%s/%s"%(ch_nums[region+"_"+ch],bkg)))
	#bkg_pdfs[bkg]       = RooExtendPdf(bkg,bkg,w.pdf("shapeBkg_%s_%s"%(bkg,ch_num)),bkg_norms[bkg])
	bkg_uncs[bkg]       = RooRealVar(normset.find("%s/%s"%(ch_nums[region+"_"+ch],bkg))).getError()

    bkg_pdfs["WW"]    = RooExtendPdf("WW","WW",w.pdf("shapeSig_aTGC_WW_%s_%s_%s"%(region,ch,ch_num)),w.function("shapeSig_aTGC_WW_%s_%s_%s__norm"%(region,ch,ch_num)))
    #bkg_norms["WW"]   = w.function("shapeSig_aTGC_WW_%s_%s_%s__norm"%(region,ch,ch_num))
    bkg_norms["WW"]   = RooRealVar("norm_%s"%bkg,"norm_%s"%bkg,normset.getRealValue("%s/aTGC_WW_%s_%s"%(ch_nums[region+"_"+ch],region,ch)))
    #bkg_pdfs["WW"]    = RooExtendPdf("WW","WW",w.pdf("shapeSig_aTGC_WW_%s_%s_%s"%(region,ch,ch_num)),bkg_norms["WW"])
    bkg_uncs["WW"]       = RooRealVar(normset.find("%s/aTGC_WW_%s_%s"%(ch_nums[region+"_"+ch],region,ch))).getError()

    bkg_pdfs["WZ"]    = RooExtendPdf("WZ","WZ",w.pdf("shapeSig_aTGC_WZ_%s_%s_%s"%(region,ch,ch_num)),w.function("shapeSig_aTGC_WZ_%s_%s_%s__norm"%(region,ch,ch_num)))
    #bkg_norms["WZ"]   = w.function("shapeSig_aTGC_WZ_%s_%s_%s__norm"%(region,ch,ch_num))
    bkg_norms["WZ"]   = RooRealVar("norm_%s"%bkg,"norm_%s"%bkg,normset.getRealValue("%s/aTGC_WZ_%s_%s"%(ch_nums[region+"_"+ch],region,ch)))
    #bkg_pdfs["WZ"]    = RooExtendPdf("WZ","WZ",w.pdf("shapeSig_aTGC_WZ_%s_%s_%s"%(region,ch,ch_num)),bkg_norms["WZ"])
    bkg_uncs["WZ"]       = RooRealVar(normset.find("%s/aTGC_WZ_%s_%s"%(ch_nums[region+"_"+ch],region,ch))).getError()

    # Print post-fit uncertainties
    total_bkg_unc   = TMath.Sqrt((bkg_uncs['WJets']*bkg_uncs['WJets']) + (bkg_uncs['TTbar']*bkg_uncs['TTbar']) + (bkg_uncs['STop']*bkg_uncs['STop']) + (bkg_uncs['WW']*bkg_uncs['WW']) + (bkg_uncs['WZ']*bkg_uncs['WZ']))
    for unc in bkg_uncs:
        print str(unc) + 'unc: ' + str(bkg_uncs[unc])
    print 'Total post-fit uncertainty in region: ' + str(total_bkg_unc)
    #raw_input("====================")

    #print("WJets norm from pdf: ", bkg_norms["WJets"].getVal(), "WJets norm from norm_fit_s: ", normset.getRealValue("%s/%s"%(ch_nums[region+"_"+ch],"WJets")))
    #print("TTbar norm from pdf: ", bkg_norms["TTbar"].getVal(), "TTbar norm from norm_fit_s: ", normset.getRealValue("%s/%s"%(ch_nums[region+"_"+ch],"TTbar")))
    #print("STop norm from pdf: ", bkg_norms["STop"].getVal(), "STop norm from norm_fit_s: ", normset.getRealValue("%s/%s"%(ch_nums[region+"_"+ch],"STop")))
    #print("WW norm from pdf: ", bkg_norms["WW"].getVal(), "WW norm from norm_fit_s: ", normset.getRealValue("%s/aTGC_WW_%s_%s"%(ch_nums[region+"_"+ch],region,ch)))
    #print("WZ norm from pdf: ", bkg_norms["WZ"].getVal(), "WZ norm from norm_fit_s: ", normset.getRealValue("%s/aTGC_WZ_%s_%s"%(ch_nums[region+"_"+ch],region,ch)))
    #raw_input("Debugging normalizations.")

    w.var(options.poi.split(':')[0]).setVal(float(options.poi.split(':')[1]))

    #cwwwtmp=w.var('cwww').getVal();ccwtmp=w.var('ccw').getVal();cbtmp=w.var('cb').getVal();
    cwwwtmp=1.59; ccwtmp=0; cbtmp=0;

    model   = RooAddPdf("model","model",RooArgList(bkg_pdfs["WW"],bkg_pdfs["WZ"],bkg_pdfs["TTbar"],bkg_pdfs["STop"],bkg_pdfs["WJets"]))
    model_norm  = float(bkg_norms["WJets"].getVal()+bkg_norms["STop"].getVal()+bkg_norms["TTbar"].getVal()+bkg_norms["WW"].getVal()+bkg_norms["WZ"].getVal())
    rrv_model_norm = RooRealVar("rrv_model_norm","rrv_model_norm",model_norm)

    w.var("cwww").setVal(0);w.var("ccw").setVal(0);w.var("cb").setVal(0);
    model_norm_tmp = float(bkg_norms["WJets"].getVal()+bkg_norms["STop"].getVal()+bkg_norms["TTbar"].getVal()+bkg_norms["WW"].getVal()+bkg_norms["WZ"].getVal())
    model.plotOn(p,RooFit.Name("WJets"),RooFit.Components("STop,WJets,TTbar,WW,WZ"),RooFit.Normalization(model_norm_tmp,RooAbsReal.NumEvent),RooFit.FillColor(colors["WJets"]),RooFit.LineColor(kBlack),RooFit.LineWidth(1),RooFit.DrawOption("F"))
    model.plotOn(p,RooFit.Name("TTbar"),RooFit.Components("STop,TTbar,WW,WZ"),RooFit.Normalization(model_norm_tmp,RooAbsReal.NumEvent),RooFit.FillColor(kOrange),RooFit.LineColor(kBlack),RooFit.LineWidth(1),RooFit.DrawOption("F"))
    model.plotOn(p,RooFit.Name("TTbar_line"),RooFit.Components("STop,TTbar,WW,WZ"),RooFit.Normalization(model_norm_tmp,RooAbsReal.NumEvent),RooFit.LineColor(kBlack),RooFit.LineWidth(1))
    model.plotOn(p,RooFit.Name("WJets_line"),RooFit.Components("STop,WJets,TTbar,WW,WZ"),RooFit.Normalization(model_norm_tmp,RooAbsReal.NumEvent),RooFit.LineColor(kBlack),RooFit.LineWidth(1))
    model.plotOn(p,RooFit.Name("WWSM"),RooFit.Components("WW,WZ,STop"),RooFit.Normalization(model_norm_tmp,RooAbsReal.NumEvent),RooFit.FillColor(colors["WW"]),RooFit.LineColor(kBlack),RooFit.LineWidth(1),RooFit.DrawOption("F"))
    model.plotOn(p,RooFit.Name("WZSM"),RooFit.Components("WZ,STop"),RooFit.Normalization(model_norm_tmp,RooAbsReal.NumEvent),RooFit.FillColor(colors["WZ"]),RooFit.LineColor(kBlack),RooFit.LineWidth(1),RooFit.DrawOption("F"))
    model.plotOn(p,RooFit.Name("WW_line"),RooFit.Components("WW,WZ,STop"),RooFit.Normalization(model_norm_tmp,RooAbsReal.NumEvent),RooFit.LineColor(kBlack),RooFit.LineWidth(1))
    model.plotOn(p,RooFit.Name("WZ_line"),RooFit.Components("WZ,STop"),RooFit.Normalization(model_norm_tmp,RooAbsReal.NumEvent),RooFit.LineColor(kBlack),RooFit.LineWidth(1))
    model.plotOn(p,RooFit.Name("STop"),RooFit.Components("STop"),RooFit.Normalization(model_norm,RooAbsReal.NumEvent),RooFit.FillColor(colors["STop"]),RooFit.LineColor(kBlack),RooFit.LineWidth(1),RooFit.DrawOption("F"))
    model.plotOn(p,RooFit.Name("STop_line"),RooFit.Components("STop"),RooFit.Normalization(model_norm,RooAbsReal.NumEvent),RooFit.LineColor(kBlack),RooFit.LineWidth(1))

    model.plotOn(p,RooFit.Name("uncBand"),RooFit.Components("STop,WJets,TTbar,WW,WZ"),RooFit.Normalization(model_norm_tmp,RooAbsReal.NumEvent),RooFit.FillColor(kBlack),RooFit.LineColor(kBlack),RooFit.LineWidth(1),RooFit.DrawOption("F"),RooFit.FillStyle(3013),RooFit.VisualizeError(fitres))

    ## Getting number of events in different bins for different backgrounds
    #string = '{:^20}   {:^20}   {:^20}   {:^30}   {:^30}\n'.format('Bin','W+jets integral','ttbar integral', 'W+jets cumulative integral','ttbar cumulative integral')
    #string += '{:^20}   {:^20}   {:^20}   {:^30}   {:^30}\n'.format('===','===============','==============', '==========================','=========================')
    ## First bin print 'Range: 900-1000'
    #rrv_x.setRange("ReducedRange",900,1000)
    #argSet=RooArgSet(rrv_x)
    #cumulIntWJets = bkg_pdfs["WJets"].createIntegral(argSet,RooFit.NormSet(argSet),RooFit.Range("ReducedRange")).getVal() * bkg_norms["WJets"].getVal()
    #cumulIntTTbar = bkg_pdfs["TTbar"].createIntegral(argSet,RooFit.NormSet(argSet),RooFit.Range("ReducedRange")).getVal() * bkg_norms["TTbar"].getVal()
    #string += '{:^20}   {:^20}   {:^20}   {:^30}   {:^30}\n'.format('900-1000',str(round(cumulIntWJets,2)),str(round(cumulIntTTbar,2)),str(round(cumulIntWJets,2)),str(round(cumulIntTTbar,2)))
    ## Others bin-width 500
    #for rangeStart in range(1000,4500,500):
    #    rrv_x.setRange("ReducedRange",rangeStart,rangeStart+500)
    #    argSet=RooArgSet(rrv_x)
    #    intWJets = bkg_pdfs["WJets"].createIntegral(argSet,RooFit.NormSet(argSet),RooFit.Range("ReducedRange")).getVal() * bkg_norms["WJets"].getVal()
    #    intTTbar = bkg_pdfs["TTbar"].createIntegral(argSet,RooFit.NormSet(argSet),RooFit.Range("ReducedRange")).getVal() * bkg_norms["TTbar"].getVal()
    #    cumulIntWJets += intWJets
    #    cumulIntTTbar += intTTbar
    #    string += '{:^20}   {:^20}   {:^20}   {:^30}   {:^30}\n'.format(str(rangeStart)+'-'+str(rangeStart+500),str(round(intWJets,2)),str(round(intTTbar,2)),str(round(cumulIntWJets,2)),str(round(cumulIntTTbar,2)))
    #print string
    #raw_input("Comparison of background events.")

    w.var('cwww').setVal(cwwwtmp);w.var('ccw').setVal(ccwtmp);w.var('cb').setVal(cbtmp);
    #if region=='sig':
    model.plotOn(p,RooFit.Name("WWWZ_atgc"),RooFit.Normalization(model_norm,RooAbsReal.NumEvent),RooFit.FillColor(colors["WW"]-6),RooFit.LineStyle(9),RooFit.LineWidth(2))
    #model.plotOn(p,RooFit.Name("WZ_atgc"),RooFit.Components("STop,WJets,TTbar,WZ"),RooFit.Normalization(model_norm,RooAbsReal.NumEvent),RooFit.FillColor(colors["WW"]-8),RooFit.LineStyle(9),RooFit.LineWidth(2))

    data_histo   = data.binnedClone("data","data").createHistogram("data",rrv_x,RooFit.Cut("CMS_channel==CMS_channel::%s"%ch_num))
    data_histo.Print()
    data_plot    = RooHist(data_histo,rrv_x.getBinWidth(0))
    data_plot.SetMarkerStyle(20)
    #data_plot.SetMarkerSize(1.5)
    alpha        = 1-0.6827
    # This is data in main plots
    for iPoint in range(data_plot.GetN()):
        N = data_plot.GetY()[iPoint]
        #print "x: " + str(data_plot.GetX()[iPoint]) + "   N: " + str(N)
        if N==0 :
            L = 0
        else:
            L = Math.gamma_quantile(alpha/2.,N,1.)
        U = Math.gamma_quantile_c(alpha/2,N+1,1.)
        data_plot.SetPointEYlow(iPoint, N-L)
        data_plot.SetPointEYhigh(iPoint, U-N)
        data_plot.SetPointEXlow(iPoint,0)
        data_plot.SetPointEXhigh(iPoint,0)
    data_plot.SetName('data')
    #raw_input("Data histogram")

    p.addPlotable(data_plot,"PE")

    return p

def plot_all(w,ch="el",reg='sig'):
   
    # These things just have to be kept in memory so that ROOT does not make them disappear in the next loop
    pads = []
    medianLines = []
    paveTexts = []
    legendsMWV = []

    # This is for median line in MWV pull plot
    ratio_style     = TH1D('ratio_style','ratio_style',36,900,4500)
    ratio_style.SetMaximum(2.9)
    ratio_style.SetMinimum(-3)
    ratio_style.GetXaxis().SetTitle('m_{WV} (GeV)')
    ratio_style.GetXaxis().SetTitleSize(0.2)
    ratio_style.GetXaxis().SetTitleOffset(0.75)
    ratio_style.GetXaxis().SetLabelSize(0.125)
    ratio_style.GetYaxis().SetNdivisions(7)
    ratio_style.GetYaxis().SetTitle('#frac{Data-Fit}{#sigma_{Data}}  ')
    ratio_style.GetYaxis().SetLabelSize(0.125)
    ratio_style.GetYaxis().SetTitleSize(0.16)
    ratio_style.GetYaxis().SetTitleOffset(0.25)
    ratio_style.SetMarkerStyle(20)
    #ratio_style.SetMarkerSize(1.5)

    canvas = TCanvas(ch,ch,800,640)
    pad1        = TPad('pad','pad',0.,0.26,1.,1.)
    pad_pull    = TPad('pad_pull','pad_pull',0.,0.,1.,0.323)
    pads.append(pad1)
    pads.append(pad_pull)
    # Main MWV plot
    p=plot(w,fitres,normset,ch,reg)
    p.GetXaxis().SetTitleSize(0)
    p.GetXaxis().SetLabelSize(0)
    p.GetYaxis().SetRangeUser(7e-2,5e3)
    p.GetYaxis().SetTitle('Events / 100 GeV')
    p.GetYaxis().SetTitleSize(0.08)
    p.GetYaxis().SetTitleOffset(0.59)
    p.GetYaxis().SetLabelSize(0.06)

    canvas.cd()
    pad1.Draw()
    pad_pull.Draw()

    pad1.cd()
    pad1.SetLogy()
    pad1.SetTicky()
    pad1.SetBottomMargin(0.1)
    p.Draw()

    pad_pull.cd()
    # MWV pull plot
    pad_pull.SetTopMargin(0)
    pad_pull.SetBottomMargin(0.35)
    pullhist = p.pullHist("data","WJets")
    ratio_style.Draw("")
    pullhist.SetLineColor(kBlack)
    pullhist.SetMarkerStyle(20)
    #pullhist.SetMarkerSize(1.5)
    pullhist.SetMarkerColor(kBlack)
    pullhist.Draw("SAME PE")
    #for i in range(pullhist.GetN()):
        #print pullhist.GetY()[i]
        #raw_input("Printed pull hist")
    
    # Lumi text and channel
    pad1.cd()
    #if reg=='sb_lo':
    CMS_lumi.lumiTextSize=0.0
    #CMS_lumi.writeExtraText=True
    CMS_lumi.cmsTextSize=0.75
    CMS_lumi.relPosY    = -0.09
    CMS_lumi.relExtraDX = 0.2
    CMS_lumi.relExtraDY = 0.24
    CMS_lumi.CMS_lumi(pad1,4,11)
    #elif reg=='sb_hi':
    CMS_lumi.cmsTextSize=0.0
    CMS_lumi.writeExtraText=False
    CMS_lumi.lumiTextSize=0.65
    CMS_lumi.lumiTextOffset=0.2
    CMS_lumi.CMS_lumi(pad1,4,11)
    #else:
    pt2 = TPaveText(0.125,0.8,0.325,0.975, "blNDC")
    pt2.SetFillStyle(0)
    pt2.SetBorderSize(0)
    pt2.SetTextAlign(13)
    pt2.SetTextSize(0.08)
    if ch=="el":
        pt2.AddText("Electron channel")
    else:
        pt2.AddText("Muon channel")
    pt2.Draw()
    paveTexts.append(pt2)

    # Legend
    legMWV=TLegend(0.525,0.4,0.88,0.875)
    legMWV.SetFillColor(0)
    legMWV.SetFillStyle(0)
    legMWV.SetBorderSize(0)
    legMWV.SetLineColor(0)
    legMWV.SetLineWidth(0)
    legMWV.SetLineStyle(0)
    legMWV.SetTextFont(42)
    if ch=='el':
        #legMWV.AddEntry(p.getObject(12),"CMS data, WV#rightarrow e#nuqq","P")
        legMWV.AddEntry(p.getObject(12),"Data","PE")
    else:
        #legMWV.AddEntry(p.getObject(12),"CMS data, WV#rightarrow #mu#nuqq","P")
        legMWV.AddEntry(p.getObject(12),"Data","PE")
    legMWV.AddEntry(p.getObject(11),"Signal c_{WWW}/#Lambda^{2}=1.59 TeV^{-2}","L")
    legMWV.AddEntry(p.getObject(0),"W+jets","F")
    legMWV.AddEntry(p.getObject(1),"t#bar{t}","F")
    legMWV.AddEntry(p.getObject(4),"WW","F")
    legMWV.AddEntry(p.getObject(5),"WZ","F")
    legMWV.AddEntry(p.getObject(8),"Single top","F")
    legMWV.AddEntry(p.getObject(10),"Post-fit unc.","F")
    legMWV.Draw()
    legendsMWV.append(legMWV)

    canvas.Update()
    canvas.SaveAs('postfit_%s_mWV_%s.pdf'%(ch,reg))
    raw_input(ch)

    for i in pads:
        i.Delete()



fileIn      = TFile.Open("workspace_simfit.root")
w           = fileIn.Get("w")
fileIn.Close()
fileIn      = TFile.Open("mlfit%s.root"%options.name)
fitres      = fileIn.Get("fit_s")
normset     = fileIn.Get("norm_fit_s")
fileIn.Close()

fitparas    = fitres.floatParsFinal()

#ROOT.gStyle.SetPaperSize(40,52)
#ROOT.gROOT.ProcessLine("Float_t* wid=new Float_t(0); Float_t* hei=new Float_t(0); gStyle->GetPaperSize(wid, hei); std::cout<<*wid<<"   "<<*hei<<std::endl;")
#ROOT.gROOT.SetBatch(kTRUE)
#plot_all(w,options.ch,options.reg)

string = '{:>40} : {:>30} / {:>30}\n'.format('>>name<<','>>pre-fit<<','>>post-fit<<')
for i in range(fitparas.getSize()):
    prefit  = str(round(w.var(fitparas.at(i).GetName()).getVal(),6)) + ' +- ' + str(round(w.var(fitparas.at(i).GetName()).getError(),6))
    postfit = str(round(fitparas.at(i).getVal(),6)) + ' +- ' + str(round(fitparas.at(i).getError(),6))
    string += '{:>40} : {:>30} / {:>30}\n'.format(fitparas.at(i).GetName(),prefit,postfit)
for i in range(fitparas.getSize()):
    w.var(fitparas.at(i).GetName()).setVal(fitparas.at(i).getVal())


plot_all(w,options.ch,options.reg)


print string

print 'Be careful, these are pre-fit values!!'
print 'expected events el: ' + str(w.var("normvar_WJets_el").getVal()+w.var("normvar_TTbar_el").getVal()+w.var("normvar_STop_el").getVal()+w.var("normvar_WW_el").getVal()+w.var("normvar_WZ_el").getVal())
print 'observed events el: ' + str(w.data("data_obs").sumEntries("CMS_channel==CMS_channel::ch1||CMS_channel==CMS_channel::ch3||CMS_channel==CMS_channel::ch5"))
print 'expected events mu: ' + str(w.var("normvar_WJets_mu").getVal()+w.var("normvar_TTbar_mu").getVal()+w.var("normvar_STop_mu").getVal()+w.var("normvar_WW_mu").getVal()+w.var("normvar_WZ_mu").getVal())
print 'observed events mu: ' + str(w.data("data_obs").sumEntries("CMS_channel==CMS_channel::ch2||CMS_channel==CMS_channel::ch4||CMS_channel==CMS_channel::ch6"))

print '\nThis is post-fit'
print 'Post-fit uncertainty example (TTbar ch1): ' + str(RooRealVar(normset.find("ch1/TTbar")).getError())

