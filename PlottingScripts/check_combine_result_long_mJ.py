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
parser.add_option('-P','--POI', dest='poi', default='cwww:0')
parser.add_option('-n', dest='name', default='')
parser.add_option('-a', '--asimov', dest='asimovDataFile', default='none')
(options,args) = parser.parse_args()
gStyle.SetOptStat(0)
gStyle.SetOptTitle(0)


def plot(w,fitres,normset,ch,region):
    
    rrv_x       = w.var("mj_%s"%region)

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

    #w.var("cwww").setVal(0);w.var("ccw").setVal(0);w.var("cb").setVal(0);
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

    #draw_error_band(model, rrv_x.GetName(), rrv_model_norm_tmp, fitres.floatParsFinal(), w, p, kBlack, "F", 3013)
    pdfCurve=p.getCurve("WJets")
    uncBand=TGraphErrors(int((rrv_x.getMax()-rrv_x.getMin())/0.25))
    for bin in range(int((rrv_x.getMax()-rrv_x.getMin())/0.25)+1):
        x = rrv_x.getMin() + (bin*0.25)
        uncBand.SetPoint(bin, x, pdfCurve.interpolate(x))
        uncBand.SetPointError(bin, 0, pdfCurve.interpolate(x) * 1 * total_bkg_unc / model_norm_tmp)
    uncBand.SetFillStyle(3013)
    p.addObject(uncBand,"E4")

    ## Getting background final yields plus uncertainties
    #rrv_x.setRange("RangeSigLo",65,85)
    #rrv_x.setRange("RangeSigHi",85,105)
    #argSet=RooArgSet(rrv_x)

    #print "\nWW and WZ Region Yields"
    #print "=========================\n"
    #valueDibosonLo  = 0.0
    #uncDibosonLo    = 0.0
    #valueDibosonHi  = 0.0
    #uncDibosonHi    = 0.0
    #for process in ["WJets","TTbar","STop","WW","WZ"]:
    #    intSigLo  = bkg_pdfs[process].createIntegral(argSet,RooFit.NormSet(argSet),RooFit.Range("RangeSigLo"))
    #    value     = intSigLo.getVal() * bkg_norms[process].getVal()
    #    shapeUnc  = intSigLo.getPropagatedError(fitres) * value
    #    normUnc   = bkg_uncs[process] * intSigLo.getVal()
    #    unc       = TMath.Sqrt((shapeUnc*shapeUnc) + (normUnc*normUnc))
    #    print process + " in WW region: " + str(value) + " +- quad(" + str(shapeUnc) + ", " + str(normUnc) + ") = " + str(unc)
    #    if process == "WW" or process == "WZ":
    #        valueDibosonLo += value
    #        uncDibosonLo = TMath.Sqrt((uncDibosonLo*uncDibosonLo) + (unc*unc))
    
    #    intSigHi  = bkg_pdfs[process].createIntegral(argSet,RooFit.NormSet(argSet),RooFit.Range("RangeSigHi"))
    #    value     = intSigHi.getVal() * bkg_norms[process].getVal()
    #    shapeUnc  = intSigHi.getPropagatedError(fitres) * value
    #    normUnc   = bkg_uncs[process] * intSigHi.getVal()
    #    unc       = TMath.Sqrt((shapeUnc*shapeUnc) + (normUnc*normUnc))
    #    print process + " in WZ region: " + str(value) + " +- quad(" + str(shapeUnc) + ", " + str(normUnc) + ") = " + str(unc)
    #    if process == "WW" or process == "WZ":
    #        valueDibosonHi += value
    #        uncDibosonHi = TMath.Sqrt((uncDibosonHi*uncDibosonHi) + (unc*unc))

    #print "Diboson in WW region: " + str(valueDibosonLo) + " +- " + str(uncDibosonLo)
    #print "Diboson in WZ region: " + str(valueDibosonHi) + " +- " + str(uncDibosonHi)

    #raw_input("\nFinal yields.\n")

    w.var('cwww').setVal(cwwwtmp);w.var('ccw').setVal(ccwtmp);w.var('cb').setVal(cbtmp);
    #model.plotOn(p,RooFit.Name("WWWZ_atgc"),RooFit.Normalization(model_norm,RooAbsReal.NumEvent),RooFit.FillColor(colors["WW"]-6),RooFit.LineStyle(9),RooFit.LineWidth(2))
    #model.plotOn(p,RooFit.Name("WZ_atgc"),RooFit.Components("STop,WJets,TTbar,WZ"),RooFit.Normalization(model_norm,RooAbsReal.NumEvent),RooFit.FillColor(colors["WW"]-8),RooFit.DrawOption("F"))

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

def plot_mJ(canvas,xlo,xhi,reg,w,fitres,normset,ch,pads,medianLines,paveTexts,legends):


    canvas.cd()
    pad     = TPad(reg,reg,xlo,0.25,xhi,1.0)
    if reg=='sb_lo':
        pad2    = TPad(reg+"_pull",reg+"_pull",xlo,0.0,xhi,0.31)
        pad.SetMargin(0.25,0,0.1,0.1)
        pad2.SetMargin(0.25,0.0,0.4,0)
    elif reg=='sb_hi':
        pad2    = TPad(reg+"_pull",reg+"_pull",xlo,0.0,xhi,0.31)
        pad.SetMargin(0,0.1,0.1,0.1)
        pad2.SetMargin(0,0.1,0.4,0)
    else:
        pad2    = TPad(reg+"_pull",reg+"_pull",xlo-0.00445,0.0,xhi+0.001,0.31)
        pad.SetMargin(0,0,0.1,0.1)
        pad2.SetMargin(0.022,0,0.4,0)
    pad.Draw()
    pad2.Draw()

    pad.cd()
    # Main MJ plot
    p       = plot(w,fitres,normset,ch,reg)
    #p.GetXaxis().SetTitleSize(0.0575)
    #p.GetXaxis().SetTitleOffset(0.7)
    p.GetXaxis().SetTitleSize(0)
    p.GetXaxis().SetLabelSize(0)
    p.GetYaxis().SetRangeUser(0,1450)
    p.GetYaxis().SetLabelSize(0.05)
    p.GetYaxis().SetTitle('Events / 5 GeV')
    p.GetYaxis().SetTitleSize(0.08)
    p.GetYaxis().SetTitleOffset(1.1)
    if reg=='sb_lo':
        p.GetYaxis().SetRangeUser(0,1535)
        p.GetXaxis().SetTitle('')
    elif reg=='sb_hi':
        p.GetYaxis().SetTitle('')
        p.GetYaxis().SetLabelSize(0)
    else:
        p.GetXaxis().SetTitle('')
        p.GetYaxis().SetTitle('')
        p.GetYaxis().SetLabelSize(0)
    p.Draw()

    # Lumi text and channel
    if reg=='sb_lo':
       CMS_lumi.lumiTextSize=0.0
       #CMS_lumi.writeExtraText=True
       CMS_lumi.cmsTextSize=0.75
       CMS_lumi.relPosY    = -0.09
       CMS_lumi.relExtraDX = 0.3
       CMS_lumi.relExtraDY = 0.24
       CMS_lumi.CMS_lumi(pad,4,11)
    elif reg=='sb_hi':
       CMS_lumi.cmsTextSize=0.0
       #CMS_lumi.writeExtraText=False
       CMS_lumi.lumiTextSize=0.65
       CMS_lumi.lumiTextOffset=0.2
       CMS_lumi.CMS_lumi(pad,4,11)

    # Channel text
    if reg=='sig':
        pt = TPaveText(0.375,0.8,0.625,0.97, "blNDC")
        pt.SetFillStyle(0)
        pt.SetBorderSize(0)
        pt.SetTextAlign(23)
        pt.SetTextSize(0.09)
        if ch=="el":
            pt.AddText("Electron channel")
        else:
            pt.AddText("Muon channel")
        pt.Draw()
        paveTexts.append(pt)

    # Legend
    if reg=='sb_hi':
        legMJ=TLegend(0.625,0.39,0.875,0.875)
        legMJ.SetFillColor(0)
        legMJ.SetFillStyle(0)
        legMJ.SetBorderSize(0)
        legMJ.SetLineColor(0)
        legMJ.SetLineWidth(0)
        legMJ.SetLineStyle(0)
        legMJ.SetTextFont(42)
        if ch=='el':
            #legMJ.AddEntry(p.getObject(11),"CMS data, WV#rightarrow e#nuqq","P")
            legMJ.AddEntry(p.getObject(11),"Data","PE")
        else:
            #legMJ.AddEntry(p.getObject(11),"CMS data, WV#rightarrow #mu#nuqq","P")
            legMJ.AddEntry(p.getObject(11),"Data","PE")
        #legMJ.AddEntry(p.getObject(10),"Signal","L")
        legMJ.AddEntry(p.getObject(0),"W+jets","F")
        legMJ.AddEntry(p.getObject(1),"t#bar{t}","F")
        legMJ.AddEntry(p.getObject(4),"WW","F")
        legMJ.AddEntry(p.getObject(5),"WZ","F")
        legMJ.AddEntry(p.getObject(8),"Single top","F")
        legMJ.AddEntry(p.getObject(10),"Post-fit unc.","F")
        legMJ.Draw()
    else:
        legMJ=[]

    pad2.cd()

    # MJ pull histograms
    if reg=='sb_lo':
        pullhist    = p.pullHist("data","WJets")
        pullhist.SetMaximum(4.9)
        pullhist.SetMinimum(-5)
        pullhist.GetXaxis().SetRangeUser(40,65)
        pullhist.GetXaxis().SetNdivisions(505)
        pullhist.GetXaxis().SetLabelSize(0.14)
        pullhist.GetYaxis().SetNdivisions(7)
        pullhist.GetYaxis().SetTitle('#frac{Data-Fit}{#sigma_{Data}}')
        pullhist.GetYaxis().SetLabelSize(0.125)
        pullhist.GetYaxis().SetTitleSize(0.2)
        pullhist.GetYaxis().SetTitleOffset(0.35)
        pullhist.SetMarkerStyle(20)
        #pullhist.SetMarkerSize(1.5)
        pullhist.SetLineColor(kBlack)
        pullhist.SetMarkerColor(kBlack)
        pullhist.Draw("AP")
        medianLine = TLine(pullhist.GetXaxis().GetXmin(),0.,pullhist.GetXaxis().GetXmax(),0.); medianLine.SetLineWidth(1); medianLine.SetLineColor(kBlue); medianLine.Draw();
        pullhist.Draw("Psame")

    elif reg=='sig':
        pullhist    = p.pullHist("data","WJets")
        pullhist.SetMaximum(4.9)
        pullhist.SetMinimum(-5)
        pullhist.GetXaxis().SetRangeUser(65,105)
        pullhist.GetXaxis().SetLabelSize(0.14)
        pullhist.GetYaxis().SetNdivisions(7)
        pullhist.GetYaxis().SetLabelSize(0)
        pullhist.GetYaxis().SetTitleSize(0.2)
        pullhist.GetYaxis().SetTitleOffset(0.32)
        pullhist.SetMarkerStyle(20)
        #pullhist.SetMarkerSize(1.5)
        pullhist.SetLineColor(kBlack)
        pullhist.SetMarkerColor(kBlack)
        pullhist.Draw("AP")
        medianLine = TLine(65.,0.,105.,0.); medianLine.SetLineWidth(1); medianLine.SetLineColor(kBlue); medianLine.Draw();
        pullhist.Draw("Psame")
    else:
        pullhist    = p.pullHist("data","WJets")
        pullhist.SetMaximum(4.9)
        pullhist.SetMinimum(-5)
        pullhist.GetXaxis().SetRangeUser(105,150)
        pullhist.GetXaxis().SetLabelSize(0.14)
        pullhist.GetXaxis().SetTitle('m_{SD} (GeV)')
        pullhist.GetXaxis().SetTitleSize(0.2)
        pullhist.GetXaxis().SetTitleOffset(0.85)
        pullhist.GetYaxis().SetNdivisions(7)
        pullhist.GetYaxis().SetLabelSize(0)
        pullhist.GetYaxis().SetTitleSize(0.2)
        pullhist.GetYaxis().SetTitleOffset(0.32)
        pullhist.SetMarkerStyle(20)
        #pullhist.SetMarkerSize(1.5)
        pullhist.SetLineColor(kBlack)
        pullhist.SetMarkerColor(kBlack)
        pullhist.Draw("AP")
        medianLine = TLine(105.,0.,150.,0.); medianLine.SetLineWidth(1); medianLine.SetLineColor(kBlue); medianLine.Draw();
        pullhist.Draw("Psame")

    canvas.Update()

    pads.append(pad)
    pads.append(pad2)
    medianLines.append(medianLine)
    legends.append(legMJ)

def plot_all(w,ch="el"):
   
    canvas = TCanvas(ch,ch,2500,500)

    offset = 0.15
    # These things just have to be kept in memory so that ROOT does not make them disappear in the next loop
    pads = []
    medianLines = []
    paveTexts = []
    legendsMJ = []
    legendsMWV = []

    # MJ main and pull plots
    plot_mJ(canvas,0+offset-0.05,5/22.+offset/2,'sb_lo',w,fitres,normset,ch,pads,medianLines,paveTexts,legendsMJ)
    plot_mJ(canvas,5/22.+offset/2,13/22.-offset/2,'sig',w,fitres,normset,ch,pads,medianLines,paveTexts,legendsMJ)
    plot_mJ(canvas,13/22.-offset/2,1-offset,'sb_hi',w,fitres,normset,ch,pads,medianLines,paveTexts,legendsMJ)

    canvas.cd()
    canvas.Draw()
    canvas.Update()
    canvas.SaveAs('postfit_%s_mJ.eps'%ch)
    os.system('epstopdf postfit_%s_mJ.eps'%ch)
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

ROOT.gStyle.SetPaperSize(40,52)
#ROOT.gROOT.ProcessLine("Float_t* wid=new Float_t(0); Float_t* hei=new Float_t(0); gStyle->GetPaperSize(wid, hei); std::cout<<*wid<<"   "<<*hei<<std::endl;")
#ROOT.gROOT.SetBatch(kTRUE)
#plot_all(w,options.ch)

string = '{:>40} : {:>30} / {:>30}\n'.format('>>name<<','>>pre-fit<<','>>post-fit<<')
for i in range(fitparas.getSize()):
    prefit  = str(round(w.var(fitparas.at(i).GetName()).getVal(),6)) + ' +- ' + str(round(w.var(fitparas.at(i).GetName()).getError(),6))
    postfit = str(round(fitparas.at(i).getVal(),6)) + ' +- ' + str(round(fitparas.at(i).getError(),6))
    string += '{:>40} : {:>30} / {:>30}\n'.format(fitparas.at(i).GetName(),prefit,postfit)
for i in range(fitparas.getSize()):
    w.var(fitparas.at(i).GetName()).setVal(fitparas.at(i).getVal())


plot_all(w,options.ch)


print string

print 'Be careful, these are pre-fit values!!'
print 'expected events el: ' + str(w.var("normvar_WJets_el").getVal()+w.var("normvar_TTbar_el").getVal()+w.var("normvar_STop_el").getVal()+w.var("normvar_WW_el").getVal()+w.var("normvar_WZ_el").getVal())
print 'observed events el: ' + str(w.data("data_obs").sumEntries("CMS_channel==CMS_channel::ch1||CMS_channel==CMS_channel::ch3||CMS_channel==CMS_channel::ch5"))
print 'expected events mu: ' + str(w.var("normvar_WJets_mu").getVal()+w.var("normvar_TTbar_mu").getVal()+w.var("normvar_STop_mu").getVal()+w.var("normvar_WW_mu").getVal()+w.var("normvar_WZ_mu").getVal())
print 'observed events mu: ' + str(w.data("data_obs").sumEntries("CMS_channel==CMS_channel::ch2||CMS_channel==CMS_channel::ch4||CMS_channel==CMS_channel::ch6"))

print '\nThis is post-fit'
print 'Post-fit uncertainty example (TTbar ch1): ' + str(RooRealVar(normset.find("ch1/TTbar")).getError())

