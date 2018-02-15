from ROOT import  *
from array import array
from optparse import OptionParser
import math as math
import os

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


def plot(w,fitres,normset,spectrum,ch,region):
    
    if spectrum == "mj":
        rrv_x       = w.var("mj_%s"%region)
    elif spectrum == "mlvj":
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
    for bkg in bkgs:
        bkg_pdfs[bkg]       = RooExtendPdf(bkg,bkg,w.pdf("shapeBkg_%s_%s"%(bkg,ch_num)),w.function("shapeBkg_%s_%s__norm"%(bkg,ch_num)))
        bkg_norms[bkg]      = w.function("shapeBkg_%s_%s__norm"%(bkg,ch_num))
        #bkg_norms[bkg]      = RooRealVar("norm_%s"%bkg,"norm_%s"%bkg,normset.getRealValue("%s/%s"%(ch_nums[region+"_"+ch],bkg)))
	#bkg_pdfs[bkg]       = RooExtendPdf(bkg,bkg,w.pdf("shapeBkg_%s_%s"%(bkg,ch_num)),bkg_norms[bkg])
    bkg_pdfs["WW"]    = RooExtendPdf("WW","WW",w.pdf("shapeSig_aTGC_WW_%s_%s_%s"%(region,ch,ch_num)),w.function("shapeSig_aTGC_WW_%s_%s_%s__norm"%(region,ch,ch_num)))
    bkg_norms["WW"]   = w.function("shapeSig_aTGC_WW_%s_%s_%s__norm"%(region,ch,ch_num))
    #bkg_norms["WW"]   = RooRealVar("norm_%s"%bkg,"norm_%s"%bkg,normset.getRealValue("%s/aTGC_WW_%s_%s"%(ch_nums[region+"_"+ch],region,ch)))
    #bkg_pdfs["WW"]    = RooExtendPdf("WW","WW",w.pdf("shapeSig_aTGC_WW_%s_%s_%s"%(region,ch,ch_num)),bkg_norms["WW"])
    bkg_pdfs["WZ"]    = RooExtendPdf("WZ","WZ",w.pdf("shapeSig_aTGC_WZ_%s_%s_%s"%(region,ch,ch_num)),w.function("shapeSig_aTGC_WZ_%s_%s_%s__norm"%(region,ch,ch_num)))
    bkg_norms["WZ"]   = w.function("shapeSig_aTGC_WZ_%s_%s_%s__norm"%(region,ch,ch_num))
    #bkg_norms["WZ"]   = RooRealVar("norm_%s"%bkg,"norm_%s"%bkg,normset.getRealValue("%s/aTGC_WZ_%s_%s"%(ch_nums[region+"_"+ch],region,ch)))
    #bkg_pdfs["WZ"]    = RooExtendPdf("WZ","WZ",w.pdf("shapeSig_aTGC_WZ_%s_%s_%s"%(region,ch,ch_num)),bkg_norms["WZ"])

    #print("WJets norm from pdf: ", bkg_norms["WJets"].getVal(), "WJets norm from norm_fit_s: ", normset.getRealValue("%s/%s"%(ch_nums[region+"_"+ch],"WJets")))
    #print("TTbar norm from pdf: ", bkg_norms["TTbar"].getVal(), "TTbar norm from norm_fit_s: ", normset.getRealValue("%s/%s"%(ch_nums[region+"_"+ch],"TTbar")))
    #print("STop norm from pdf: ", bkg_norms["STop"].getVal(), "STop norm from norm_fit_s: ", normset.getRealValue("%s/%s"%(ch_nums[region+"_"+ch],"STop")))
    #print("WW norm from pdf: ", bkg_norms["WW"].getVal(), "WW norm from norm_fit_s: ", normset.getRealValue("%s/aTGC_WW_%s_%s"%(ch_nums[region+"_"+ch],region,ch)))
    #print("WZ norm from pdf: ", bkg_norms["WZ"].getVal(), "WZ norm from norm_fit_s: ", normset.getRealValue("%s/aTGC_WZ_%s_%s"%(ch_nums[region+"_"+ch],region,ch)))
    #raw_input("Debugging normalizations.")

    w.var(options.poi.split(':')[0]).setVal(float(options.poi.split(':')[1]))

    cwwwtmp=w.var('cwww').getVal();ccwtmp=w.var('ccw').getVal();cbtmp=w.var('cb').getVal()

    model   = RooAddPdf("model","model",RooArgList(bkg_pdfs["WW"],bkg_pdfs["WZ"],bkg_pdfs["TTbar"],bkg_pdfs["STop"],bkg_pdfs["WJets"]))
    model_norm  = float(bkg_norms["WJets"].getVal()+bkg_norms["STop"].getVal()+bkg_norms["TTbar"].getVal()+bkg_norms["WW"].getVal()+bkg_norms["WZ"].getVal())
    rrv_model_norm = RooRealVar("rrv_model_norm","rrv_model_norm",model_norm)

    if spectrum == "mj":
        #model.plotOn(p,RooFit.Name("WWWZ_atgc"),RooFit.Normalization(model_norm,RooAbsReal.NumEvent),RooFit.FillColor(colors["WW"]-8),RooFit.DrawOption("F"))
        #model.plotOn(p,RooFit.Name("WW_atgc"),RooFit.Components("STop,WJets,TTbar,WW"),RooFit.Normalization(model_norm,RooAbsReal.NumEvent),RooFit.FillColor(colors["WW"]-6),RooFit.DrawOption("F"))

        #w.var('cwww').setVal(0);w.var('ccw').setVal(0);w.var('cb').setVal(0);
        model_norm_tmp = float(bkg_norms["WJets"].getVal()+bkg_norms["STop"].getVal()+bkg_norms["TTbar"].getVal()+bkg_norms["WW"].getVal()+bkg_norms["WZ"].getVal())
        model.plotOn(p,RooFit.Name("WZSM"),RooFit.Normalization(model_norm_tmp,RooAbsReal.NumEvent),RooFit.FillColor(colors["WZ"]),RooFit.DrawOption("F"))
        model.plotOn(p,RooFit.Name("WWSM"),RooFit.Components("WJets,TTbar,STop,WW"),RooFit.Normalization(model_norm_tmp,RooAbsReal.NumEvent),RooFit.FillColor(colors["WW"]),RooFit.DrawOption("F"))
        model.plotOn(p,RooFit.Name("WZ_line"),RooFit.Normalization(model_norm_tmp,RooAbsReal.NumEvent),RooFit.LineColor(kBlack),RooFit.LineWidth(1))
        model.plotOn(p,RooFit.Name("WW_line"),RooFit.Components("WJets,TTbar,STop,WW"),RooFit.Normalization(model_norm_tmp,RooAbsReal.NumEvent),RooFit.LineColor(kBlack),RooFit.LineWidth(1))
        w.var('cwww').setVal(cwwwtmp);w.var('ccw').setVal(ccwtmp);w.var('cb').setVal(cbtmp);


        model.plotOn(p,RooFit.Name("TTbar"),RooFit.Components("WJets,STop,TTbar"),RooFit.Normalization(model_norm,RooAbsReal.NumEvent),RooFit.FillColor(colors["TTbar"]),RooFit.DrawOption("F"))
        model.plotOn(p,RooFit.Name("STop"),RooFit.Components("WJets,STop"),RooFit.Normalization(model_norm,RooAbsReal.NumEvent),RooFit.FillColor(colors["STop"]),RooFit.DrawOption("F"))
        model.plotOn(p,RooFit.Name("WJets"),RooFit.Components("WJets"),RooFit.Normalization(model_norm,RooAbsReal.NumEvent),RooFit.FillColor(colors["WJets"]),RooFit.DrawOption("F"))

        model.plotOn(p,RooFit.Name("TTbar_line"),RooFit.Components("WJets,STop,TTbar"),RooFit.Normalization(model_norm,RooAbsReal.NumEvent),RooFit.LineColor(kBlack),RooFit.LineWidth(1))
        model.plotOn(p,RooFit.Name("STop_line"),RooFit.Components("WJets,STop"),RooFit.Normalization(model_norm,RooAbsReal.NumEvent),RooFit.LineColor(kBlack),RooFit.LineWidth(2))
        model.plotOn(p,RooFit.Name("WJets_line"),RooFit.Components("WJets"),RooFit.Normalization(model_norm,RooAbsReal.NumEvent),RooFit.LineColor(kBlack),RooFit.LineWidth(2))

    elif spectrum == "mlvj":
        #model.plotOn(p,RooFit.Name("WWWZ_atgc"),RooFit.Normalization(model_norm,RooAbsReal.NumEvent),RooFit.FillColor(colors["WW"]-6),RooFit.DrawOption("F"))
        #model.plotOn(p,RooFit.Name("WZ_atgc"),RooFit.Components("STop,WJets,TTbar,WZ"),RooFit.Normalization(model_norm,RooAbsReal.NumEvent),RooFit.FillColor(colors["WW"]-8),RooFit.DrawOption("F"))

        #w.var("cwww").setVal(0);w.var("ccw").setVal(0);w.var("cb").setVal(0);
        model_norm_tmp = float(bkg_norms["WJets"].getVal()+bkg_norms["STop"].getVal()+bkg_norms["TTbar"].getVal()+bkg_norms["WW"].getVal()+bkg_norms["WZ"].getVal())
        model.plotOn(p,RooFit.Name("WJets"),RooFit.Components("STop,WJets,TTbar,WW,WZ"),RooFit.Normalization(model_norm_tmp,RooAbsReal.NumEvent),RooFit.FillColor(colors["WJets"]),RooFit.DrawOption("F"))
        model.plotOn(p,RooFit.Name("TTbar"),RooFit.Components("STop,TTbar,WW,WZ"),RooFit.Normalization(model_norm_tmp,RooAbsReal.NumEvent),RooFit.FillColor(kOrange),RooFit.DrawOption("F"))
        model.plotOn(p,RooFit.Name("TTbar_line"),RooFit.Components("STop,TTbar,WW,WZ"),RooFit.Normalization(model_norm_tmp,RooAbsReal.NumEvent),RooFit.LineColor(kBlack),RooFit.LineWidth(1))
        model.plotOn(p,RooFit.Name("WJets_line"),RooFit.Components("STop,WJets,TTbar,WW,WZ"),RooFit.Normalization(model_norm_tmp,RooAbsReal.NumEvent),RooFit.LineColor(kBlack),RooFit.LineWidth(1))
        model.plotOn(p,RooFit.Name("WWSM"),RooFit.Components("WW,WZ,STop"),RooFit.Normalization(model_norm_tmp,RooAbsReal.NumEvent),RooFit.FillColor(colors["WW"]),RooFit.DrawOption("F"))
        model.plotOn(p,RooFit.Name("WZSM"),RooFit.Components("WZ,STop"),RooFit.Normalization(model_norm_tmp,RooAbsReal.NumEvent),RooFit.FillColor(colors["WZ"]),RooFit.DrawOption("F"))
        model.plotOn(p,RooFit.Name("WW_line"),RooFit.Components("WW,WZ,STop"),RooFit.Normalization(model_norm_tmp,RooAbsReal.NumEvent),RooFit.LineColor(kBlack),RooFit.LineWidth(1))
        model.plotOn(p,RooFit.Name("WZ_line"),RooFit.Components("WZ,STop"),RooFit.Normalization(model_norm_tmp,RooAbsReal.NumEvent),RooFit.LineColor(kBlack),RooFit.LineWidth(1))
        w.var('cwww').setVal(cwwwtmp);w.var('ccw').setVal(ccwtmp);w.var('cb').setVal(cbtmp);

        model.plotOn(p,RooFit.Name("STop"),RooFit.Components("STop"),RooFit.Normalization(model_norm,RooAbsReal.NumEvent),RooFit.FillColor(colors["STop"]),RooFit.DrawOption("F"))
        model.plotOn(p,RooFit.Name("STop_line"),RooFit.Components("STop"),RooFit.Normalization(model_norm,RooAbsReal.NumEvent),RooFit.LineColor(kBlack),RooFit.LineWidth(1))

    data_histo   = data.binnedClone("data","data").createHistogram("data",rrv_x,RooFit.Cut("CMS_channel==CMS_channel::%s"%ch_num))
    data_histo.Print()
    data_plot    = RooHist(data_histo,rrv_x.getBinWidth(0))
    data_plot.SetMarkerStyle(20)
    data_plot.SetMarkerSize(1)
    alpha        = 1-0.6827
    for iPoint in range(data_plot.GetN()):
        N = data_plot.GetY()[iPoint]
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
    p.addPlotable(data_plot,"PE")

    return p

def make_pull(canvas,xlo,xhi,reg,w,fitres,normset,ch,pads):


    canvas.cd()
    pad     = TPad(reg,reg,xlo,0.125,xhi,0.5)
    pad2    = TPad(reg+"_pull",reg+"_pull",xlo,0.,xhi,0.125)
    pad.SetMargin(0,0,0,0.1)
    pad2.SetMargin(0,0,0.1,0)
    if reg=='sb_hi':
        pad.SetRightMargin(0.1)
        pad2.SetRightMargin(0.1)
    pad.Draw()
    pad2.Draw()

    pad.cd()
    p       = plot(w,fitres,normset,'mj',ch,reg)
    p.GetXaxis().SetTitleSize(0)
    p.GetYaxis().SetRangeUser(0,1450)
    p.Draw()

    pad2.cd()

    pullhist    = p.pullHist("data","WZSM")
    pullhist.SetMaximum(2.9)
    pullhist.SetMinimum(-3)
    pullhist.GetYaxis().SetNdivisions(7)
    pullhist.GetYaxis().SetTitle('#frac{MC-Fit}{error}')
    pullhist.GetYaxis().SetLabelSize(0.125)
    pullhist.GetYaxis().SetTitleSize(0.2)
    pullhist.GetYaxis().SetTitleOffset(0.2)
    pullhist.SetMarkerStyle(20)
    pullhist.SetMarkerSize(1)
    pullhist.SetLineColor(kBlack)
    pullhist.SetMarkerColor(kBlack)
    pullhist.Draw("")
    canvas.Update()

    pads.append(pad)
    pads.append(pad2)


def plot_all(w,ch="el",name='test.png'):
    
    canvas = TCanvas(ch,ch,1250,600)

    offset = 0.15
    pads = []

    make_pull(canvas,0+offset,5/22.+offset/2,'sb_lo',w,fitres,normset,ch,pads)
    make_pull(canvas,5/22.+offset/2,13/22.-offset/2,'sig',w,fitres,normset,ch,pads)
    make_pull(canvas,13/22.-offset/2,1-offset,'sb_hi',w,fitres,normset,ch,pads)

    canvas.cd()
    canvas.Draw()
    canvas.Update()

    regs    = ['sb_lo','sig','sb_hi']

    ratio_style     = TH1D('ratio_style','ratio_style',36,900,4500)
    ratio_style.SetMaximum(2.9)
    ratio_style.SetMinimum(-3)
    ratio_style.GetYaxis().SetNdivisions(7)
    ratio_style.GetYaxis().SetTitle('#frac{MC-Fit}{error}')
    ratio_style.GetYaxis().SetLabelSize(0.125)
    ratio_style.GetYaxis().SetTitleSize(0.2)
    ratio_style.GetYaxis().SetTitleOffset(0.2)
    ratio_style.SetMarkerStyle(20)
    ratio_style.SetMarkerSize(1)

    for i in range(3):
        pad1        = TPad('pad%s'%i,'pad%s'%i,i*0.33,0.645,(i+1)*0.33,1.)
        pad_pull    = TPad('pad_pull%s'%i,'pad_pull%s'%i,i*0.33,0.52,(i+1)*0.33,0.645)
        pads.append(pad1)
        pads.append(pad_pull)
        p=plot(w,fitres,normset,'mlvj',ch,regs[i])
        p.GetYaxis().SetRangeUser(7e-2,5e3)

        canvas.cd()
        pad1.Draw()
        pad_pull.Draw()

        pad1.cd()
        pad1.SetLogy()
        pad1.SetTicky()
        pad1.SetBottomMargin(0)
        if regs[i]!='sg':
         p.Draw()

         pad_pull.cd()   
         pad_pull.SetTopMargin(0)
         pullhist = p.pullHist("data","WJets")
         ratio_style.Draw("")
         pullhist.SetLineColor(kBlack)
         pullhist.SetMarkerStyle(20)
         pullhist.SetMarkerSize(1)
         pullhist.SetMarkerColor(kBlack)
         pullhist.Draw("SAME PE")
	 #for i in range(pullhist.GetN()):
		#print pullhist.GetY()[i]
	 #raw_input("Printed pull hist")

        canvas.Update()
    canvas.SaveAs(name)
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

plot_all(w,options.ch,'prefit_%s.png'%options.ch)

string = '{:>40} : {:>30} / {:>30}\n'.format('>>name<<','>>pre-fit<<','>>post-fit<<')
for i in range(fitparas.getSize()):
    prefit  = str(round(w.var(fitparas.at(i).GetName()).getVal(),6)) + ' +- ' + str(round(w.var(fitparas.at(i).GetName()).getError(),6))
    postfit = str(round(fitparas.at(i).getVal(),6)) + ' +- ' + str(round(fitparas.at(i).getError(),6))
    string += '{:>40} : {:>30} / {:>30}\n'.format(fitparas.at(i).GetName(),prefit,postfit)
for i in range(fitparas.getSize()):
    w.var(fitparas.at(i).GetName()).setVal(fitparas.at(i).getVal())


plot_all(w,options.ch,'postfit_%s.png'%options.ch)


print string
print 'expected events el: ' + str(w.var("normvar_WJets_el").getVal()+w.var("normvar_TTbar_el").getVal()+w.var("normvar_STop_el").getVal()+w.var("normvar_WW_el").getVal()+w.var("normvar_WZ_el").getVal())
print 'observed events el: ' + str(w.data("data_obs").sumEntries("CMS_channel==CMS_channel::ch1||CMS_channel==CMS_channel::ch3||CMS_channel==CMS_channel::ch5"))
print 'expected events mu: ' + str(w.var("normvar_WJets_mu").getVal()+w.var("normvar_TTbar_mu").getVal()+w.var("normvar_STop_mu").getVal()+w.var("normvar_WW_mu").getVal()+w.var("normvar_WZ_mu").getVal())
print 'observed events mu: ' + str(w.data("data_obs").sumEntries("CMS_channel==CMS_channel::ch2||CMS_channel==CMS_channel::ch4||CMS_channel==CMS_channel::ch6"))


