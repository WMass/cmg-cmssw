#!/usr/bin/env python
import sys, os
from math import *
from array import array
from CMGTools.WMass.plotter.mcPlots import doSpam

## safe batch mode
import sys
args = sys.argv[:]
sys.argv = ['-b']
import ROOT
ROOT.gROOT.SetBatch(True)
ROOT.PyConfig.IgnoreCommandLineOptions = True

def doTinyCmsPrelimStandalone(textLeft="_default_",textRight="_default_",hasExpo=False,textSize=0.033,lumi=None, xoffs=0, options=None):
    if textLeft  == "_default_": textLeft  = "#bf{CMS} #it{Preliminary}"
    if textRight == "_default_": textRight = "%(lumi) (8 TeV)"
    if lumi      == None       : lumi      = 19.7
    if   lumi > 3.54e+1: lumitext = "%.0f fb^{-1}" % lumi
    elif lumi > 3.54e+0: lumitext = "%.1f fb^{-1}" % lumi
    elif lumi > 3.54e-1: lumitext = "%.2f fb^{-1}" % lumi
    elif lumi > 3.54e-2: lumitext = "%.0f pb^{-1}" % (lumi*1000)
    elif lumi > 3.54e-3: lumitext = "%.1f pb^{-1}" % (lumi*1000)
    else               : lumitext = "%.2f pb^{-1}" % (lumi*1000)
    lumitext = "%.1f fb^{-1}" % lumi
    textLeft = textLeft.replace("%(lumi)",lumitext)
    textRight = textRight.replace("%(lumi)",lumitext)
    if textLeft not in ['', None]:
        doSpam(textLeft, (.28 if hasExpo else .17)+xoffs, .955, .60+xoffs, .995, align=12, textSize=textSize)
    if textRight not in ['', None]:
        doSpam(textRight,.68+xoffs, .955, .99+xoffs, .995, align=32, textSize=textSize)

def parseTable(filename):
    results = {}
    inputfile = open(filename)
    lines = inputfile.readlines()
    for i,line in enumerate(lines):
        l=line.rstrip()
        fields=l.split()
        ptrange = fields[1:3]
        etarange = "_".join(fields[4:6])
        datafit = [fields[6],fields[8]]; mcfit = [fields[9],fields[11]]; diff = [fields[12],fields[14]]
        if etarange not in results: results[etarange] = []
        results[etarange].append((ptrange,datafit,mcfit,diff))
    return results

if __name__ == '__main__':

    if len(args) == 0:
        print "You must provide at least one file to parse"
        exit(0)

    ROOT.gROOT.ProcessLine(".x %s/src/CMGTools/WMass/python/plotter/tdrstyle.cc" %  os.environ['CMSSW_BASE'])
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetOptTitle(0)

    pdgMZ0 = 91.1876
    results = parseTable(args[1])
    histos = {}
    for key,etabin in results.iteritems():
        print "---- eta bin: ",key," ----"
        ptedges = []; diffs = []; diffs_e = []
        for b,(ptrange,datafit,mcfit,diff) in enumerate(etabin):
            print "%s - %s GeV   &   %s $\pm$ %s  &  %s $\pm$ %s   &  %.2f $\pm$ %.2f" % (ptrange[0],ptrange[1],datafit[0],datafit[1],mcfit[0],mcfit[1],100*float(diff[0])/pdgMZ0/sqrt(2),100*float(diff[1])/pdgMZ0/sqrt(2))
            for pt in ptrange:
                if float(pt) not in ptedges: ptedges.append(float(pt))
            diffs.append(float(diff[0])/pdgMZ0/sqrt(2))
            diffs_e.append(float(diff[1])/pdgMZ0/sqrt(2))
        histos[key] = ROOT.TH1F("scale_%s" % key,"",len(ptedges)-1,array('f',ptedges))
        print "bins = ",ptedges
        print "vals = ",diffs
        for b in xrange(len(diffs)):
            histos[key].SetBinContent(b+1,diffs[b])
            histos[key].SetBinError(b+1,diffs_e[b])
        histos[key].GetXaxis().SetTitle("p_{T} [GeV]")
        histos[key].GetYaxis().SetTitle("p (data-sim)/true")

    c1 = ROOT.TCanvas("canvas","",600,600)
    leg = ROOT.TLegend(.70, .15, .90, .35)
    leg.SetFillColor(0)
    leg.SetShadowColor(0)
    leg.SetLineColor(0)
    leg.SetTextFont(42)
    leg.SetTextSize(0.035)
    ids = ["0.000_1.479","1.479_2.500"]
    labels = ["barrel","endcap"]
    mkColors = [ROOT.kRed+2,ROOT.kGreen+2]
    mkStyles = [ROOT.kOpenCircle,ROOT.kOpenSquare]
    for key,label,color,style in zip(ids,labels,mkColors,mkStyles):
        histos[key].SetMarkerStyle(style)
        histos[key].SetMarkerColor(color); histos[key].SetLineColor(color)
        histos[key].GetYaxis().SetRangeUser(-0.003,0.003)
        histos[key].Draw("same pe" if "endcap" in label else "pe")
        leg.AddEntry(histos[key],label,"pe")

    lineUp = ROOT.TLine(histos["0.000_1.479"].GetXaxis().GetXmin(),1e-3,histos["0.000_1.479"].GetXaxis().GetXmax(),1e-3)
    lineDn = ROOT.TLine(histos["0.000_1.479"].GetXaxis().GetXmin(),-1e-3,histos["0.000_1.479"].GetXaxis().GetXmax(),-1e-3)
    lineUp.SetLineWidth(2); lineUp.SetLineColor(58);
    lineDn.SetLineWidth(2); lineDn.SetLineColor(58);
    lineUp.Draw("L"); lineDn.Draw("L");

    doTinyCmsPrelimStandalone(lumi=19.7)
    leg.Draw()
    
    [c1.SaveAs("electronscale_datamc_pt.%s" % ext) for ext in ["pdf","png"]]
