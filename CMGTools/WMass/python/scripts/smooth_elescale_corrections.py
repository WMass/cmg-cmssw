#!/usr/bin/env python
from math import *
import re
import os, os.path
from array import array
from CMGTools.WMass.plotter.mcPlots import *

## safe batch mode
import sys
args = sys.argv[:]
sys.argv = ['-b']
import ROOT
sys.argv = args
ROOT.gROOT.SetBatch(True)
ROOT.PyConfig.IgnoreCommandLineOptions = True

XMIN=20
XMAX=80

def parseCorrectionsFile(filename):
    corrections = {}
    inputfile = open(filename)
    lines = inputfile.readlines()
    for line in lines:
        l=line.rstrip()
        etarange = "_".join((l.split("-")[0]).split("_")[-2:])
        print l
        elecat = l.split("-")[1]
        etrange = ((l.split(' ')[0].rstrip()).split("-")[2]).split("_")[-2:]
        correction = l.split()[4:6]
        #print etarange, "   ",elecat,"   ",etrange,"   ",correction
        category = (etarange,elecat)
        if category not in corrections: corrections[category] = []
        corrections[category].append(etrange+correction)
    return corrections

def doFit(corrections,category,xmin,xmax,closure=False):
    ROOT.gROOT.ProcessLine(".x %s/src/CMGTools/WMass/python/plotter/tdrstyle.cc" %  os.environ['CMSSW_BASE'])
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetOptTitle(0)
    corr_bycat = corrections[category]
    edges = [float(c[0]) for c in corr_bycat]
    edges.append(100.)
    print edges
    graph = ROOT.TGraphErrors(len(edges)-2)
    for b,corr in enumerate(corr_bycat):
        if b==0: continue # skip the 10 GeV bin
        graph.SetPoint(b-1,(edges[b+1]+edges[b])/2.,float(corr[2]))
        graph.SetPointError(b-1,(edges[b+1]-edges[b])/2.,float(corr[3]))
    graph.GetXaxis().SetTitle("electron p_{T} [GeV]")
    graph.GetYaxis().SetTitle("correction")
    
    outputName="eta_"+category[0]+"_"+category[1]
    if closure: outputName = outputName+"_closure"
    c1 = ROOT.TCanvas(outputName+"_canvas", outputName,600,600)
    graph.Draw('ape')

    pol = ROOT.TF1("f1","[0]/x+[1]/sqrt(x)+[2]",xmin,xmax)
    pol.SetParLimits(0,0,10)
    graph.Fit('f1','R')
    pdir = "plots/corrections/v0/step5/"
    for ext in [".pdf",".png"]:
        c1.Print(pdir+outputName+ext)
    return pol

def dumpPolynomial(pol,xstep,etab,elecat,outfile):
    catPfx = 'absEta_'+etab+'-'+elecat+'-Et_'
    catlen = max(30,len(catPfx))
    catpatt = "%%-%ds" % catlen
    catmin = catPfx+('0_%d' % XMIN)
    outfile.write((catpatt % catmin) +'   runNumber   3   999999   {CORR:.4f} 0.0   0.0   0.0\n'.format(CORR=pol.Eval(XMIN)))
    for i in range(XMIN,XMAX,xstep):
        icat = catPfx+('%d_%d' % (i,i+xstep))
        outfile.write((catpatt % icat) +'   runNumber   3   999999   {CORR:.4f} 0.0   0.0   0.0\n'.format(CORR=pol.Eval(i+xstep/2.)))
    catmax = catPfx+('%d_10000' % XMAX)
    outfile.write((catpatt % catmax) +'   runNumber   3   999999   {CORR:.4f} 0.0   0.0   0.0\n'.format(CORR=pol.Eval(XMAX)))
        

if __name__ == "__main__":
    corr_file = ("%s/src/EgammaAnalysis/ElectronTools/data/WMass_Winter17_ResidualCorrections_ele_scales.dat" % os.environ['CMSSW_BASE']);
    corrs = parseCorrectionsFile(corr_file)
    etabins=['0_1','1_1.4442','1.566_2','2_2.5']
    elecats=['gold','bad']

    corr_file_fine = ("%s/src/EgammaAnalysis/ElectronTools/data/WMass_Winter17_ResidualCorrections_ele_scales_fine.dat" % os.environ['CMSSW_BASE']);
    outfile = open(corr_file_fine,'w')

    for etab in etabins:
        for cat in elecats:
            pol = doFit(corrs,(etab,cat),XMIN,XMAX)
            dumpPolynomial(pol,1,etab,cat,outfile)

    # corrs_fine = parseCorrectionsFile(corr_file_fine)
    # for etab in etabins:
    #     for cat in elecats:
    #         doFit(corrs_fine,(etab,cat),XMIN,XMAX,closure=True)
