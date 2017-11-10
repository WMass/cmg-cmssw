from CMGTools.TTHAnalysis.treeReAnalyzer import *
from CMGTools.WMass.tools.standaloneElectronCalibrator import ElectronCalibrator
from CMGTools.WMass.tools.PileUpReWeighter import PileUpReWeighter
from CMGTools.WMass.tools.eventVars_wmass import SimpleVBoson
import types, os, math
import numpy as np

class DrellYanAngularReweighter:
    def __init__(self,hname,verbose=False):
        self.verbose = verbose
        plotfile = "%s/src/CMGTools/WMass/data/theory/angularcoeffs/cms_data_8TeV_distributions.root" % os.environ['CMSSW_BASE']
        f_plots = ROOT.TFile.Open(plotfile)
        histo_data = f_plots.Get(hname+"_data")
        histo_bkg = f_plots.Get(hname+"_background")
        self.histo_sig = (f_plots.Get(hname+"_signal")).Clone("weights")
        self.histo_sig.SetDirectory(None)
        histo_num = histo_data-histo_bkg
        self.myvals = self.load(self.histo_sig)
        self.targetvals = self.load(histo_num)
        def w2(t,m):
            if t == 0: return (0 if m else 1)
            return (t/m if m else 1)
        self.weights = [ w2(t,m) for (m,t) in zip(self.myvals,self.targetvals) ]
        if self.verbose:
            print "Raw weights for variable %s, with a number of bins equal to %d; max %.3f, min %.3f, avg %.3f" % (
                hname, len(self.weights), max(self.weights), min(self.weights), sum(self.weights)/len(self.weights) )
            print self.weights
        self.fixLargeWeights()
        if self.verbose:
            print "Initalized weights for variable %s, with a number of bins equal to %d; max %.3f, min %.3f, avg %.3f" % (
                hname, len(self.weights), max(self.weights), min(self.weights), sum(self.weights)/len(self.weights) )
            print self.weights
    def fixLargeWeights(self,maxshift=0.0025,hardmax=1.5):
        def checkIntegral(weights):
            myint  = sum(a*b for (a,b) in zip(weights,     self.myvals))
            refint = sum(a*b for (a,b) in zip(self.weights,self.myvals))
            return (myint-refint)/refint
        maxw = min(max(self.weights),5)
        while maxw > hardmax:
            cropped = [ min(maxw,w) for w in self.weights ]
            check = checkIntegral(cropped)
            if self.verbose:
                print "For maximum weight %.3f: integral match: %.5f" % (maxw, check)
            if abs(check) > maxshift:
                break
            maxw *= 0.9
        maxw /= 0.9
        cropped = [ min(maxw,w) for w in self.weights ]
        normshift = checkIntegral(cropped)
        recalibrated = [ c*(1-normshift) for c in cropped ]
        if self.verbose:
            print "Cropped weights up to maximum %d. Normalization shift %.5f, corrected overall to %g" % (maxw,normshift,checkIntegral(recalibrated))
        self.weights = recalibrated
    def load(self,hist,norm=True):
        vals = [ max(hist.GetBinContent(i),0.) for i in xrange(0,hist.GetSize()) ]
        if self.verbose:
            print "Normalization of ",hist.GetName(),": ",sum(vals)
        if norm: 
            scale = 1.0/sum(vals)
            vals = [ v*scale for v in vals ]
        return vals
    def __call__(self,valx,valy=None):
        xmin = self.histo_sig.GetXaxis().GetXmin()+1e-3
        xmax = self.histo_sig.GetXaxis().GetXmax()-1e-3
        bin = self.histo_sig.FindBin(min(max(valx,xmin),xmax))
        if valy!=None:
            if not self.histo_sig.InheritsFrom("TH2"):
                print "ERROR! The weights input istogram for signal: ",histo_sig.GetName()," is not a 2D histogram! Returning weight 0!"
                return 0.
            ymin = self.histo_sig.GetYaxis().GetXmin()+1e-3
            ymax = self.histo_sig.GetYaxis().GetXmax()-1e-3
            bin = self.histo_sig.FindBin(min(max(valx,xmin),xmax),min(max(valy,ymin),ymax))
        return self.weights[bin]

class ZtoWTransport:
    def __init__(self,hname,plotfile=None):
        if plotfile == None: plotfile="%s/src/CMGTools/WMass/data/theory/angularcoeffs/woverz.root" % os.environ['CMSSW_BASE']
        tf = ROOT.TFile.Open(plotfile)
        self.ratio = tf.Get(hname).Clone("histo")
    def __call__(self,valx):
        xmin = self.ratio.GetXaxis().GetXmin()+1e-3
        xmax = self.ratio.GetXaxis().GetXmax()-1e-3
        bin = self.ratio.FindBin(min(max(valx,xmin),xmax))
        return self.ratio.GetBinContent(bin)

class KinematicVars:
    def __init__(self):
        ROOT.gSystem.Load("libFWCoreFWLite.so")
        ROOT.AutoLibraryLoader.enable()
        ROOT.gSystem.Load("libCMGToolsWMass.so")
        self.zpt_weight = DrellYanAngularReweighter("scaledptZ",True)
        self.aipi_weight = DrellYanAngularReweighter("y_vs_sumAiPi",True)
        self.cth_weight = DrellYanAngularReweighter("y_vs_ctheta")
        self.phi_weight = DrellYanAngularReweighter("y_vs_phi")
        self.ZtoW = ZtoWTransport("gen_scaledptv_WZ_ratio")
    def initSample(self,nevents,dataset):
        self.isMC = not any(x in dataset for x in "DoubleMu DoubleEl DoubleEG MuEG MuonEG SingleMu SingleEl".split())
        aifilename = "%s/src/CMGTools/WMass/data/theory/angularcoeffs/Ai_8TeV_atlas.root" % os.environ['CMSSW_BASE']
        tfile = ROOT.TFile.Open(aifilename)
        etabins = ["00","10","20","35"]
        self.Ai = {}
        for i in range(5):
            for e in range(len(etabins)-1):
                aistr = "A%d_8TeV_atlas_eta%s_%s" % (i,etabins[e],etabins[e+1])
                self.Ai[(i,e)] = tfile.Get(aistr)
    def listBranches(self):
        mylist = [ ("nAP","I"), ("AP","F",8,"nAP"),
                   ("z_y","F"), ("costheta_cs","F"), ("phi_cs","F"),
                   ("zpt_w","F"), ("zpt_w_up","F"), ("zpt_w_dn","F"),
                   ("aipi_w","F"), ("cth_w","F"), ("phi_w","F")]
        self.branches = mylist
        print "self.branches = ",self.branches[:]
        return self.branches[:]
    def CSFrame(self,dilepton):
        pMass = 0.938272
        beamE = 4000;
        sign = np.sign(dilepton.Z())
        proton1 = ROOT.TLorentzVector(0.,0.,sign*beamE,hypot(beamE,pMass));  proton2 = ROOT.TLorentzVector(0.,0.,-sign*beamE,hypot(beamE,pMass))
        proton1.Boost(-dilepton.BoostVector()); proton2.Boost(-dilepton.BoostVector())
        CSAxis = (proton1.Vect().Unit()-proton2.Vect().Unit()).Unit()
        yAxis = (proton1.Vect().Unit()).Cross((proton2.Vect().Unit()));
        yAxis = yAxis.Unit();
        xAxis = yAxis.Cross(CSAxis);
        xAxis = xAxis.Unit();
        return (xAxis,yAxis,CSAxis)
    def cosThetaCS(self,lplus,lminus):
        dilep = lplus + lminus
        boostedLep = ROOT.TLorentzVector(lminus)
        boostedLep.Boost(-dilep.BoostVector())
        csframe = self.CSFrame(dilep)
        return cos(boostedLep.Angle(csframe[2]))
    def phiCS(self,lplus,lminus):
        dilep = lplus + lminus
        boostedLep = ROOT.TLorentzVector(lminus)
        boostedLep.Boost(-dilep.BoostVector())
        csframe = self.CSFrame(dilep)
        phi = atan2((boostedLep.Vect()*csframe[1]),(boostedLep.Vect()*csframe[0]))
        if(phi<0): return phi + 2*ROOT.TMath.Pi()
        else: return phi
    def P_i(self,theta,phi,i):
        ct = cos(theta); st = sqrt(1-pow(ct,2))
        cp = cos(phi); sp = sqrt(1-pow(cp,2))
        if i==0: return 0.5*(1-3*ct*ct)
        elif i==1: return 2*ct*st*cp
        elif i==2: return 0.5*st*st*(cp*cp-sp*sp)
        elif i==3: return st*cp
        elif i==4: return ct
        elif i==5: return 2*st*st*sp*cp
        elif i==6: return 2*st*ct*sp
        elif i==7: return st*sp
        else: return -999.
    def A_i(self,i,ptz,yz):
        if abs(yz)<1: histo = self.Ai[(i,0)]
        elif abs(yz)<2: histo = self.Ai[(i,1)]
        else: histo = self.Ai[(i,2)]
        if histo.GetEntries() == 0: return 0
        xmin = histo.GetXaxis().GetXmin()+1e-3
        xmax = histo.GetXaxis().GetXmax()-1e-3
        xbin = histo.GetXaxis().FindBin(min(max(ptz,xmin),xmax))
        return histo.GetBinContent(xbin)
    def PtEtaPhiM4V(self,pt,eta,phi,pdgId,motherId=None):
        tlv = ROOT.TLorentzVector()
        m = 0.000511 if abs(pdgId)==11 else 0.1057
        tlv.SetPtEtaPhiM(pt,eta,phi,m)
        tlv.charge = np.sign(pdgId)
        if motherId!=None: tlv.motherId = motherId
        return tlv
    def __call__(self,event):
        # prepare output
        ret = {}; 
        genp = [p for p in Collection(event,"GenP6StatusThree","nGenP6StatusThree")] if self.isMC else []
        leps = [l for l in Collection(event,"LepGood","nLepGood")]

        leps_4v=[ self.PtEtaPhiM4V(l.pt,l.eta,l.phi,l.pdgId) for l in leps ]
        genleps_4v=[ self.PtEtaPhiM4V(p.pt,p.eta,p.phi,p.pdgId,p.motherId) for p in genp if ( (abs(p.pdgId)>=11 and abs(p.pdgId)<=14) and (p.motherId==23 or abs(p.motherId)==24) and abs(p.grandmaId)!=6 ) ]
        genleps_4v.sort(key = lambda p: p.Pt(), reverse = True)

        ret["nAP"] = 5 if len(leps)>=2 else 0
        ret["AP"]=[]
        if len(leps)>=2:
            (lplus,lminus) = (leps_4v[0],leps_4v[1]) if leps[0].charge>0 else (leps_4v[1],leps_4v[0])
            Z = SimpleVBoson(leps_4v[:2])
            costheta_cs = self.cosThetaCS(lplus,lminus); phi_cs = self.phiCS(lplus,lminus);
            sumAiPi = 1+costheta_cs*costheta_cs+sum([self.A_i(i,Z.pt(),Z.y())*self.P_i(acos(costheta_cs),phi_cs,i) for i in range(ret["nAP"])])            
            ret["costheta_cs"] = costheta_cs; ret["phi_cs"] = phi_cs
            ret["z_y"] = Z.y()
            for i in range(ret["nAP"]):
                ret["AP"].append(self.A_i(i,Z.pt(),Z.y())*self.P_i(acos(costheta_cs),phi_cs,i))
        else:
            ret["costheta_cs"] = ret["phi_cs"] = ret["z_y"] = -999

        # kinematic reweighting
        if len(genleps_4v)==2:
            (lplus,lminus) = (genleps_4v[0],genleps_4v[1]) if genleps_4v[1].charge<0 else (genleps_4v[1],genleps_4v[0]) # Z or W- case: use l-
            if genleps_4v[0].motherId==24: (lplus,lminus) = (genleps_4v[0],genleps_4v[1]) if genleps_4v[1].charge==0 else (genleps_4v[1],genleps_4v[0]) # W+ case: use l+
            Z = SimpleVBoson(genleps_4v[:2])
            costheta_cs = self.cosThetaCS(lplus,lminus); phi_cs = self.phiCS(lplus,lminus);
            sumAiPi = 1+costheta_cs*costheta_cs+sum([self.A_i(i,Z.pt(),Z.y())*self.P_i(acos(costheta_cs),phi_cs,i) for i in range(ret["nAP"])])
            ret["aipi_w"] = self.aipi_weight(sumAiPi,abs(Z.y()))
            ret["cth_w"] = self.cth_weight(costheta_cs,abs(Z.y()))
            ret["phi_w"] = self.phi_weight(phi_cs,abs(Z.y()))
            ptOverM = Z.pt()/Z.mll()*91.1876
            ret["zpt_w"] = self.zpt_weight(ptOverM)
            ret["zpt_w_up"] = self.ZtoW(ptOverM)
            ret["zpt_w_dn"] = self.ZtoW(ptOverM)
        else:
            ret["aipi_w"] = ret["cth_w"] = ret["phi_w"] = ret["zpt_w"] = 1.0

        ### return
        return ret

if __name__ == '__main__':
    from sys import argv
    file = ROOT.TFile(argv[1])
    tree = file.Get("tree")
    tree.vectorTree = True
    class Tester(Module):
        def __init__(self, name):
            Module.__init__(self,name,None)
            self.sf = EventVarsWmass()
        def analyze(self,ev):
            #if ev.metNoMu_pt < 200: return True
            print "\nrun %6d lumi %4d event %d: metNoMu %d" % (ev.run, ev.lumi, ev.evt, ev.metNoMu_pt)
            print self.sf(ev)
    el = EventLoop([ Tester("tester") ])
    el.loop([tree], maxEvents = 50)

