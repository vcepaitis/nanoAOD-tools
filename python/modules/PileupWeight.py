import os
import sys
import math
import json
import ROOT
import random
import numpy

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module

class PileupWeight(Module):
    def __init__(
        self,
        outputName ="puweight",
        processName=None,
        globalOptions={"isData":False, "year":2016}
    ):
        self.outputName = outputName
        self.globalOptions = globalOptions
        self.processName = processName

    def beginJob(self):
        if not self.globalOptions["isData"]:
            if self.globalOptions["year"] == 2016:
                #read in up and down files
                self.dataFile = os.path.expandvars("$CMSSW_BASE/src/PhysicsTools/NanoAODTools/data/pu/2016/data_pileup_2016_69200.root")
                self.dataFile_up = os.path.expandvars("$CMSSW_BASE/src/PhysicsTools/NanoAODTools/data/pu/2016/data_pileup_2016_72500.root")
                self.dataFile_down = os.path.expandvars("$CMSSW_BASE/src/PhysicsTools/NanoAODTools/data/pu/2016/data_pileup_2016_65500.root")
                self.mcFile =  os.path.expandvars("$CMSSW_BASE/src/PhysicsTools/NanoAODTools/data/pu/2016/pileup.root")
            elif self.globalOptions["year"] == 2017:
                self.dataFile = os.path.expandvars("$CMSSW_BASE/src/PhysicsTools/NanoAODTools/data/pu/2017/data_pileup_2017_69200.root")
                self.dataFile_up = os.path.expandvars("$CMSSW_BASE/src/PhysicsTools/NanoAODTools/data/pu/2017/data_pileup_2017_72500.root")
                self.dataFile_down = os.path.expandvars("$CMSSW_BASE/src/PhysicsTools/NanoAODTools/data/pu/2017/data_pileup_2017_65500.root")
                self.mcFile =  os.path.expandvars("$CMSSW_BASE/src/PhysicsTools/NanoAODTools/data/pu/2017/pileup.root")
            elif self.globalOptions["year"] == 2018:
                self.dataFile = os.path.expandvars("$CMSSW_BASE/src/PhysicsTools/NanoAODTools/data/pu/2018/data_pileup_2018_69200.root")
                self.dataFile_up = os.path.expandvars("$CMSSW_BASE/src/PhysicsTools/NanoAODTools/data/pu/2018/data_pileup_2018_72500.root")
                self.dataFile_down = os.path.expandvars("$CMSSW_BASE/src/PhysicsTools/NanoAODTools/data/pu/2018/data_pileup_2018_65500.root")
                self.mcFile =  os.path.expandvars("$CMSSW_BASE/src/PhysicsTools/NanoAODTools/data/pu/2018/pileup.root")
            else:
                print "wrong year selected", year
                sys.exit(1)

            self.mcHistPerProcess = {}
            fMC = ROOT.TFile(self.mcFile)
            if not fMC:
                print "ERROR: Cannot find pileup file: ",self.mcFile
                sys.exit(1)
            for k in fMC.GetListOfKeys():
                self.mcHistPerProcess[k.GetName()] = fMC.Get(k.GetName()).Clone(k.GetName()+str(random.random()))
                self.mcHistPerProcess[k.GetName()].SetDirectory(0)
            fMC.Close()

            #add up and down hists in a loop
            for var in ["_up", "", "_down"]:

                fData = ROOT.TFile(getattr(self, "dataFile"+var))

                if not fData:
                    print "ERROR: Cannot find pileup file: ",getattr(self, "dataFile"+var)
                    sys.exit(1)

                setattr(self, "dataHist"+var, fData.Get("pileup").Clone())
                getattr(self, "dataHist"+var).SetDirectory(0)

    def getWeight(self,nTrueInteractions):

        mcBin = self.mcHist.FindBin(nTrueInteractions)
        w = []
        #add w_up and down
        for var in ["_up", "", "_down"]:
            dataBin = getattr(self, "dataHist"+var).FindBin(nTrueInteractions)
            w.append(getattr(self, "dataHist"+var).GetBinContent(dataBin)/(self.mcHist.GetBinContent(mcBin)+self.mcHist.Integral()*0.0001))
            if w[-1]>5.:
                w[-1] = 0

        #w_up >= w >= w_down
        w = numpy.sort(w)[::-1]

        return w

    def endJob(self):
        pass

    def normHist(self,hist):
        #normalization makes weight independent of binning scheme/range of histograms
        hist.Scale(1./hist.Integral())
        for ibin in range(hist.GetNbinsX()):
            w = hist.GetBinWidth(ibin+1)
            c = hist.GetBinContent(ibin+1)
            hist.SetBinContent(ibin+1,c/w)

    '''
    def interpolateHist(self,hist,binning):
        newHist = ROOT.TH1F("new"+hist.GetName()+str(random.random()),"",
            len(binning)-1,binning
        )
        newHist.SetDirectory(0)
        for ibin in range(1,newHist.GetNbinsX()-1):
            oldBin = hist.GetXaxis().FindBin(newHist.GetBinCenter(ibin+1))
            if newHist.GetBinCenter(ibin+1)>hist.GetBinCenter(oldBin):
                leftC = hist.GetBinContent(oldBin)
                leftP = hist.GetBinCenter(oldBin)

                rightC = hist.GetBinContent(oldBin+1)
                rightP = hist.GetBinCenter(oldBin+1)
            else:
                leftC = hist.GetBinContent(oldBin-1)
                leftP = hist.GetBinCenter(oldBin-1)

                rightC = hist.GetBinContent(oldBin)
                rightP = hist.GetBinCenter(oldBin)

            frac = (newHist.GetBinCenter(ibin+1)-leftP)/(rightP-leftP)
            interC = frac*leftC+(1-frac)*rightC
            newHist.SetBinContent(interC)
    '''

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        if not self.globalOptions["isData"]:
            self.mcHist = None
            for process in self.mcHistPerProcess.keys():
                processName = inputFile.GetName()
                if self.processName!=None:
                    processName = self.processName
                if processName.find(process)>=0:
                    self.mcHist = self.mcHistPerProcess[process]
                    break
            if self.mcHist==None:
                print "ERROR: Cannot find pileup profile for file: "+inputFile.GetName()
                sys.exit(1)

            self.normHist(self.mcHist)

            for var in ["_up", "", "_down"]:
                self.out.branch(self.outputName+var,"F")
                self.normHist(getattr(self, "dataHist"+var))

            self.sum2 = 0
            self.sum = 0
            self.n = 0

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        if not self.globalOptions["isData"] and self.n>0 and (self.sum2/(1.*self.n))>(self.sum**2/(1.*self.n**2)):
            avg = 1.*self.sum/self.n
            sig = math.sqrt(self.sum2/(1.*self.n)-self.sum**2/(1.*self.n**2))

            print "Average pileup weight (%s): %6.3f +- %6.3f"%(self.outputName,avg,sig)

    def analyze(self, event):
        if not self.globalOptions["isData"]:
            puWeight = numpy.ones(3)
            puWeight = self.getWeight(event.Pileup_nTrueInt)
            self.n += 1
            self.sum+=puWeight[1]
            self.sum2+=puWeight[1]**2

            for i, var in enumerate(["_up", "", "_down"]):
                self.out.fillBranch(self.outputName+var,puWeight[i])

        return True
