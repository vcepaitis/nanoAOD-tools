import os
import sys
import math
import json
import ROOT
import random

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module

class XsecWeight(Module):
    def __init__(self,yields,lumi,globalOptions={"isData":False}, outputName=None):
        self.globalOptions=globalOptions
        self.xsecs = {
            "TTToHadronic": 377.96,
            "TTToSemiLeptonic": 365.34,
            "TTTo2L2Nu": 88.29,
            "ST_tW_top_5f": 30.11,
            "ST_tW_antitop": 30.09,
            "ST_t-channel_top_4f": 136.02,
            "ST_t-channel_antitop_4f": 80.95,
            "QCD_Pt-15to20_MuEnriched": 3819570.,
            "QCD_Pt-20to30_MuEnriched": 2960198.,
            "QCD_Pt-30to50_MuEnriched": 1652471.,
            "QCD_Pt-50to80_MuEnriched": 438000.,
            "QCD_Pt-80to120_MuEnriched": 106000.,
            "QCD_Pt-120to170_MuEnriched": 25200.,
            "QCD_Pt-170to300_MuEnriched": 8654.,
            "QCD_Pt-300to470_MuEnriched": 797.,
            "QCD_Pt-470to600_MuEnriched": 78.91,
            "QCD_Pt-600to800_MuEnriched": 25.095,
            "QCD_Pt-800to1000_MuEnriched": 4.707,
            "QCD_Pt-1000toInf_MuEnriched": 1.621,
            "QCD_Pt-20to30_EMEnriched": 5352959.,
            "QCD_Pt-30to50_EMEnriched": 9928000.,
            "QCD_Pt-50to80_EMEnriched": 2890800.,
            "QCD_Pt-80to120_EMEnriched": 350000.,
            "QCD_Pt-120to170_EMEnriched": 62964.,
            "QCD_Pt-170to300_EMEnriched": 18810.,
            "QCD_Pt-300toInf_EMEnriched": 1350.,
            "WToLNu_0J": 50131.98259,
            "WToLNu_1J": 8426.094336,
            "WToLNu_2J": 3172.958208,
            "WJetsToLNu_0J": 50131.98259,
            "WJetsToLNu_1J": 8426.094336,
            "WJetsToLNu_2J": 3172.958208,
            "DYJetsToLL_M-50": 5941.,
            "DYJetsToLL_M-10to50": 18810.,
            "TGJets": 2.967,
            "TTGJets": 3.697,
            "WGToLNuG": 405.271,
            "WZTo3LNu": 4.4297,
            "ZGToLLG": 131.1
        }
        self.yields = yields
        self.lumi = lumi
        
    def beginJob(self):
        pass
        
    def endJob(self):
        pass
        
    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        self.out.branch("xsecweight", "F")
        self.out.branch("isData", "F")
        self.out.branch("processId", "F")
        
        self.processXsec = None
        self.processYield = None
        self.proccessId = 0
        
        for i,processName in enumerate(sorted(self.xsecs.keys())):
            if processName in inputFile.GetName():
                if self.processXsec != None:
                    print "ERROR - ambigious process name '%s' vs '%s' for file '%s' in xsec dict found!"%(
                        processName,self.processXsec,inputFile.GetName()
                    )
                    sys.exit(1)
                self.processXsec = self.xsecs[processName]
                self.proccessId = i+1
        
        for processName in self.yields.keys():
            if processName in inputFile.GetName():
                if self.processYield != None:
                    print "ERROR - ambigious process name '%s' vs '%s' for file '%s' in yield dict found!"%(
                        processName,self.processYield,inputFile.GetName()
                    )
                    sys.exit(1)
                self.processYield = self.yields[processName]['weighted']
        print "Xsec: ",self.processXsec
        print "Yield: ",self.processYield
        print "Lumi: ",self.lumi
            
    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass
        
        
    def analyze(self, event):
        if self.globalOptions['isData']:
            self.out.fillBranch("xsecweight", 1.)
            self.out.fillBranch("isData", 1.)
            self.out.fillBranch("processId", 0.)
        else:
            self.out.fillBranch("xsecweight", event.Generator_weight*self.lumi*self.processXsec/self.processYield)
            self.out.fillBranch("isData", 0.)
            self.out.fillBranch("processId", self.proccessId)
            
        return True
            
            
