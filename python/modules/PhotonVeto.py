import os
import sys
import math
import json
import ROOT
import random
from utils import deltaR

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module

class PhotonVeto(Module):

    def __init__(
        self,
        inputCollection = lambda event: Collection(event, "Photon"),
        outputName = "vetoPhotons",
        photonMinPt = 10.,
        photonMaxEta = 2.5,
        globalOptions={"isData":False}
    ):
        self.globalOptions = globalOptions
        self.inputCollection = inputCollection
        self.outputName = outputName
        self.photonMinPt = photonMinPt
        self.photonMaxEta = photonMaxEta
 
    def beginJob(self):
        pass
        
    def endJob(self):
        pass
        
    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        self.out.branch("n"+self.outputName,"I")
        
    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass
        
    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""
        photons = self.inputCollection(event)
        muons = Collection(event, "Muon")
        electrons = Collection(event, "Electron")
 
        selectedphotons = []
        unselectedphotons = []
        
        for photon in photons:
            if photon.pt>self.photonMinPt and math.fabs(photon.eta)<self.photonMaxEta and photon.cutBased>0:
                dr_mu = [deltaR(photon, muon) for muon in muons]
                dr_e = [deltaR(photon, electron) for electron in electrons]
                if len(dr_mu) > 0:
                    min_drmu = min(dr_mu)
                    if min_drmu < 0.4:
                        unselectedphotons.append(photon)
                        continue
                if len(dr_e) > 0:
                    min_dre = min(dr_e)
                    if min_dre < 0.4:
                        unselectedphotons.append(photon)
                        continue
                selectedphotons.append(photon)
 
            else:
                unselectedphotons.append(photon)
  
        self.out.fillBranch("n"+self.outputName,len(selectedphotons))
        
        setattr(event,self.outputName,selectedphotons)
        setattr(event,self.outputName+"_unselected",unselectedphotons)

        return True
        
