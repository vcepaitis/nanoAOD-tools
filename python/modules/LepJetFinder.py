import os
import sys
import math
import json
import ROOT
import random


import utils 
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module


class LepJetFinder(Module):
    def __init__(
        self,
        jetCollection,
        leptonCollection,
        outputName="lepJet",
    ):
        self.jetCollection = jetCollection
        self.leptonCollection = leptonCollection
        self.outputName = outputName

    def beginJob(self):
        pass

    def endJob(self):
        pass

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        self.out.branch(self.outputName+"_lepJets_index",
                        "F", lenVar="n"+self.outputName)
        self.out.branch(self.outputName+"_pt",
                        "F", lenVar="n"+self.outputName)
        self.out.branch(self.outputName+"_eta",
                        "F", lenVar="n"+self.outputName)
        self.out.branch(self.outputName+"_phi",
                        "F", lenVar="n"+self.outputName)
        self.out.branch(self.outputName+"_deltaR",
                        "F", lenVar="n"+self.outputName)
        self.out.branch(self.outputName+"_deltaR_jet",
                        "F", lenVar="n"+self.outputName)
        self.out.branch(self.outputName+"_deltaPhi_jet",
                        "F", lenVar="n"+self.outputName)

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def analyze(self, event):
        jetCollection = self.jetCollection(event)
        leptonCollection = self.leptonCollection(event)

        jet_pts = []
        jet_etas = []
        jet_phis = []
        jet_deltaRs = []
        jet_deltaRsjets = []
        jet_deltaPhisjets = []
        lepJets = []
	lepJets_index = []

        for lepton in leptonCollection:
            jet = jetCollection[0]
            deltaR = lepton.p4().DeltaR(jet.p4())
	    jet_index = 0
            for i , _jet in enumerate(jetCollection):
                _deltaR = lepton.p4().DeltaR(_jet.p4())
                if _deltaR < deltaR:
                    jet = _jet
                    deltaR = _deltaR
		    jet_index = i 

            lepJets.append(jet)
	    lepJets_index.append(jet_index)
            jet_deltaRs.append(lepton.p4().DeltaR(jet.p4()))
            jet_pts.append(jet.pt)
            jet_etas.append(jet.eta)
            jet_phis.append(jet.phi)

	    
	    lepjet = lepJets[0]
	    lepjet_ind = lepJets_index[0]	
	    deltaR_jet = 10.
	    deltaPhi_jet = 10.
       	    for i2 , jet2 in enumerate(jetCollection): 
	     if i2 != lepjet_ind:  
		_deltaR = lepjet.p4().DeltaR(jet2.p4()) 
		_deltaPhi_jet = utils.deltaPhi(lepjet.phi , jet2.phi)
	   	if  _deltaR < deltaR_jet: 
			deltaR_jet  = _deltaR
		if _deltaPhi_jet < deltaPhi_jet: 
			deltaPhi_jet = _deltaPhi_jet

            jet_deltaRsjets.append(deltaR_jet)
	    jet_deltaPhisjets.append(deltaPhi_jet) 
		 	
        self.out.fillBranch(self.outputName+"_lepJets_index", lepJets_index)
        self.out.fillBranch(self.outputName+"_pt", jet_pts)
        self.out.fillBranch(self.outputName+"_eta", jet_etas)
        self.out.fillBranch(self.outputName+"_phi", jet_phis)
        self.out.fillBranch(self.outputName+"_deltaR", jet_deltaRs)
        self.out.fillBranch(self.outputName+"_deltaR_jet", jet_deltaRsjets)
        self.out.fillBranch(self.outputName+"_deltaPhi_jet", jet_deltaPhisjets)
        setattr(event, self.outputName, lepJets)

        return True
