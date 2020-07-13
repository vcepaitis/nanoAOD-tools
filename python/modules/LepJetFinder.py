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
        storeKinematics=['pt', 'eta', 'phi', 'jetId', 'deltaR', 'nConstituents', 'jetIdx' , 'deltaR_jet', 'deltaPhi_jet'],
    ):
        self.jetCollection = jetCollection
        self.leptonCollection = leptonCollection
        self.outputName = outputName
        self.storeKinematics = storeKinematics

    def beginJob(self):
        pass

    def endJob(self):
        pass

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree

        for variable in self.storeKinematics:
            self.out.branch(self.outputName+"_"+variable, "F", lenVar="n"+self.outputName)

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def analyze(self, event):
        jetCollection = self.jetCollection(event)
        leptonCollection = self.leptonCollection(event)

        lepJets = []

        if len(jetCollection) == 0:
            setattr(event, self.outputName, lepJets)
            return True

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
	    
	    lepjet = lepJets[0]
	    deltaR_jet = 10.
	    deltaPhi_jet = 10.
       	    for i2 , jet2 in enumerate(jetCollection): 
	     if i2 != jet_index:  
		_deltaR = lepjet.p4().DeltaR(jet2.p4()) 
		_deltaPhi_jet = utils.deltaPhi(lepjet.phi , jet2.phi)
	   	if  _deltaR < deltaR_jet: 
			deltaR_jet  = _deltaR
		if _deltaPhi_jet < deltaPhi_jet: 
			deltaPhi_jet = _deltaPhi_jet

            setattr(jet, "deltaR", deltaR)
            setattr(jet, "jetIdx", jet_index)
            setattr(jet , "deltaR_jet" , deltaR_jet)
            setattr(jet , "deltaPhi_jet", deltaPhi_jet)

        for variable in self.storeKinematics:
            self.out.fillBranch(self.outputName+"_"+variable,
                                map(lambda jet: getattr(jet, variable), lepJets))
        setattr(event, self.outputName, lepJets)

        return True
