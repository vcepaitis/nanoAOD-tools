import os
import sys
import math
import json
import ROOT
import random

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module

from utils import deltaR, getCtauLabel


class JetTaggerResult(Module):

    def __init__(
        self,
        inputCollection="selectedJets_nominal",
        taggerName="llpdnnx",
        profiledLabels = ['LLP_Q','LLP_QE','LLP_QMU'],
        globalOptions={"isData": False}
    ):
        self.inputCollection = inputCollection
        self.taggerName = taggerName
        self.profiledLabels = profiledLabels
        self.globalOptions = globalOptions

    def beginJob(self):
        pass

    def endJob(self):
        pass

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        self.out.branch(self.inputCollection+"_"+self.taggerName+"_parameter","F",lenVar="n"+self.inputCollection)
        for label in self.profiledLabels:
            self.out.branch(self.inputCollection+"_"+self.taggerName+"_"+label,"F",lenVar="n"+self.inputCollection)

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""

        jets = getattr(event, self.inputCollection)

        jetPredictionsPerLabel = {}
        jetParameters = []
        for label in self.profiledLabels:
            jetPredictionsPerLabel[label] = []

        for ijet, jet in enumerate(jets):
            for label in self.profiledLabels:
                if not hasattr(jet, self.taggerName):
                    print "WARNING - jet ", jet, " has no ", self.taggerName, " set to minus 1"
                    jetPredictionsPerLabel[label].append(-1.)
                    jetParameters.append(-999.)
                else:
                    jetPredictionsPerLabel[label].append(getattr(jet, self.taggerName)[label])
                    jetParameters.append(getattr(jet, self.taggerName)['parameter'])

        self.out.fillBranch(self.inputCollection+"_"+self.taggerName+"_parameter", jetParameters)

        for label in self.profiledLabels:
            self.out.fillBranch(self.inputCollection+"_"+self.taggerName+"_"+label, jetPredictionsPerLabel[label])

        return True
