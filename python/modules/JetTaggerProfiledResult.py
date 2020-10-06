import os
import sys
import math
import json
import ROOT
import numpy
import random

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module

from utils import deltaR, getCtauLabel


class JetTaggerProfiledResult(Module):

    def __init__(
        self,
        inputCollection=lambda event: Collection(event, "Jet"),
        taggerName="llpdnnx",
        outputName="selectedJets",
        profiledLabels = ['LLP_Q','LLP_E','LLP_MU','LLP_TAU','LLP_QE','LLP_QMU','LLP_QTAU'],
        globalOptions={"isData": False}
    ):
        self.globalOptions = globalOptions
        self.taggerName = taggerName
        self.outputName = outputName
        self.inputCollection = inputCollection
        self.profiledLabels = profiledLabels


    def beginJob(self):
        pass

    def endJob(self):
        pass

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        for label in self.profiledLabels+['parameter']:
            self.out.branch(self.outputName+"_"+self.taggerName+"_"+label,"F",lenVar="n"+self.outputName)

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""

        jets = self.inputCollection(event)

        taggerResults = {label: [-1.]*len(jets) for label in self.profiledLabels}
        taggerResults['parameter'] = [-10]*len(jets)
        hasTagger = False
        for ijet, jet in enumerate(jets):
            if not hasattr(jet, self.taggerName):
                continue
            hasTagger = True
            predictions = getattr(jet,self.taggerName)
            for label in self.profiledLabels+['parameter']:
                taggerResults[label][ijet] = predictions[label]

        if not hasTagger:
            print "WARNING - no jet in the event has the ", self.taggerName, " result stored"

        #print len(jets)
        for k in self.profiledLabels+['parameter']:
            #print self.outputName+"_"+self.taggerName+"_"+k,taggerResults[k]
            self.out.fillBranch(self.outputName+"_"+self.taggerName+"_"+k,taggerResults[k])
        #print
        return True
