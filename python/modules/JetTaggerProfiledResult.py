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
        for label in self.profiledLabels:
            for k in ['single','ratio','avg']:
                self.out.branch(self.outputName+"_"+self.taggerName+"_"+k+"_"+label,"F",lenVar="n"+self.outputName)
                self.out.branch(self.outputName+"_"+self.taggerName+"_"+k+"_"+label+"_param","F",lenVar="n"+self.outputName)

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""

        jets = self.inputCollection(event)

        taggerResults = {}
        taggerParameters = {}
        for k in ['single','ratio','avg']:
            taggerResults[k] = {label: [-1.]*len(jets) for label in self.profiledLabels}
            taggerParameters[k] = {label: [-1.]*len(jets) for label in self.profiledLabels}
        
        hasTagger = False
        for ijet, jet in enumerate(jets):
            if not hasattr(jet, self.taggerName):
                continue
            hasTagger = True
            predictions = getattr(jet,self.taggerName)
            for label in self.profiledLabels:
                for k in ['single','ratio','avg']:
                    taggerResults[k][label][ijet] = predictions[k][label]['output']
                    taggerParameters[k][label][ijet] = predictions[k][label]['parameter']

        if not hasTagger:
            print "WARNING - no jet in the event has the ", self.taggerName, " result stored"

        #print len(jets)
        for label in self.profiledLabels:
            for k in ['single','ratio','avg']:
                self.out.fillBranch(self.outputName+"_"+self.taggerName+"_"+k+"_"+label,taggerResults[k][label])
                self.out.fillBranch(self.outputName+"_"+self.taggerName+"_"+k+"_"+label+"_param",taggerParameters[k][label])
        #print
        return True
