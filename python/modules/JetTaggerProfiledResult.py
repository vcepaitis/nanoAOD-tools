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
        kinds = ['single','ratio','avg'],
        maxOnly = True,
        globalOptions={"isData": False}
    ):
        self.globalOptions = globalOptions
        self.taggerName = taggerName
        self.outputName = outputName
        self.inputCollection = inputCollection
        self.profiledLabels = profiledLabels
        self.kinds = kinds
        self.maxOnly = maxOnly

    def beginJob(self):
        pass

    def endJob(self):
        pass

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        if not self.maxOnly:
            self.out.branch("n"+self.outputName+"_"+self.taggerName,"I")
        for label in self.profiledLabels:
            for k in self.kinds:
                if self.maxOnly:
                    self.out.branch(self.outputName+"_"+self.taggerName+"_"+k+"_"+label,"F")
                    self.out.branch(self.outputName+"_"+self.taggerName+"_"+k+"_"+label+"_param","F")

                else:
                    self.out.branch(self.outputName+"_"+self.taggerName+"_"+k+"_"+label,"F",lenVar="n"+self.outputName+"_"+self.taggerName)
                    self.out.branch(self.outputName+"_"+self.taggerName+"_"+k+"_"+label+"_param","F",lenVar="n"+self.outputName+"_"+self.taggerName)

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""

        jets = self.inputCollection(event)

        taggerResults = {}
        taggerParameters = {}
        for k in self.kinds:
            taggerResults[k] = {label: [-1.]*len(jets) for label in self.profiledLabels}
            taggerParameters[k] = {label: [-1.]*len(jets) for label in self.profiledLabels}
        
        hasTagger = False
        for ijet, jet in enumerate(jets):
            if not hasattr(jet, self.taggerName):
                continue
            hasTagger = True
            predictions = getattr(jet,self.taggerName)
            for label in self.profiledLabels:
                for k in self.kinds:
                    taggerResults[k][label][ijet] = predictions[k][label]['output']
                    taggerParameters[k][label][ijet] = predictions[k][label]['parameter']

        if len(jets)>0 and not hasTagger:
            print "WARNING - no jet in the event has the ", self.taggerName, " result stored"

        
        if self.maxOnly:
            #save only the maximum over all jets
            for label in self.profiledLabels:
                for k in self.kinds:
                    pMax = -1
                    paramMax = -10
                    for i in range(len(taggerResults[k][label])):
                        if pMax<taggerResults[k][label][i]:
                            pMax = taggerResults[k][label][i]
                            paramMax = taggerParameters[k][label][i]
                    self.out.fillBranch(self.outputName+"_"+self.taggerName+"_"+k+"_"+label,pMax)
                    self.out.fillBranch(self.outputName+"_"+self.taggerName+"_"+k+"_"+label+"_param",paramMax)
        else:
            #save all jets
            self.out.fillBranch("n"+self.outputName+"_"+self.taggerName,len(jets))
            for label in self.profiledLabels:
                for k in self.kinds:
                    self.out.fillBranch(self.outputName+"_"+self.taggerName+"_"+k+"_"+label,taggerResults[k][label])
                    self.out.fillBranch(self.outputName+"_"+self.taggerName+"_"+k+"_"+label+"_param",taggerParameters[k][label])
            
        return True
