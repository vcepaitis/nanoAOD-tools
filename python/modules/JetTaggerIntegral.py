import os
import sys
import math
import json
import ROOT
import random
import numpy as np

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module

from utils import getAbscissasAndWeights

class JetTaggerIntegral(Module):

    def __init__(
        self,
        inputCollection=lambda event: Collection(event, "Jet"),
        taggerName="llpdnnx",
        outputName="selectedJets",
        predictionLabels=["B", "C", "UDS", "G", "PU", "isLLP_QMU_QQMU", "isLLP_Q_QQ"],
        globalOptions={"isData": False},
        integrateDisplacementOrder = 3,
    ):
        self.globalOptions = globalOptions
        self.taggerName = taggerName
        self.outputName = outputName
        self.inputCollection = inputCollection
        self.predictionLabels = predictionLabels
        self.abscissas, self.weights = getAbscissasAndWeights(integrateDisplacementOrder)

    def beginJob(self):
        pass

    def endJob(self):
        pass

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        self.evalDict = {}

        file_path = "PhysicsTools/NanoAODTools/data/hnl/L0.json"
        with open(file_path) as json_file:
            L0_values = json.load(json_file)
        for sample, L0 in L0_values.iteritems():
            shortName = sample.replace('HeavyNeutrino_lljj_', '').replace('M-', 'M').replace('V-', 'V').replace('.', 'p')
            print(shortName)
            self.evalDict[shortName] = []
            for abscissa in self.abscissas:
                logDisplacement = math.log10(abscissa/L0)
                self.evalDict[shortName].append(logDisplacement)
            self.evalDict[shortName] = np.array(self.evalDict[shortName], dtype=np.float32)


            for label in self.predictionLabels:
                self.out.branch(self.outputName+"_"+shortName+"_"+label,"F",lenVar="n"+self.outputName)

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""

        jets = self.inputCollection(event)

        taggerResults = {name: {className: [-1.]*len(jets) for className in self.predictionLabels} for name in self.evalDict.keys()}
        for ijet, jet in enumerate(jets):
            if not hasattr(jet, self.taggerName):
                print "WARNING - jet ", jet, " has no ", self.taggerName, " result stored -> skip"
                continue
            predictions = getattr(jet, self.taggerName)
            for label in self.predictionLabels:
                taggerResultsPerModel = {}
                for name, logDisplacements in self.evalDict.iteritems():
                    taggerResultsPerModel[name] = 0
                    for i, logDisplacement in enumerate(logDisplacements):
                        value = predictions[logDisplacement][label]
                        taggerResultsPerModel[name] += value*self.weights[i]
                    taggerResults[name][label][ijet] = taggerResultsPerModel[name]
  
        for name in self.evalDict.keys():
            for label in self.predictionLabels:
                self.out.fillBranch(self.outputName+"_"+name+"_"+label,
                        taggerResults[name][label])

        return True
