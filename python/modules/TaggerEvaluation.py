import os
import sys
import math
import json
import ROOT
import random
import numpy as np
import time
import imp

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module

from utils import getCtauLabel, getAbscissasAndWeights

class TaggerEvaluation(Module):

    def __init__(
        self,
        modelPath,
        featureDictFile,
        inputCollections = [lambda event: Collection(event, "Jet")],
        taggerName = "llpdnnx",
        evalValues = range(-1, 4),
        integrateDisplacementOrder = 2,
        globalOptions = {"isData":False},
    ):
        self.globalOptions = globalOptions
        self.inputCollections = inputCollections
        self.predictionLabels = predictionLabels
        self.evalValues = evalValues
        self.logctau = np.array(evalValues,dtype=np.float32)

        self.modelPath = os.path.expandvars(modelPath)
        feature_dict_module = imp.load_source(
            'feature_dict',
            os.path.expandvars(featureDictFile)
        )
        
        self.featureDict = feature_dict_module.featureDict
        self.predictionLabels = feature_dict_module.predictionLabels
        
        self.taggerName = taggerName

    def beginJob(self):
        pass

    def endJob(self):
        pass

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        self.setup(inputTree)


    def setupTFEval(self,tree,modelFile):
        print "Building TFEval object"
        tfEval = ROOT.TFEval()
        print "Succesfully built TFEval object"

        if (not tfEval.loadGraph(modelFile)):
            sys.exit(1)
        tfEval.addOutputNodeName("prediction")
        print "--- Model: ",modelFile," ---"
        for groupName,featureCfg in self.featureDict.iteritems():
            if featureCfg.has_key("max"):
                print "building group ... %s, shape=[%i,%i], length=%s"%(groupName,featureCfg["max"],len(featureCfg["branches"]),featureCfg["length"])
                lengthBranch = ROOT.TFEval.createAccessor(tree.arrayReader(featureCfg["length"]))
                featureGroup = ROOT.TFEval.ArrayFeatureGroup(
                    groupName,
                    len(featureCfg["branches"]),
                    featureCfg["max"],
                    lengthBranch
                )
                for branchName in featureCfg["branches"]:
                    print branchName
                    featureGroup.addFeature(ROOT.TFEval.createAccessor(tree.arrayReader(branchName)))
                tfEval.addFeatureGroup(featureGroup)
            else:
                print "building group ... %s, shape=[%i]"%(groupName,len(featureCfg["branches"]))
                featureGroup = ROOT.TFEval.ValueFeatureGroup(
                    groupName,
                    len(featureCfg["branches"])
                )
                for branchName in featureCfg["branches"]:
                    print branchName
                    featureGroup.addFeature(ROOT.TFEval.createAccessor(tree.arrayReader(branchName)))
                tfEval.addFeatureGroup(featureGroup)

        return tfEval

    def setup(self,tree):
        self.tfEvalParametric = self.setupTFEval(tree,self.modelPath)
        print "setup model successfully..."

        genFeatureGroup = ROOT.TFEval.ValueFeatureGroup("gen",1)
        self.nJets = 0
            
        genFeatureGroup.addFeature(ROOT.TFEval.PyAccessor(lambda: self.nJets, lambda jetIndex, batchIndex: self.evalValues[batchIndex%len(self.evalValues)]))
        self.tfEvalParametric.addFeatureGroup(genFeatureGroup)

        self._ttreereaderversion = tree._ttreereaderversion


    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def analyze(self, event):

        jetglobal = Collection(event, "global")
        jetglobal_indices = [global_jet.jetIdx for global_jet in jetglobal]

        jetOriginIndices = set() # superset of all indices to evaluate

        for jetCollection in self.inputCollections:
            jets = jetCollection(event)
            for ijet, jet in enumerate(jets):
                try:
                    global_jet_index = jetglobal_indices.index(jet._index)
                except ValueError:
                    print "WARNING: jet (pt: %.2f, eta: %.2f) does not have a matching global jet --> tagger cannot be evaluated!" % (jet.pt, jet.eta)
                    continue
                else:
                    global_jet = jetglobal[global_jet_index]
                    jetp4 = jet.p4()
                    #note: p4 can change in jet selection if subtracted pT<15 GeV
                    if hasattr(jet,"unselectedP4"):
                        jetp4 = jet.unselectedP4
                    if math.fabs(jetp4.Eta() - global_jet.eta) > 0.01 or \
                       math.fabs(deltaPhi(jetp4.Phi(),global_jet.phi)) > 0.01:
                           print "Warning ->> jet might be mismatched! (phi: %.2f, eta: %.2f) != (phi: %.2f, eta: %.2f)"%(jet.phi, jet.eta, global_jet.phi, global_jet.eta)
                    jetOriginIndices.add(global_jet_index)
                    setattr(jet, "globalIdx", global_jet_index)

        jetOriginIndices = list(jetOriginIndices)


        evaluationIndices = []
        for index in jetOriginIndices:
            evaluationIndices.extend([index]*len(self.evalValues))

        if event._tree._ttreereaderversion > self._ttreereaderversion:
            self.setup(event._tree)

        self.nJets = len(jetglobal)

        if len(jetOriginIndices)==0:
            return True

        evaluationIndices = np.array(evaluationIndices,np.int64)
        result = self.tfEvalParametric.evaluate(
            evaluationIndices.shape[0],
            evaluationIndices
        )

        predictionsPerIndexAndValue = {}

        for ijet,jetIndex in enumerate(jetOriginIndices):
            predictionsPerIndexAndValue[jetIndex] = {}

            for ictau,ctau in enumerate(self.evalValues):
                predictionIndex = ijet*len(self.evalValues)+ictau
                predictionsPerIndexAndCtau[jetIndex][ctau] = result.get("prediction",predictionIndex)
        for jetCollection in self.inputCollections:
            jets = jetCollection(event)

            for ijet, jet in enumerate(jets):
                taggerOutput = {}

                for ivalue, value in enumerate(self.evalValues):
                    taggerOutput[self.evalValues[ivalue]] = {}

                    for iclass, classLabel in enumerate(self.predictionLabels):
                        if hasattr(jet, "globalIdx"):
                            taggerOutput[self.evalValues[ictau]][classLabel] = \
                                    predictionsPerIndexAndCtau[jet.globalIdx][ctau][iclass]
                        else:
                            taggerOutput[self.evalValues[ictau]][classLabel] = -1.0

                setattr(jet, self.taggerName, taggerOutput)
        return True
