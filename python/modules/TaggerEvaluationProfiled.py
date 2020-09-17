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

class TaggerEvaluationProfiled(Module):
    def __init__(
        self,
        modelPath,
        featureDictFile,
        inputCollections = [lambda event: Collection(event, "Jet")],
        taggerName = "llpdnnx",
        evalValues = range(-1, 4),
        profiledLabelDict = {
                            'LLP_Q': ['LLP_Q', 'LLP_TAU', 'LLP_QTAU'],
                            'LLP_QE': ['LLP_E', 'LLP_QE'],
                            'LLP_QMU': ['LLP_MU', 'LLP_QMU']
                            },
        globalOptions = {"isData":False},
    ):
        self.globalOptions = globalOptions
        self.inputCollections = inputCollections
        
        self.evalValues = list(evalValues)
        self.nEvalValues = len(evalValues)
        self.profiledLabelDict = profiledLabelDict
        self.profiledLabels = []

        self.modelPath = os.path.expandvars(modelPath)
        
        feature_dict_module = imp.load_source(
            'feature_dict',
            os.path.expandvars(featureDictFile)
        )
        
        self.featureDict = feature_dict_module.featureDict
        self.predictionLabels = feature_dict_module.predictionLabels
        for _, profiledLabels in self.profiledLabelDict.iteritems():
            for profiledSubLabel in profiledLabels:
                self.profiledLabels.append(profiledSubLabel)
                if profiledSubLabel not in self.predictionLabels:
                    print "ERROR (TaggerEvaluationProfiled) - label for profiling '%s' not in predicted label list: %s"%(
                        profiledSubLabel,
                        str(self.predictionLabels)
                    )
                    sys.exit(1)
        
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
                    #print branchName
                    featureGroup.addFeature(ROOT.TFEval.createAccessor(tree.arrayReader(branchName)))
                tfEval.addFeatureGroup(featureGroup)
            else:
                print "building group ... %s, shape=[%i]"%(groupName,len(featureCfg["branches"]))
                featureGroup = ROOT.TFEval.ValueFeatureGroup(
                    groupName,
                    len(featureCfg["branches"])
                )
                for branchName in featureCfg["branches"]:
                    #print branchName
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
                    print "WARNING: jet (pt: %s, eta: %s) does not have a matching global jet --> tagger cannot be evaluated!" % (jet.pt, jet.eta)
                    continue
                else:
                    global_jet = jetglobal[global_jet_index]
                    if abs(jet.eta - global_jet.eta) > 0.01 or \
                       abs(jet.phi - global_jet.phi) > 0.01:
                           print "Warning ->> jet might be mismatched!"
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
            for jetCollection in self.inputCollections:
                jets = jetCollection(event)
                for ijet, jet in enumerate(jets):
                    taggerOutput = {}
                    for ilabel, label in enumerate(self.predictionLabels):
                        print("Jet output set to -1!")
                        taggerOutput[label] = {
                            'output': -1.0,
                            'parameter': min(self.evalValues)
                        }
                    setattr(jet, self.taggerName, taggerOutput)
            return True

        evaluationIndices = np.array(evaluationIndices,np.int64)

        result = self.tfEvalParametric.evaluate(
            evaluationIndices.shape[0],
            evaluationIndices
        )

        outputPerIndex = {}

        for ijet,jetIndex in enumerate(jetOriginIndices):
            

            maxSinglePrediction = -1
            ivalueAtMaxPrediction = 0
            
            for ilabel,label in enumerate(self.predictionLabels):
                if label not in self.profiledLabels:
                    continue
                
                for ivalue, value in enumerate(self.evalValues):
                    predictionIndex = ijet*len(self.evalValues)+ivalue
                    predictions = result.get("prediction",predictionIndex)
                    singlePrediction = predictions[ilabel]
                    if singlePrediction>maxSinglePrediction:
                        maxSinglePrediction = singlePrediction
                        ivalueAtMaxPrediction = ivalue
            
            
            maxPredictionIndex = ijet*len(self.evalValues)+ivalueAtMaxPrediction
            maxPredictions = result.get("prediction",maxPredictionIndex)

            outputPerIndex[jetIndex] = {}
            outputPerIndex[jetIndex]['parameter'] = self.evalValues[ivalueAtMaxPrediction]
            for ilabel,label in enumerate(self.predictionLabels):
                if label not in self.profiledLabels:
                    continue
                outputPerIndex[jetIndex][label] = maxPredictions[ilabel]


        # group some class outputs together
        # to be removed after retraining
        outputPerIndexGrouped = {}

        for index, output in outputPerIndex.iteritems():
            outputGrouped = {}
            for label, sublabels in self.profiledLabelDict.iteritems():
                outputGrouped[label] = 0.
                for sublabel in sublabels:
                    outputGrouped[label] += output[sublabel]
            outputGrouped['parameter'] = output['parameter']
            outputPerIndexGrouped[index] = outputGrouped

        for jetCollection in self.inputCollections:
            jets = jetCollection(event)

            for ijet, jet in enumerate(jets):
                if hasattr(jet, "globalIdx"):
                    setattr(jet, self.taggerName, outputPerIndexGrouped[jet.globalIdx])
                else:
                    print("Jet output set to -1!")
                    taggerOutput = {k:-1 for k in self.profiledLabels}
                    taggerOutput['parameter'] = min(self.evalValues)
                    setattr(jet, self.taggerName, taggerOutput)
                
        return True
