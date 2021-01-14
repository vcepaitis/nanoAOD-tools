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

class TaggerEvaluationCache(Module):

    def __init__(
        self,
        modelPath,
        featureDictFile,
        taggerName = "llpdnnx",
        evalValues = np.linspace(-1.9,1.9,5*4),
    ):
        self.evalValues = evalValues

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
        
        self.out.branch("n"+self.taggerName+"_jets",'I')
        self.out.branch(self.taggerName+"_jets_pt","F",lenVar="n"+self.taggerName+"_jets")
        self.out.branch(self.taggerName+"_jets_eta","F",lenVar="n"+self.taggerName+"_jets")
        self.out.branch(self.taggerName+"_jets_phi","F",lenVar="n"+self.taggerName+"_jets")
        self.out.branch(self.taggerName+"_jets_jetIdx","I",lenVar="n"+self.taggerName+"_jets")
        
        self.out.branch("n"+self.taggerName+"_cached",'I')
        for label in self.predictionLabels:
            self.out.branch(self.taggerName+"_cached"+"_"+label,"F",lenVar="ncached_"+self.taggerName)


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
        self.nJets = min(6,len(jetglobal))
        
        
        evaluationIndices = []
        for ijet in range(self.nJets):
            evaluationIndices.extend([ijet]*len(self.evalValues))

        if event._tree._ttreereaderversion > self._ttreereaderversion:
            self.setup(event._tree)
        
        if self.nJets==0:
            self.out.fillBranch("n"+self.taggerName+"_jets",0)
            self.out.fillBranch(self.taggerName+"_jets_pt",[])
            self.out.fillBranch(self.taggerName+"_jets_eta",[])
            self.out.fillBranch(self.taggerName+"_jets_phi",[])
            self.out.fillBranch(self.taggerName+"_jets_jetIdx",[])
            
            self.out.fillBranch("n"+self.taggerName+"_cached",self.nJets*len(self.evalValues))
            for label in self.predictionLabels:
                self.out.fillBranch(self.taggerName+"_cached"+"_"+label,[])
            return True

        evaluationIndices = np.array(evaluationIndices,np.int64)
        result = self.tfEvalParametric.evaluate(
            evaluationIndices.shape[0],
            evaluationIndices
        )
        
        branchJetPt = np.zeros((self.nJets,),dtype=np.float32)
        branchJetEta = np.zeros((self.nJets,),dtype=np.float32)
        branchJetPhi = np.zeros((self.nJets,),dtype=np.float32)
        branchJetIdx = np.zeros((self.nJets,),dtype=np.int32)
        
        for i in range(self.nJets):
            branchJetPt[i] = jetglobal[i].pt
            branchJetPt[i] = jetglobal[i].eta
            branchJetPt[i] = jetglobal[i].phi
            branchJetPt[i] = jetglobal[i].jetIdx
            

        self.out.fillBranch("n"+self.taggerName+"_jets",self.nJets)
        self.out.fillBranch(self.taggerName+"_jets_pt",branchJetPt)
        self.out.fillBranch(self.taggerName+"_jets_eta",branchJetEta)
        self.out.fillBranch(self.taggerName+"_jets_phi",branchJetPhi)
        self.out.fillBranch(self.taggerName+"_jets_jetIdx",branchJetIdx)
        
        self.out.fillBranch("n"+self.taggerName+"_cached",self.nJets*len(self.evalValues))


        for ilabel,label in enumerate(self.predictionLabels):
            branchData = np.zeros((self.nJets*len(self.evalValues)))
            for ijet in range(self.nJets):
                for ivalue in range(len(self.evalValues)):
                    predictionIndex = ijet*len(self.evalValues)+ivalue
                    branchData[ijet*len(self.evalValues)+ivalue] = result.get("prediction",predictionIndex)[ilabel]
            self.out.fillBranch(self.taggerName+"_cached"+"_"+label,branchData)
            
        return True
