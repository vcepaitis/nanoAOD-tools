import os
import sys
import math
import json
import ROOT
import random

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from utils import deltaR, getCtauLabel

class EventCategorization(Module):
    def __init__(
        self,
        globalOptions={"isData":False, "isSignal":False},
        outputName=None,
        tightLeptons=None,
        looseLeptons=None,
        jetsCollection=None,
        taggerName="llpdnnx",
        maxDeltaR=1.3,
        jetLabels=['LLP_Q','LLP_QE','LLP_QMU'],
        flags={
            'isPrompt_MU': ['isPrompt_MU'],
            'isPrompt_E': ['isPrompt_E'],
            'isPrompt_TAU': ['isPrompt_TAU'],
            'isB': ['isB', 'isBB', 'isGBB', 'isLeptonic_B', 'isLeptonic_C'],
            'isC': ['isC', 'isCC', 'isGCC'],
            'isUDS': ['isS', 'isUD'],
            'isG': ['isG'],
            'isPU': ['isPU'],
            'isLLP_Q': ['isLLP_RAD','isLLP_Q','isLLP_QQ'],
            'isLLP_MU': ['isLLP_MU'],
            'isLLP_QMU': ['isLLP_QMU', 'isLLP_QQMU'],
            'isLLP_E': ['isLLP_E'],
            'isLLP_QE': ['isLLP_QE','isLLP_QQE'],
            'isLLP_TAU': ['isLLP_TAU'],
            'isLLP_QTAU': ['isLLP_QTAU','isLLP_QQTAU'],
            'isUndefined': ['isUndefined'] ,
            'isLLP_B':['isLLP_B' , 'isLLP_BB' ,'isLLP_BMU' ,'isLLP_BE', 'isLLP_BBMU', 'isLLP_BBE'],
        },
    ):
        self.globalOptions = globalOptions
        self.outputName = outputName
        self.tightLeptons = tightLeptons
        self.looseLeptons = looseLeptons
        self.jetsCollection = jetsCollection
        self.taggerName = taggerName
        self.jetLabels = jetLabels
        self.maxDeltaR = maxDeltaR
        self.flags = flags

    def beginJob(self):
        pass

    def endJob(self):
        pass

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        self.out.branch("n"+self.outputName, "I")
        self.out.branch(self.outputName+"_taggerBestOutputValue" , "F", lenVar="n"+self.outputName )
        self.out.branch(self.outputName+"_taggerBestOutputLabel", "F", lenVar="n"+self.outputName )
        self.out.branch(self.outputName+"_taggerBestOutputParameter", "F" , lenVar="n"+self.outputName)
        self.out.branch(self.outputName+"_taggerBestJet_deltaR", "F" , lenVar="n"+self.outputName)
        self.out.branch(self.outputName+"_outputSum", "F", lenVar="n"+self.outputName)
        self.out.branch(self.outputName+"_taggerJetsPt", "F", lenVar="n"+self.outputName)
        self.out.branch(self.outputName+"_allCategories", "I")
        self.out.branch(self.outputName+"_WCandidateMass", "F")
        self.out.branch(self.outputName+"_HNLCandidateMass", "F")
        if self.globalOptions["isSignal"]:
            self.out.branch(self.outputName+"_taggerBestOutputTruth","F",lenVar="n"+self.outputName)
            self.out.branch(self.outputName+"_index_truth", "I")
            self.out.branch("n"+self.outputName+"_TrueJets", "I")
            self.out.branch("n"+self.outputName+"_TrueLepJets", "I")
            self.out.branch("n"+self.outputName+"_TrueResJets", "I")

        self.out.branch("n"+self.outputName, "I")

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def analyze(self, event):
        ## you need to add the variables you want to store and to fill the branches afterword.
        jets = self.jetsCollection(event)
        looseLeptons = self.looseLeptons(event)
        tightLeptons = self.tightLeptons(event)

        '''
        dict = {'LLP_Q': 0,
                'LLP_MU': 1,
                'LLP_QMU': 2,
                'LLP_E': 3,
                'LLP_QE': 4,
                'LLP_TAU': 5,
                'LLP_QTAU': 6,
               }
        '''

        dict = {
                'LLP_Q'  : 0,
                'LLP_QMU': 1,
                'LLP_QE' : 2,
               }

        dictTruth = {
                      'isLLP_Q': 0,
                      'isLLP_MU': 1,
                      'isLLP_QMU': 2,
                      'isLLP_E': 3,
                      'isLLP_QE': 4,
                      'isLLP_TAU': 5,
                      'isLLP_QTAU': 6,
                      'isPU': 7,
                      'isUDS': 8,
                      'isG': 9,
                      'isB': 10,
                      'isC': 11,
                      'isPrompt_MU': 12,
                      'isPrompt_E': 13,
                      'isPrompt_TAU': 14,
                      'isUndefined': 15,
                      'isLLP_B' : 16 ,
                    }


        # only for MC: truth labels
        if self.globalOptions["isSignal"]:
            jetOrigin = Collection(event, "jetorigin")
            indexFlags = []
            indexFlags.append(-10.)


        ## highest probability jet categorisation
        # Take leading jet by default and first ordered label

        bestJets = []
        bestIndices = []
      	bestValues = []
      	bestParams = []
      	bestDicts = []
      	bestJets.append(-10.)

        bestIndices.append(-10.)
        bestValues.append(-10.)
        bestParams.append(-10.)
        bestDicts.append(-10.)

        for ijet, jet in enumerate(jets):
            taggerOutput = getattr(jet, self.taggerName)
            for label in self.jetLabels:
                 if bestValues[0] < taggerOutput[label]:
                    bestIndices[0] = dict[label]
                    bestValues[0] = taggerOutput[label]
                    bestParams[0] = taggerOutput['parameter']
                    bestJets[0] = jet
        bestDicts[0] = taggerOutput

        # looking for second llp proba  jet.
        secondValue = -10.
        secondIndex = -10
        secondParam = -10.
        secondJet = -10.
        secondDict = {}
        if len(jets) > 1:
            for ijet, jet in enumerate(jets):
                if jet != bestJets[0]:
                    taggerOutput = getattr(jet, self.taggerName)
                    for label in self.jetLabels:
                        if secondValue < taggerOutput[label] and taggerOutput[label] < bestValues[0]:
                            secondIndex = dict[label]
                            secondValue = taggerOutput[label]
                            secondParam = taggerOutput['parameter']
                            secondJet  = jet
                            secondDict = taggerOutput
        bestValues.append(secondValue)
        bestIndices.append(secondIndex)
        bestParams.append(secondParam)
        bestJets.append(secondJet)
        bestDicts.append(secondDict)


        # looking for the truth labeling for each best jet per label. only for MC
        if self.globalOptions["isSignal"]:
            for k in sorted(self.flags.keys()):
                # loop over general truth flags
                # loop over sub true flags

                for originFlag in self.flags[k]:
                    if len(bestJets) == 1 and bestJets[0]!= -10.:
                        flagValue = getattr(jetOrigin[bestJet[0]._index], originFlag)
                        if (flagValue > 0.5):
                            indexFlags[0] = dictTruth[k]
                    elif len(bestJets) > 1 and bestJets[1]!= -10.:
                        flagValue = getattr(jetOrigin[bestJets[0]._index], originFlag)
                        flagValue2 = getattr(jetOrigin[bestJets[1]._index], originFlag)
                        if (flagValue > 0.5):
                            indexFlags[0] = dictTruth[k]
                        if (flagValue2 > 0.5):
                            indexFlags.append(dictTruth[k])

        # sum of all llp proba for the best LLP jets.
        outputSums = []

        for i, jet in enumerate(bestJets):
            output = 0.
            if not jet == -10.:
                for label in self.jetLabels:
           	        output += getattr(jet, self.taggerName)[label]
                outputSums.append(output)


        nLLP = 0
        deltaRllpJet = []
        if len(looseLeptons) == 0:
           for l in bestJets:
              deltaRllpJet.append( -1.)
        else:
           for j in bestJets:
             if not j == -10.:
                 deltaRllpJet.append(deltaR(looseLeptons[0], j))

        for i, o in enumerate(outputSums):
            if o > 0.5:
              nLLP +=1

        nleptons = 0
        if len(looseLeptons) == 0 :
            nleptons = 1
        else:
            nleptons = 2

        resolved = 0
        merged = 0

        if nLLP == 1:
          if bestIndices[0] == 0: resolved = 1
          elif bestIndices[0] == 1 or bestIndices[0] == 2: merged = 1

        elif nLLP == 2:
          sum_ = 2
          for i in bestIndices:
            sum_ += i
          resolved = sum_

        # binning the categories.
        xbin = -1
        if nleptons == 2 and nLLP == 1 and  resolved == 1:
            xbin = 1

        elif nleptons == 2 and nLLP == 1 and  merged == 1:
            xbin = 2

        elif nleptons == 2 and nLLP == 2 and  resolved  == 2:
            #xbin  = 3   # merging the 2 jets with 1 jet diplepton events
    	    xbin = 2

        elif nleptons == 2 and nLLP == 2 and (resolved  == 3 or resolved == 4):
            #xbin = 4
            xbin = 1

        elif nleptons == 1 and nLLP == 1 and  resolved  == 1:
            xbin = 3

        elif nleptons == 1 and nLLP == 2 and  resolved  == 2:
            xbin = 4

        taggerJetsPt = []
        for index, b in enumerate(bestJets) :
            if not b < -9.:
     		     taggerJetsPt.append(b.pt)

        self.out.fillBranch(self.outputName+"_allCategories", xbin)
        self.out.fillBranch("n"+self.outputName, nLLP)
        self.out.fillBranch(self.outputName+"_taggerJetsPt", taggerJetsPt)
        self.out.fillBranch(self.outputName+"_taggerBestJet_deltaR", deltaRllpJet)
        self.out.fillBranch(self.outputName+"_taggerBestOutputValue", bestValues)
        self.out.fillBranch(self.outputName+"_taggerBestOutputParameter",bestParams)
        self.out.fillBranch(self.outputName+"_taggerBestOutputLabel",bestIndices)
        self.out.fillBranch(self.outputName+"_outputSum",outputSums)
        if self.globalOptions["isSignal"]:
            self.out.fillBranch(self.outputName+"_taggerBestOutputTruth",indexFlags)

            nTrueResJets = 0
            nTrueLepJets = 0
            category_index_truth = -1

            for jet in jets:
                if jet.isLLP_Q > 0 or jet.isLLP_QTAU > 0 or jet.isLLP_TAU > 0:
                    nTrueResJets += 1
                if jet.isLLP_QMU > 0 or jet.isLLP_QE > 0 or jet.isLLP_MU > 0 or jet.isLLP_E > 0:
                    nTrueLepJets += 1

            nTrueJets = nTrueLepJets+nTrueResJets

            if nTrueJets == 1:
                if nTrueResJets == 1:
                    category_index_truth = 1
                elif nTrueLepJets == 1:
                    category_index_truth = 2
            elif nTrueJets > 1:
                if nTrueLepJets > 0:
                    category_index_truth = 4
                else:
                    category_index_truth = 3
            self.out.fillBranch("n"+self.outputName+"_TrueJets", nTrueJets)
            self.out.fillBranch("n"+self.outputName+"_TrueLepJets", nTrueLepJets)
            self.out.fillBranch("n"+self.outputName+"_TrueResJets", nTrueResJets)
            self.out.fillBranch(self.outputName+"_index_truth", category_index_truth)

        HNLCandidateLorentzVector = ROOT.TLorentzVector(0,0,0,0)
        WCandidateLorentzVector = ROOT.TLorentzVector(0,0,0,0)

        if xbin > 0:
            HNLCandidateLorentzVector += bestJets[0].p4()
            WCandidateLorentzVector += bestJets[0].p4()
        for lepton in looseLeptons:
            HNLCandidateLorentzVector += lepton.p4()
            WCandidateLorentzVector += lepton.p4()
        for lepton in tightLeptons:
            WCandidateLorentzVector += lepton.p4()

        WCandidateMass = WCandidateLorentzVector.M()

        if xbin > 2:
            HNLCandidateMass = bestJets[0].mass
        else:
            HNLCandidateMass = HNLCandidateLorentzVector.M()

        self.out.fillBranch(self.outputName+"_WCandidateMass", WCandidateMass)
        self.out.fillBranch(self.outputName+"_HNLCandidateMass", HNLCandidateMass)

        return True
