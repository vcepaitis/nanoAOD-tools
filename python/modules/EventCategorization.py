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
        muonsTight = None , 
        electronsTight = None , 
        muonsLoose = None , 
        electronsLoose = None , 
        outputName=None ,  
        looseLeptons = None , 
        jetsCollection = None ,
        taggerName="llpdnnx",
        jetLabels=['LLP_Q','LLP_E','LLP_MU','LLP_TAU','LLP_QE','LLP_QMU','LLP_QTAU'], 
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
        self.globalOptions=globalOptions
        self.outputName=outputName
        self.muonsTight = muonsTight  
        self.electronsTight = electronsTight
        self.muonsLoose =  muonsLoose
        self.electronsLoose =  electronsLoose
        self.looseLeptons = looseLeptons
        self.jetsCollection = jetsCollection
        self.taggerName = taggerName 
        self.jetLabels = jetLabels
        self.flags = flags
       
    def beginJob(self):
        pass
        
    def endJob(self):
        pass
        
    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        self.out.branch("n"+self.outputName, "I") 	
        self.out.branch(self.outputName+"_lepton_dxy_sig","F")
        self.out.branch(self.outputName+"_taggerBestJet_deltaR", "F")
        self.out.branch(self.outputName+"_taggerBestOutputValue" , "F", lenVar="n"+self.outputName ) 
        self.out.branch(self.outputName+"_taggerBestOutputLabel", "F", lenVar="n"+self.outputName )      
        self.out.branch(self.outputName+"_taggerBestOutputParameter", "F" , lenVar="n"+self.outputName)     
	self.out.branch(self.outputName+"_outputSum", "F", lenVar="n"+self.outputName) 
        self.out.branch(self.outputName+"_muonmuon", "I") 
        self.out.branch(self.outputName+"_electronelectron", "I")
        self.out.branch(self.outputName+"_muonelectron", "I")
        self.out.branch(self.outputName+"_electronmuon", "I")
        self.out.branch(self.outputName+"_muonjets", "I")
        self.out.branch(self.outputName+"_electronjets", "I")
        self.out.branch(self.outputName+"_allCategories", "I")
        if self.globalOptions["isSignal"]:
            self.out.branch(self.outputName+"_taggerBestOutputTruth","F",lenVar="n"+self.outputName)    
        self.out.branch("n"+self.outputName, "I") 
    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass
     
    def analyze(self, event):
        ## you need to add the variables you want to store and to fill the branches afterword.
        jets = self.jetsCollection(event) 
        tightMuon = self.muonsTight(event) 
        tightElectron = self.electronsTight(event)
        looseMuons = self.muonsLoose(event)
        looseElectrons = self.electronsLoose(event)
        looseLeptons = self.looseLeptons(event)
        muonmuon = 0 
        electronelectron = 0
        muonelectron = 0 
        electronmuon = 0
        muonjets = 0 
        electronjets = 0
        deltaRllpJet = -1.

        dict = {'LLP_Q': 0,
                'LLP_MU': 1,
                'LLP_QMU': 2,
                'LLP_E': 3,
                'LLP_QE': 4,
                'LLP_TAU': 5,
                'LLP_QTAU': 6,
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
        
        ## flavour categorisation :
        if len(tightMuon) == 1 and len(looseMuons) == 1:
            muonmuon = 1
             
        elif len(tightElectron) == 1 and len(looseElectrons) == 1:
            electronelectron = 1 

        elif len(tightMuon) == 1 and len(looseElectrons) == 1:
            muonelectron = 1 

        elif len(tightElectron) and len(looseMuons) == 1:
            electronmuon = 1

        elif len(tightMuon) == 1 and len(tightElectron) == 1:
            if tightMuon[0].pt > tightElectron[0].pt:
                muonelectron = 1
            else:
                electronmuon = 1

        elif len(tightMuon) == 1 and len(tightElectron) == 0 and len(looseMuons) == 0 and len(looseElectrons) == 0:
            muonjets = 1 
        elif len(tightMuon) == 0 and len(tightElectron) == 1 and len(looseMuons) == 0 and len(looseElectrons) == 0:
            electronjets = 1


        if muonmuon or muonelectron or muonjets:
            if not event.IsoMuTrigger_flag:
                return False
        elif electronelectron or electronmuon or electronjets:
            if not event.IsoElectronTrigger_flag:
                return False



        # only for MC: truth labels
        if self.globalOptions["isSignal"]:
            jetOrigin = Collection(event, "jetorigin")
            indexFlag = []
            indexFlag.append(-10.)
            

        ## highest probability jet categorisation

        # Take leading jet by default and first ordered label
         
	bestJet = []
        bestIndex = []
	bestValue = []
	bestParam = [] 
	bestDict = []
	bestJet.append(jets[0])
        bestIndex.append( dict[self.jetLabels[0]])
	bestValue.append(getattr(bestJet[0], self.taggerName)[self.jetLabels[0]] )
	bestParam.append(getattr(bestJet[0], self.taggerName)['parameter'])
	bestDict.append(getattr(bestJet[0], self.taggerName))
    	
        for ijet, jet in enumerate(jets):
            taggerOutput = getattr(jet, self.taggerName)
            for label in self.jetLabels:
                 if bestValue[0] < taggerOutput[label]:
                    bestIndex[0] = dict[label]
                    bestValue[0] = taggerOutput[label] 
                    bestParam[0] = taggerOutput['parameter']
                    bestJet[0] = jet 
		    bestDict[0] = taggerOutput 

	# looking for second llp proba  jet. 
	secondValue = 0.
	secondIndex = -10
	secondParam = -10.
	secondJet = 0. 
	secondDict = {}
	if len(jets) > 1 :
           for ijet, jet in enumerate(jets):
             if jet !=  bestJet[0]:
               taggerOutput = getattr(jet, self.taggerName)
               for label in self.jetLabels:
                 if secondValue < taggerOutput[label] and taggerOutput[label] < bestValue[0] :
                    secondIndex= dict[label]
                    secondValue = taggerOutput[label] 
                    secondParam = taggerOutput['parameter']
                    secondJet  = jet 
                    secondDict = taggerOutput 
           bestValue.append(secondValue)
           bestIndex.append(secondIndex)
           bestParam.append(secondParam)
           bestJet.append(secondJet)
           bestDict.append(secondDict)


        # looking for the truth labeling for each best jet per label. only for MC
	if self.globalOptions["isSignal"]:
           for k in sorted(self.flags.keys()):
                # loop over general truth flags 
                # loop over sub true flags 

                for originFlag in self.flags[k]:
		    if  len(bestJet) == 1 : 
                         flagValue = getattr(jetOrigin[bestJet[0]._index], originFlag)
                         if (flagValue > 0.5):
                              indexFlag[0]  = dictTruth[k]
		    elif len(bestJet) > 1 : 
                      flagValue = getattr(jetOrigin[bestJet[0]._index], originFlag)
                      flagValue2 = getattr(jetOrigin[bestJet[1]._index], originFlag)
                      if (flagValue > 0.5):
                          indexFlag[0]  = dictTruth[k]
                      if(flagValue2  > 0.5) :
			indexFlag.append(dictTruth[k])
	   
                               
	# sum of all llp proba for the best LLP jets.
        outputSum = []

        for i , jet in enumerate(bestJet) : 
           output = 0.
           for label in self.jetLabels : 
		output += getattr(jet, self.taggerName)[label]
	   outputSum.append(output)
                  

        # higher level categorisation 
        # number of leptons per event. 
        nleptons = 0 
	if len(looseLeptons) == 0 :
           nleptons = 1 
        else : 
           nleptons = 2 

        nLLP = 0 
        for o in outputSum : 
            if o > 0.5 :
             nLLP +=1 

        # resolved value will depend on the  nb of LLP and their labels 

        # resolved =
        # 1  // 1 llp_Q 
        # 2  // 2 LLP_Q 
        # 3  // 1 LLP_Q + 1 LLP_MU 
        # 4  // 1 LLP_Q + 1 LLP_QMU. 

        resolved = 0
        merged = 0
        bkgd = 0 

        if nLLP == 1 : 
          if bestIndex[0] == 0 : 
              resolved = 1 
          elif bestIndex[0] == 2 or bestIndex[0] == 4 or bestIndex[0] == 6 : 
	      merged = 1
          else : 
              bkgd = 1 

        elif  nLLP == 2 : 
           sum = 2 
           for i in bestIndex : 
               sum += 1 
           resolved = sum 

        # binning the categories. 
        xbin = -1 
        if nleptons == 2 and nLLP == 1 and  resolved == 1 : 
      
          xbin = 1 

        elif nleptons == 2 and nLLP == 1 and  merged  == 1 : 

          xbin  = 2 
 
             
        elif nleptons == 2 and nLLP == 2 and  resolved  == 2 : 

          xbin  = 3 

        elif nleptons == 2 and nLLP == 2 and  resolved  == 3 : 

          xbin  = 4 


        elif nleptons == 2 and nLLP == 2  and  resolved  == 4  : 

          xbin  = 5 

        elif nleptons == 1 and nLLP == 1 and  resolved  == 1 : 

          xbin  = 6 

        elif nleptons == 1 and nLLP == 2 and  resolved  == 2 : 

          xbin  = 7 
    
     
	# you need now the deltaR of the second muon with the llp jet.  0 and >1 loose lepton categories
        if len(looseLeptons) == 0: 
           deltaRllpJet = -1. 
        else :  
           deltaRllpJet = deltaR(looseLeptons[0], bestJet[0])  
        
        self.out.fillBranch(self.outputName+"_muonmuon", muonmuon) 
        self.out.fillBranch(self.outputName+"_electronelectron", electronelectron)
        self.out.fillBranch(self.outputName+"_muonelectron", muonelectron)
        self.out.fillBranch(self.outputName+"_electronmuon", electronmuon)
        self.out.fillBranch(self.outputName+"_muonjets", muonjets)
        self.out.fillBranch(self.outputName+"_electronjets", electronjets)
        self.out.fillBranch(self.outputName+"_electronjets", electronjets)
        self.out.fillBranch(self.outputName+"_allCategories" , xbin)
 
        if len(looseLeptons) == 0 or math.fabs(looseLeptons[0].dxyErr) < 1e-6 :
             self.out.fillBranch(self.outputName+"_lepton_dxy_sig",-1
)
        elif math.fabs(looseLeptons[0].dxyErr) > 1e-6 and len(looseLeptons) > 0 :
             self.out.fillBranch(self.outputName+"_lepton_dxy_sig", math.fabs(looseLeptons[0].dxy)/math.fabs(looseLeptons[0].dxyErr))

        self.out.fillBranch("n"+self.outputName, len(bestJet))  
        self.out.fillBranch(self.outputName+"_taggerBestOutputValue", bestValue)
        self.out.fillBranch(self.outputName+"_taggerBestOutputParameter",bestParam)
        self.out.fillBranch(self.outputName+"_taggerBestOutputLabel",bestIndex)
        self.out.fillBranch(self.outputName+"_outputSum",outputSum)
        if self.globalOptions["isSignal"]:
        	self.out.fillBranch(self.outputName+"_taggerBestOutputTruth",indexFlag)
        return True 
	

