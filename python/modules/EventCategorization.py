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
        outputName=None ,  
        looseLeptons = None , 
        jetsCollection = None ,
        taggerName="llpdnnx",
        #jetLabels=['LLP_Q','LLP_E','LLP_MU','LLP_TAU','LLP_QE','LLP_QMU','LLP_QTAU'], 
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
        self.globalOptions=globalOptions
        self.outputName=outputName
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
        self.out.branch(self.outputName+"_taggerBestOutputValue" , "F", lenVar="n"+self.outputName ) 
        self.out.branch(self.outputName+"_taggerBestOutputLabel", "F", lenVar="n"+self.outputName )      
        self.out.branch(self.outputName+"_taggerBestOutputParameter", "F" , lenVar="n"+self.outputName)     
        self.out.branch(self.outputName+"_taggerBestJet_deltaR", "F" , lenVar="n"+self.outputName)
        self.out.branch(self.outputName+"_outputSum", "F", lenVar="n"+self.outputName) 
        self.out.branch(self.outputName+"_taggerJetsPt", "F", lenVar="n"+self.outputName) 
        self.out.branch(self.outputName+"_allCategories", "I")
        if self.globalOptions["isSignal"]:
            self.out.branch(self.outputName+"_taggerBestOutputTruth","F",lenVar="n"+self.outputName)    
        self.out.branch("n"+self.outputName, "I") 
    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass
     
    def analyze(self, event):
        ## you need to add the variables you want to store and to fill the branches afterword.
        jets = self.jetsCollection(event) 
        looseLeptons = self.looseLeptons(event)
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
                'LLP_Q'  : 0 , 
                'LLP_QMU': 1 , 
                'LLP_QE' : 2 ,
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
            indexFlag = []
            indexFlag.append(-10.)
            

        ## highest probability jet categorisation

        # Take leading jet by default and first ordered label
         
        bestJet = []
        bestIndex = []
      	bestValue = []
      	bestParam = [] 
      	bestDict = []
      	bestJet.append(-10.)
       
        bestIndex.append(-10.)
        bestValue.append(-10. )
        bestParam.append(-10.)
        bestDict.append(-10.)
        
        for ijet, jet in enumerate(jets):
            taggerOutput = getattr(jet, self.taggerName)
            for label in self.jetLabels:
                 if bestValue[0] < taggerOutput[label] :
                    bestIndex[0] = dict[label]
                    bestValue[0] = taggerOutput[label] 
                    bestParam[0] = taggerOutput['parameter']
                    bestJet[0] = jet 
		    bestDict[0] = taggerOutput 

	# looking for second llp proba  jet. 
	secondValue = -10.
	secondIndex = -10
	secondParam = -10.
	secondJet = -10. 
	secondDict = {}
	if len(jets) > 1 :
           for ijet, jet in enumerate(jets):
             if jet !=  bestJet[0]:
               taggerOutput = getattr(jet, self.taggerName)
               for label in self.jetLabels:
                 if secondValue < taggerOutput[label] and taggerOutput[label] < bestValue[0]  :
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
		    if  len(bestJet) == 1  and bestJet[0]!= -10. : 
                         flagValue = getattr(jetOrigin[bestJet[0]._index], originFlag)
                         if (flagValue > 0.5):
                              indexFlag[0]  = dictTruth[k]
		    elif len(bestJet) > 1 and bestJet[1]!= -10.  : 
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
           if not jet == -10. : 
              for label in self.jetLabels : 
	       	  output += getattr(jet, self.taggerName)[label]
	      outputSum.append(output)
             

 
        nLLP = 0 
        deltaRllpJet = [] 
        if len(looseLeptons) == 0 : 
           for l in bestJet : 
              deltaRllpJet.append( -1.)
        else : 
           for j in bestJet : 
             if not j == -10. : 
                 deltaRllpJet.append(deltaR(looseLeptons[0], j))



        for i , o in enumerate(outputSum) : 
            if o > 0.5  :
              nLLP +=1 
       
 
        nleptons = 0 
	if len(looseLeptons) == 0 :
           nleptons = 1 
        else : 
           nleptons = 2 

        resolved = 0
        merged = 0

        if nLLP == 1 : 
          if bestIndex[0] == 0 : resolved = 1 
          elif bestIndex[0] == 1 or bestIndex[0] == 2 : merged = 1

        elif  nLLP == 2 :
          sum_ = 2 
          for i in bestIndex : 
            sum_ += i 
          resolved = sum_

        # binning the categories. 
        xbin = -1 
        if nleptons == 2 and nLLP == 1 and  resolved == 1 : 
      
          xbin = 1 

        elif nleptons == 2 and nLLP == 1 and  merged  == 1 : 

          xbin  = 2 
 
             
        elif nleptons == 2 and nLLP == 2 and  resolved  == 2 : 

          #xbin  = 3   # merging the 2 jets with 1 jet diplepton events 
 	  xbin   = 2

        elif nleptons == 2 and nLLP == 2 and   (resolved  == 3 or resolved == 4) :

          #xbin = 4
	  xbin = 1
 
        elif nleptons == 1 and nLLP == 1 and  resolved  == 1 : 

          xbin  = 3

        elif nleptons == 1 and nLLP == 2 and  resolved  == 2 : 

          xbin  = 4

	tagggerJetsPt = [] 
	for index , b in enumerate(bestJet) : 
          if not b < -9. : 
 		tagggerJetsPt.append(b.pt)
      
        self.out.fillBranch(self.outputName+"_allCategories" , xbin) 
        if len(looseLeptons) == 0 or math.fabs(looseLeptons[0].dxyErr) < 1e-6 :
             self.out.fillBranch(self.outputName+"_lepton_dxy_sig",-1)
        elif math.fabs(looseLeptons[0].dxyErr) > 1e-6 and len(looseLeptons) > 0 :
             self.out.fillBranch(self.outputName+"_lepton_dxy_sig", math.fabs(looseLeptons[0].dxy)/math.fabs(looseLeptons[0].dxyErr))

        self.out.fillBranch("n"+self.outputName, nLLP)  
        self.out.fillBranch(self.outputName+"_taggerJetsPt", tagggerJetsPt)
        self.out.fillBranch(self.outputName+"_taggerBestJet_deltaR", deltaRllpJet)
        self.out.fillBranch(self.outputName+"_taggerBestOutputValue", bestValue)
        self.out.fillBranch(self.outputName+"_taggerBestOutputParameter",bestParam)
        self.out.fillBranch(self.outputName+"_taggerBestOutputLabel",bestIndex)
        self.out.fillBranch(self.outputName+"_outputSum",outputSum)
        if self.globalOptions["isSignal"]:
        	self.out.fillBranch(self.outputName+"_taggerBestOutputTruth",indexFlag)
        return True 
	

