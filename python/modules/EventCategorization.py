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
        globalOptions={"isData":False}, 
        muonsTight = None , 
        electronsTight = None , 
        muonsLoose = None , 
        electronsLoose = None , 
        outputName=None ,  
	looseLeptons = None , 
        jetsCollection = None ,
        taggerName="llpdnnx",
        jetLabels=['LLP_Q','LLP_MU','LLP_E','LLP_TAU'], 
	flags={
            'isB': ['isB', 'isBB', 'isGBB', 'isLeptonic_B', 'isLeptonic_C'],
            'isC': ['isC', 'isCC', 'isGCC'],
            'isUDS': ['isS', 'isUD'],
            'isG': ['isG'],
            'isPU': ['isPU'],
            'isLLP_Q': ['isLLP_RAD','isLLP_Q','isLLP_QQ'],
            'isLLP_MU': ['isLLP_MU','isLLP_QMU','isLLP_QQMU'],
            'isLLP_E': ['isLLP_E','isLLP_QE','isLLP_QQE'],
            'isLLP_TAU': ['isLLP_TAU','isLLP_QTAU','isLLP_QQTAU'],
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
	self.out.branch(self.outputName+"_TaggerBestJet_deltaR", "F")
        self.out.branch(self.outputName+"_TaggerBestOutputValue_truth","F" )	
	self.out.branch(self.outputName+"_TaggerBestOutputValue" , "F" ) 
	self.out.branch(self.outputName+"_TaggerBestOutputLabel", "F" )		 
        self.out.branch(self.outputName+"_muonmuon1jet", "F") 
	self.out.branch(self.outputName+"_muonmuonNjets", "F")
	self.out.branch(self.outputName+"_electronelectron1jet", "F")
	self.out.branch(self.outputName+"_electronelectronNjets", "F")
	self.out.branch(self.outputName+"_muonelectron1jet", "F")
	self.out.branch(self.outputName+"_muonelectronNjets", "F")
	self.out.branch(self.outputName+"_electronmuon1jet", "F")
	self.out.branch(self.outputName+"_electronmuonNjets", "F")
	self.out.branch(self.outputName+"_TaggerBestOutputValue" , "F" ) 
	self.out.branch(self.outputName+"_TaggerBestOutputLabel", "F" )		 
        for k in sorted(self.flags.keys()):
            for originFlag in self.flags[k]:
		self.out.branch(self.outputName+"_truth_"+originFlag+"_flag", "F" )		 	
        for label in self.jetLabels:
        	self.out.branch(self.outputName+"_"+label+"_lepton_deltaR","F")
		self.out.branch(self.outputName+"_"+label+"_lepton_dxy_sig","F")
        	self.out.branch(self.outputName+"_"+label+"_jet_pt","F")
        	self.out.branch(self.outputName+"_"+label+"_jet_output","F")
        	self.out.branch(self.outputName+"_"+label+"_jet_parameter","F")
        	## declare branches here.
         
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
	muonmuon1jet = 0 
	muonmuonNjets = 0
	electronelectron1jet = 0
	electronelectronNjets = 0
	muonelectron1jet = 0 
	electronmuon1jet = 0
	muonelectronNjets = 0
	electronmuonNjets = 0
	deltaRBestLabel = -1.
	dict = {'LLP_Q' : 0 ,'LLP_MU' : 1,'LLP_E': 2 ,'LLP_TAU' : 3 }
	dictTruth = { 'isB': 0,
            'isC': 7,
            'isUDS': 6,
            'isG': 5,
            'isPU': 4,
            'isLLP_Q': 0,
            'isLLP_MU': 1,
            'isLLP_E': 2,
            'isLLP_TAU': 3,
       }
	
       	## flavour categorisation :
	
	if len(tightMuon)  == 1 and len(looseMuons) == 1 and  len(jets) == 1 :
		muonmuon1jet = 1
	elif len(tightMuon) == 1 and len(looseMuons) == 1  and len(jets) > 1 :
		muonmuonNjets = 1
	elif len(tightElectron) == 1 and len(looseElectrons) == 1 and len(jets) == 1 : 
		electronelectron1jet = 1 
	elif len(tightElectron) == 1 and len(looseElectrons) == 1 and  len(jets) > 1 : 
		electronelectronNjets = 1 
	elif len(tightMuon) == 1 and len(looseElectrons) == 1 and len(jets) == 1 : 
		muonelectron1jet = 1 
	elif len(tightMuon) == 1 and len(tightElectron) == 1  and len(jets) == 1 : 
		if tightMuon[0].pt > tightElectron[0].pt :
			muonelectron1jet = 1 
		else : 
			electronmuon1jet = 1

	elif len(tightMuon) == 1 and  len(looseElectrons) == 1  and len(jets) > 1 : 
		muonelectronNjets = 1
	elif len(tightMuon) == 1 and  len(tightElectron) == 1 and len(jets) > 1 : 
		if tightMuon[0].pt > tightElectron[0].pt :
			muonelectronNjets = 1
		else: 
			electronmuonNjets = 1
 
	## highest probability jet categorisation
	jetOrigin = Collection(event, "jetorigin")	
        bestJetsPerLabel = {}
	bestIndexPerLabel = {}
	indexFlag = {}

        for label in self.jetLabels:
            bestJetsPerLabel[label] = jets[0]
	
        for ijet , jet in enumerate(jets):
            taggerOutput = getattr(jet, self.taggerName)
            for label in self.jetLabels:
                if getattr(bestJetsPerLabel[label],self.taggerName)[label]['output']<taggerOutput[label]['output']:
                    bestJetsPerLabel[label] = jet
	bestResult = 0.00
	for label in self.jetLabels:
            jet = bestJetsPerLabel[label]
            ### looking for the truth labeling for each best jet per label. 
            for k in sorted(self.flags.keys()):
                   flavorFlag = 0.
                   for originFlag in self.flags[k]:
                     flagValue = getattr(jetOrigin[jet._index], originFlag)
                     if (flagValue > 0.5):
                        flavorFlag = 1.
			indexFlag[label] = dictTruth[k]
                     self.out.fillBranch(self.outputName+"_truth_"+originFlag+"_flag",flavorFlag)
	    ### selecting best LLP jet  per event. You loop over the best jet per label list and you peack the best one.  
            taggerResult = getattr(jet, self.taggerName)[label]
	    if taggerResult['output'] > bestResult :
		bestResult = taggerResult['output']
		bestLabel = dict[label]
		bestLabelTruth = indexFlag[label]
		deltaRBestLabel = deltaR(looseLeptons[0],jet)
		 
            closestLepton = None
            minDeltaR = 100.
            for looseLepton in looseLeptons:
                dR = deltaR(looseLepton,jet)
                if dR<minDeltaR:
                    minDeltaR = dR
                    closestLepton = looseLepton
            if closestLepton==None:
                self.out.fillBranch(self.outputName+"_"+label+"_lepton_deltaR",-1)
                self.out.fillBranch(self.outputName+"_"+label+"_lepton_dxy_sig",-1)
            else:
                self.out.fillBranch(self.outputName+"_"+label+"_lepton_deltaR",minDeltaR)
                if math.fabs(closestLepton.dxyErr)<1e-6:
                    self.out.fillBranch(self.outputName+"_"+label+"_lepton_dxy_sig",-1)
                else:
                    self.out.fillBranch(self.outputName+"_"+label+"_lepton_dxy_sig",math.fabs(closestLepton.dxy)/math.fabs(closestLepton.dxyErr))
            self.out.fillBranch(self.outputName+"_"+label+"_jet_pt",jet.pt)
            self.out.fillBranch(self.outputName+"_"+label+"_jet_output",taggerResult['output'])
            self.out.fillBranch(self.outputName+"_"+label+"_jet_parameter",taggerResult['parameter'])
	## Best label will give you the index in the dictionary. 
	self.out.fillBranch(self.outputName+"_TaggerBestJet_deltaR", deltaRBestLabel )
        self.out.fillBranch(self.outputName+"_TaggerBestOutputValue_truth", bestLabelTruth ) 
 	self.out.fillBranch(self.outputName+"_TaggerBestOutputValue" , bestResult )
	self.out.fillBranch(self.outputName+"_TaggerBestOutputLabel", bestLabel )		 
        self.out.fillBranch(self.outputName+"_muonmuon1jet", muonmuon1jet) 
	self.out.fillBranch(self.outputName+"_muonmuonNjets", muonmuonNjets)
	self.out.fillBranch(self.outputName+"_electronelectron1jet", electronelectron1jet)
	self.out.fillBranch(self.outputName+"_electronelectronNjets", electronelectronNjets)
	self.out.fillBranch(self.outputName+"_muonelectron1jet", muonelectron1jet)
	self.out.fillBranch(self.outputName+"_muonelectronNjets", muonelectronNjets)
	self.out.fillBranch(self.outputName+"_electronmuon1jet", electronmuon1jet)
	self.out.fillBranch(self.outputName+"_electronmuonNjets", electronmuonNjets)
	return True 


    
            


          
