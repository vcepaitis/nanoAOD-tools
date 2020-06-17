import os
import sys
import math
import json
import ROOT
import random

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module

class EventCategorization(Module):
    def __init__(
        self,
        globalOptions={"isData":False}, 
        outputName=None ,  
        muonsTight = None , 
        electronsTight = None , 
        muonsLoose = None , 
        electronsLoose = None , 
        jets = None 
    ):
        self.globalOptions=globalOptions
        self.outputName=outputName
        self.muonsTight = muonsTight  
	self.electronsTight = electronsTight
	self.muonsLoose =  muonsLoose
	self.electronsLoose =  electronsLoose
	self.jets = jets
        
    def beginJob(self):
        pass
        
    def endJob(self):
        pass
        
    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        if self.outputName is not None:
            self.out.branch(self.outputName, "I")
            ## declare branches here.
            
    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass
        
    def analyze(self, event):
	##self.out.fillBranch(self.outputName, self.passFilters(event))
	## you need to add the variables you want to store and to fill the branches afterword.
	jets = self.jets(event) 
        tightMuon = self.muonsTight(event) 
        tightElectron = self.electronsTight(event)
	looseMuons = self.muonsLoose(event)
	looseElectrons = self.electronsLoose(event)
	for jet in jets : 
		print "jet tagger output : " , jet.llpdnnx
	if len(tightMuon)  == 1 and len(looseMuons) == 1 and  len(jets) == 1 :
		muonmuon1jet = 1
	elif len(tightMuon) == 1 and len(looseMuons) == 1  and len(jets) > 1 :
		muonmuonNjets = 1
	elif len(tightElectron) == 1 and len(looseElectrons) == 1 and len(jets) == 1 : 
		electronelectron1jet = 1 
	elif len(tightElectron) == 1 and len(looseElectrons) == 1 and  len(jets) > 1 : 
		electronelectronNjets = 1 
	elif len(tightMuon) == 1 and (len(tightElectron) == 1 or len(looseElectrons) == 1) and len(jets) == 1 : 
		if tightMuon[0].pt > tightElectron[0].pt :
			muonelectron1jet = 1 
		else : 
			electronmuon1jet = 1

	## you can't write else because you have the 0 jet catogories that we do not want to have it.  
	elif len(tightMuon) == 1 and ( len(tightElectron) == 1 or len(looseElectrons) == 1 ) and len(jets) > 1 : 
		if tightMuon[0].pt > tightElectron[0].pt :
			muonelectronNjets = 1
		else: 
			electronmuonNjets = 1 
		
	for muon in looseMuons : 
		print "loose muon pt is :  ", muon.pt 
		 
            


    
            


          
