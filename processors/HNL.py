import os
import sys
import math
import json
import argparse
import random
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
from importlib import import_module

from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import PostProcessor
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection, Object
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.modules import *


#parser = argparse.ArgumentParser()
#parser.add_argument('--isData', dest='isData', action='store_true',default=False)
#parser.add_argument('--year', dest='year', action='store',type=int, default=2016)
#parser.add_argument('--input', dest='inputFiles', action='append',default=[])
#parser.add_argument('output', nargs=1)

#args = parser.parse_args()
#print "isData:",args.isData
#print "year:",args.year
#print "inputs:",len(args.inputFiles)

isData = False
year=2016
input_filename = sys.argv[1]
year  = sys.argv[2]
isData = sys.argv[3]
outpout = sys.argv[5]
os.environ["X509_USER_PROXY"] = sys.argv[4]
outpout = sys.argv[5]
print os.environ["X509_USER_PROXY"]

#output 


#for inputFile in args.inputFiles:
print ("The name of the file is ", input_filename)

input_file = open(input_filename, 'r')
for File in input_file:

   FilesVector = []
   File = File.rstrip()
   print "file name is : ", File  , "\n"
   FilesVector.append(File)
   rootFile = ROOT.TFile.Open(File)
   if not rootFile:
	print "yes same issue",rootFile
        print "CRITICAL - file '"+str(input_file)+"' not found!"
        sys.exit(1)
   tree = rootFile.Get("Events")
   if not tree:
     #  print "CRITICAL - 'Events' tree not found in file '"+inputFile+"'!"
        sys.exit(1)
   #print " - ",inputFile,", events=",tree.GetEntries()
   #print "output directory:",args.output[0]

   globalOptions = {
      "isData":isData,
      "year":2016
   }


   muonSelection = [
      MuonSelection(
        outputName="tightMuons",
        storeKinematics=['pt','eta', 'dxy', 'dxyErr', 'dz', 'dzErr', 'phi'],
        storeWeights=True,
        muonMinPt = 25.,
        muonMaxDxy = 0.002,
        muonMaxDz = 0.01,
        muonID = MuonSelection.TIGHT,
        muonIso = MuonSelection.TIGHT,
        globalOptions=globalOptions
     ),
     MuonSelection(
        inputCollection = lambda event: event.tightMuons_unselected,
        outputName="looseMuons",
        storeKinematics=['pt','eta', 'dxy', 'dxyErr', 'dz', 'dzErr', 'phi'],
        storeWeights=False,
        muonMinPt = 5.,
        muonID = MuonSelection.LOOSE,
        muonIso = MuonSelection.NONE,
        globalOptions=globalOptions
     ),
    
     SingleMuonTriggerSelection(
        inputCollection=lambda event: event["tightMuons"],
        outputName="IsoMuTrigger",
        storeWeights=False,
        globalOptions=globalOptions
     ),
     EventSkim(selection=lambda event: event.IsoMuTrigger_flag==1),
     EventSkim(selection=lambda event: event.ntightMuons>0),
     EventSkim(selection=lambda event: event.nlooseMuons>0)
   ]

   analyzerChain = []

   analyzerChain.extend(muonSelection)

   analyzerChain.append(
      JetSelection(
      )
   )


   storeVariables = [
      [lambda tree: tree.branch("genweight","F"),lambda tree,event: tree.fillBranch("genweight",event.Generator_weight)],
   ]


   analyzerChain.append(EventInfo(storeVariables=storeVariables))

   analyzerChain.append(EventSkim(selection=lambda event: len(event.selectedJets)>0))


   analyzerChain.append(
     InvariantSystem(
        inputCollection = lambda event: [event.looseMuons[0], event.tightMuons[0]],
        outputName = "dimuon"
     )
   )

   analyzerChain.append(
     LepJetFinder(
        jetCollection = lambda event: event.selectedJets,
        leptonCollection = lambda event: event.looseMuons,
     )
   )
 
   analyzerChain.append(
     InvariantSystem(
        inputCollection = lambda event: [event.tightMuons[0], event.lepJet[0]],
        outputName = "lepjet_muon"
     )
   )

   analyzerChain.append( 
     MetFilter(
       outputName ="MET_filter"
    )
   )

   analyzerChain.append(
     EventObservables(
       jetCollection = lambda event: event.selectedJets,
       leptonCollection = lambda event: event.tightMuons[0] 
     )


   ) 

   analyzerChain.append(
     TaggerEvaluation(
        modelPath="PhysicsTools/NanoAODTools/data/nn/weight2016_75.pb",
        logctauValues = [1.74],
        inputCollections=[
            lambda event: event.lepJet
        ],
        taggerName="llpdnnx",
     )
   )


   analyzerChain.append(
     JetTaggerResult(
        inputCollection = lambda event: event.lepJet,
        taggerName = "llpdnnx",
        outputName = "lepJet",
        logctauValues = [1.74],
        predictionLabels = ["LLP_Q", "LLP_QMU"],
     )
   )


   analyzerChain.append(
     JetTruthFlags(inputCollection= lambda event: event.selectedJets,
        outputName="selectedJets"
     )
   )

   analyzerChain.append(
     JetTruthFlags(inputCollection= lambda event: event.lepJet,
        outputName="lepJet"
     )
   )


   p=PostProcessor(
    sys.argv[5],
    [FilesVector],
    modules=analyzerChain,
    maxEvents=-1,
    friend=True
   )
   p.run()
