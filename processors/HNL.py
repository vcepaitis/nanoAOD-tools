import os
import sys
import math
import argparse
import random
import ROOT
import numpy as np

from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor \
    import PostProcessor
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel \
    import Collection, Object
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.modules import *

parser = argparse.ArgumentParser()

parser.add_argument('--testMode', dest='testMode', action='store_true', default=False)
parser.add_argument('--isData', dest='isData',
                    action='store_true', default=False)
parser.add_argument('--isSignal', dest='isSignal',
                    action='store_true', default=False)
parser.add_argument('--year', dest='year',
                    action='store', type=int, default=2016)
parser.add_argument('--input', dest='inputFiles', action='append', default=[])
parser.add_argument('--noTagger', dest='noTagger', action='store_true', default=False)
parser.add_argument('output', nargs=1)

args = parser.parse_args()

testMode = args.testMode
print "isData:",args.isData
print "isSignal:",args.isSignal
print "inputs:",len(args.inputFiles)
print "Running tagger", not args.noTagger

for inputFile in args.inputFiles:
    if "-2016" in inputFile or "Run2016" in inputFile:
        year = 2016
    elif "-2017" in inputFile or "Run2017" in inputFile:
        year = 2017
    elif "-2018" in inputFile or "Run2018" in inputFile:
        year = 2018
    else:
        year = args.year
    rootFile = ROOT.TFile.Open(inputFile)
    if not rootFile:
        print "CRITICAL - file '"+inputFile+"' not found!"
        sys.exit(1)
    tree = rootFile.Get("Events")
    if not tree:
        print "CRITICAL - 'Events' tree not found in file '"+inputFile+"'!"
        sys.exit(1)
    print " - ", inputFile, ", events=", tree.GetEntries()

print "year:", year
print "output directory:", args.output[0]

globalOptions = {
    "isData": args.isData,
    "isSignal": args.isSignal,
    "year": year
}

isMC = not args.isData

minMuonPt = {2016: 25., 2017: 28., 2018: 25.}
minElectronPt = {2016: 29., 2017: 34., 2018: 34.}

if isMC:
    jecTags = {2016: 'Summer16_07Aug2017_V11_MC',
               2017: 'Fall17_17Nov2017_V32_MC',
               2018: 'Autumn18_V19_MC'
               }

    jerTags = {2016: 'Summer16_25nsV1_MC',
               2017: 'Fall17_V3_MC',
               2018: 'Autumn18_V7_MC'
               }

if args.isData:
    jecTags = {2016: 'Summer16_07Aug2017All_V11_DATA',
               2017: 'Fall17_17Nov2017_V32_DATA',
               2018: 'Autumn18_V19_DATA'
               }

met_variable = {2016: lambda event: Object(event, "MET"),
                2017: lambda event: Object(event, "METFixEE2017"),
                2018: lambda event: Object(event, "MET")
                }


leptonSelection = [
    EventSkim(selection=lambda event: event.nTrigObj > 0),
    MuonSelection(
        outputName="tightMuon",
        storeKinematics=['pt', 'eta', 'dxy', 'dxyErr', 'dz',
                         'dzErr', 'phi', 'pfRelIso04_all', 'charge'],
        storeWeights=True,
        muonMinPt=minMuonPt[globalOptions["year"]],
        triggerMatch=True,
        muonID=MuonSelection.TIGHT,
        muonIso=MuonSelection.TIGHT,
        selectLeadingOnly=True,
        globalOptions=globalOptions
    ),
    ElectronSelection(
        outputName="tightElectron",
        storeKinematics=['pt', 'eta', 'dxy', 'dxyErr', 'dz',
                         'dzErr', 'phi','pfRelIso03_all', 'charge'],
        electronMinPt=minElectronPt[globalOptions["year"]],
        electronID="Iso_WP90",
        storeWeights=True,
        triggerMatch=True,
        selectLeadingOnly=True,
        globalOptions=globalOptions
    ),
    EventSkim(selection=lambda event: event.ntightMuon + event.ntightElectron > 0),
    SingleMuonTriggerSelection(
        inputCollection=lambda event: event.tightMuon,
        outputName="IsoMuTrigger",
        storeWeights=True,
        globalOptions=globalOptions
    ),
    SingleElectronTriggerSelection(
        inputCollection=lambda event: event.tightElectron,
        outputName="IsoElectronTrigger",
        storeWeights=False,
        globalOptions=globalOptions
    ),
    EventSkim(selection=lambda event: event.IsoMuTrigger_flag + event.IsoElectronTrigger_flag > 0),
    MuonSelection(
        inputCollection=lambda event: event.tightMuon_unselected,
        outputName="looseMuons",
        storeKinematics=['pt', 'eta', 'dxy', 'dxyErr', 'dz',
                         'dzErr', 'phi', 'pfRelIso04_all',
                         'looseId', 'mediumId', 'tightId', 'charge'],
        muonMinPt=5.,
        muonID=MuonSelection.LOOSE,
        muonIso=MuonSelection.NONE,
        globalOptions=globalOptions
    ),
    ElectronSelection(
        inputCollection=lambda event: event.tightElectron_unselected,
        outputName="looseElectrons",
        storeKinematics=['pt', 'eta', 'dxy', 'dxyErr', 'dz',
                         'dzErr', 'phi','pfRelIso03_all', 'charge',
                         'mvaFall17V2noIso_WP80', 'mvaFall17V2noIso_WP90', 'mvaFall17V2noIso_WPL',
                         'mvaFall17V2Iso_WP80', 'mvaFall17V2Iso_WP90', 'mvaFall17V2Iso_WPL',
                         'cutBased'],

        electronMinPt=5.,
        electronID="Custom",
        globalOptions=globalOptions
    ),


    LeptonCollecting(
        tightMuonCollection=lambda event:event.tightMuon,
        tightElectronCollection=lambda event:event.tightElectron,
        looseMuonCollection=lambda event:event.looseMuons,
        looseElectronCollection=lambda event:event.looseElectrons,
        outputName = "Leptons"
    )

]

analyzerChain = []

analyzerChain.extend(leptonSelection)


analyzerChain.append(
    InvariantSystem(
        inputCollection= lambda event:
            sorted(event.tightMuon+event.looseMuons+event.tightElectron+event.looseElectrons,key=lambda x: -x.pt)[:2],
        outputName="dilepton"
    )
)

'''
# left for debugging
from PhysicsTools.NanoAODTools.postprocessing.modules.jme.jetmetHelperRun2       import * 
jetmetCorrector = createJMECorrector(isMC=isMC, dataYear=year, runPeriod="B", jesUncert="All", redojec=True)
#jetmetCorrector = createJMECorrector(isMC=False, dataYear=2017, runPeriod="E", metBranchName="METFixEE2017")  
analyzerChain.append(jetmetCorrector())
'''


#featureDictFile = "${CMSSW_BASE}/src/PhysicsTools/NanoAODTools/data/nn/200311/feature_dict.py"
featureDictFile = "${CMSSW_BASE}/src/PhysicsTools/NanoAODTools/data/nn/200720/feature_dict.py"
modelPath = {
    2016: "${CMSSW_BASE}/src/PhysicsTools/NanoAODTools/data/nn/200720/weight2016.pb",
    2017: "${CMSSW_BASE}/src/PhysicsTools/NanoAODTools/data/nn/200720/weight2017.pb",
    2018: "${CMSSW_BASE}/src/PhysicsTools/NanoAODTools/data/nn/200720/weight2018.pb"
}
modelGunPath = {
    2016: "${CMSSW_BASE}/src/PhysicsTools/NanoAODTools/data/nn/200720/weightGun2016.pb",
    2017: "${CMSSW_BASE}/src/PhysicsTools/NanoAODTools/data/nn/200720/weightGun2017.pb",
    2018: "${CMSSW_BASE}/src/PhysicsTools/NanoAODTools/data/nn/200720/weightGun2018.pb"
}

jesUncertaintyFile = {
    2016: "${CMSSW_BASE}/src/PhysicsTools/NanoAODTools/data/jme/Summer16_07Aug2017_V11_MC_Uncertainty_AK4PFchs.txt",
    2017: "${CMSSW_BASE}/src/PhysicsTools/NanoAODTools/data/jme/Fall17_17Nov2017_V32_MC_Uncertainty_AK4PFchs.txt",
    2018: "${CMSSW_BASE}/src/PhysicsTools/NanoAODTools/data/jme/Autumn18_V19_MC_Uncertainty_AK4PFchs.txt"
}
jerResolutionFile = {
    2016: "${CMSSW_BASE}/src/PhysicsTools/NanoAODTools/data/jme/Summer16_25nsV1_MC_PtResolution_AK4PFchs.txt",
    2017: "${CMSSW_BASE}/src/PhysicsTools/NanoAODTools/data/jme/Fall17_V3_MC_PtResolution_AK4PFchs.txt",
    2018: "${CMSSW_BASE}/src/PhysicsTools/NanoAODTools/data/jme/Autumn18_V7_MC_PtResolution_AK4PFchs.txt"
}

jerSFUncertaintyFile = {
    2016: "${CMSSW_BASE}/src/PhysicsTools/NanoAODTools/data/jme/Summer16_25nsV1_MC_SF_AK4PFchs.txt",
    2017: "${CMSSW_BASE}/src/PhysicsTools/NanoAODTools/data/jme/Fall17_V3_MC_SF_AK4PFchs.txt",
    2018: "${CMSSW_BASE}/src/PhysicsTools/NanoAODTools/data/jme/Autumn18_V7_MC_SF_AK4PFchs.txt"
}


analyzerChain.append(
     MetFilter(
        globalOptions=globalOptions,
        outputName="MET_filter"
     )
)


if isMC:
    analyzerChain.append(
        JetMetUncertainties(
            metInput=met_variable[year],
            rhoInput = lambda event: event.fixedGridRhoFastjetAll,
            jetCollection = lambda event: Collection(event,"Jet"),
            lowPtJetCollection = lambda event: Collection(event,"CorrT1METJet"),
            genJetCollection = lambda event: Collection(event,"GenJet"),
            muonCollection = lambda event: Collection(event,"Muon"),
            electronCollection = lambda event: Collection(event,"Electron"),
            jesUncertaintyFile=jesUncertaintyFile[year],
            jerResolutionFileName=jerResolutionFile[year],
            jerSFUncertaintyFileName=jerSFUncertaintyFile[year],
            propagateJER = False,
            jetKeys = ['pt', 'eta', 'phi' , 'jetId', 'nConstituents'],
        )
    )

    for systName, jetCollection in [
        ("nominal", lambda event: event.jets_nominal),
        ("jerUp", lambda event: event.jets_jerUp),
        ("jerDown", lambda event: event.jets_jerDown),
        ("jesTotalUp", lambda event: event.jets_jesTotalUp),
        ("jesTotalDown", lambda event: event.jets_jesTotalDown),
    ]:

        analyzerChain.append(
            JetSelection(
                inputCollection=jetCollection,
                leptonCollectionDRCleaning=lambda event: event.leadingLeptons,
                leptonCollectionP4Subraction=lambda event: event.subleadingLeptons,
                jetMinPt=15.,
                jetMaxEta=2.399, #TODO: change to 2.4
                jetMinNConstituents=3,
                jetId=JetSelection.LOOSE,
                outputName="selectedJets_"+systName,
                globalOptions=globalOptions
            )
        )

        analyzerChain.append(
            JetSelection(
                inputCollection=jetCollection,
                leptonCollectionDRCleaning=lambda event: event.leadingLeptons,
                jetMinPt=30.,
                jetMinEta=2.4,
                jetMaxEta=5.,
                jetId=JetSelection.LOOSE,
                storeKinematics=[],
                outputName="selectedFwdJets_"+systName,
                globalOptions=globalOptions
            )
        )

        analyzerChain.append(
            JetSelection(
                inputCollection=jetCollection,
                jetMinPt=100.,
                jetMinEta=2.25,
                jetMaxEta=3.0,
                jetId=JetSelection.LOOSE,
                storeKinematics=[],
                outputName="selectedL1PreFiringJets_"+systName,
                globalOptions=globalOptions
            )
        )

        analyzerChain.append(
            JetTruthFlags(
                inputCollection=lambda event, systName=systName: getattr(event, "selectedJets_"+systName),
                outputName="selectedJets_"+systName,
                globalOptions=globalOptions
            )
        )

        analyzerChain.append(
            LepJetFinder(
                jetCollection=lambda event, systName=systName: getattr(event, "selectedJets_"+systName),
                leptonCollection=lambda event: event.subleadingLeptons,
                outputName="lepJet_"+systName
            )
        )

    analyzerChain.append(
        EventSkim(selection=lambda event: \
            getattr(event, "nselectedJets_nominal") > 0 or
            getattr(event, "nselectedJets_jesTotalUp") > 0 or
            getattr(event, "nselectedJets_jesTotalDown") > 0 or
            getattr(event, "nselectedJets_jerUp") > 0 or
            getattr(event, "nselectedJets_jerDown") > 0
        )
    )

    
    analyzerChain.append(
        TaggerEvaluationProfiled(
            modelPath=modelPath[year],
            featureDictFile=featureDictFile,
            inputCollections=[
                lambda event: event.selectedJets_nominal[:4],
                lambda event: event.selectedJets_jesTotalUp[:4],
                lambda event: event.selectedJets_jesTotalDown[:4],
                lambda event: event.selectedJets_jerUp[:4],
                lambda event: event.selectedJets_jerDown[:4]
            ],
            taggerName="llpdnnx",
            globalOptions=globalOptions,
            evalValues = np.linspace(-3,2,5*5+1),
        )
    )
   
    for systName, jetCollection in [
        ("nominal", lambda event: event.selectedJets_nominal[:4]),
        ("jerUp", lambda event: event.selectedJets_jerUp[:4]),
        ("jerDown", lambda event: event.selectedJets_jerDown[:4]),
        ("jesTotalUp", lambda event: event.selectedJets_jesTotalUp[:4]),
        ("jesTotalDown", lambda event: event.selectedJets_jesTotalDown[:4]),
    ]:

        analyzerChain.append(
           EventCategorization(	
                looseLeptons=lambda event: event.subleadingLeptons,
                jetsCollection=jetCollection,
                taggerName="llpdnnx",
                outputName="category_"+systName,
                globalOptions=globalOptions
           )
        )
    

    for systName, jetCollection, metObject in [
        ("nominal", lambda event: event.selectedJets_nominal,
            lambda event: event.met_nominal),
        ("jerUp", lambda event: event.selectedJets_jerUp,
           lambda event: event.met_jerUp),
        ("jerDown", lambda event: event.selectedJets_jerDown,
            lambda event: event.met_jerDown),
        ("jesTotalUp", lambda event: event.selectedJets_jesTotalUp,
            lambda event: event.met_jesTotalUp),
        ("jesTotalDown", lambda event: event.selectedJets_jesTotalDown,
            lambda event: event.met_jesTotalDown),
        ("unclEnUp", lambda event: event.selectedJets_nominal,
            lambda event: event.met_unclEnUp),
        ("unclEnDown", lambda event: event.selectedJets_nominal,
            lambda event: event.met_unclEnDown),
    ]:
        analyzerChain.extend([
            WbosonReconstruction(
                leptonCollectionName='leadingLeptons',
                metObject=metObject,
                globalOptions=globalOptions,
                outputName=systName
            )
        ])

        analyzerChain.append(
            EventObservables(
                jetCollection=jetCollection,
                leptonCollection=lambda event: event.leadingLeptons[0],
                metInput=metObject,
                globalOptions=globalOptions,
                outputName="EventObservables_"+systName
            )
        )

        
        analyzerChain.append(
            XGBEvaluation(
                systName=systName,
                jetCollection=jetCollection
            )
        )
        

else:
    analyzerChain.append(
        JetSelection(
            inputCollection=lambda event: Collection(event, "Jet"),
            leptonCollectionDRCleaning=lambda event: event.leadingLeptons,
            leptonCollectionP4Subraction=lambda event: event.subleadingLeptons,
            jetMinPt=15.,
            jetMaxEta=2.399, #TODO: change to 2.4
            jetMinNConstituents=3,
            jetId=JetSelection.LOOSE,
            outputName="selectedJets_nominal",
            globalOptions=globalOptions
        )
    )

    analyzerChain.append(
        JetSelection(
            inputCollection=lambda event: Collection(event, "Jet"),
            leptonCollectionDRCleaning=lambda event: event.leadingLeptons,
            jetMinPt=30.,
            jetMinEta=2.4,
            jetMaxEta=5.,
            jetId=JetSelection.LOOSE,
            storeKinematics=[],
            outputName="selectedFwdJets_nominal",
            globalOptions=globalOptions
        )
    )

    analyzerChain.append(
        JetSelection(
            inputCollection=lambda event: Collection(event, "Jet"),
            jetMinPt=100.,
            jetMinEta=2.25,
            jetMaxEta=3.0,
            jetId=JetSelection.LOOSE,
            storeKinematics=[],
            outputName="selectedL1PreFiringJets_nominal",
            globalOptions=globalOptions
        )
    )

    analyzerChain.append(
        LepJetFinder(
            jetCollection=lambda event: event.selectedJets_nominal,
            leptonCollection=lambda event: event.subleadingLeptons,
            outputName="lepJet_nominal"
        )
    )

    analyzerChain.append(
        EventSkim(
            selection=lambda event: event.nselectedJets_nominal > 0
        )
    )

    
    analyzerChain.append(
        TaggerEvaluationProfiled(
            modelPath=modelPath[year],
            featureDictFile=featureDictFile,
            inputCollections=[
                lambda event: event.selectedJets_nominal[:4]
            ],
            taggerName="llpdnnx",
            globalOptions=globalOptions,
            evalValues = np.linspace(-3,2,5*5+1)
        )
    )
    
    
    analyzerChain.append(
	EventCategorization(
            looseLeptons=lambda event: event.subleadingLeptons,
            jetsCollection=lambda event: event.selectedJets_nominal[:4],
            taggerName="llpdnnx",
            outputName="category_nominal",
            globalOptions=globalOptions
       
       )
    )
    analyzerChain.extend([
        WbosonReconstruction(
            leptonCollectionName='leadingLeptons',
            metObject=met_variable[year],
            globalOptions=globalOptions,
            outputName="nominal"
        )
    ])

    analyzerChain.append(
        EventObservables(
            jetCollection=lambda event: event.selectedJets_nominal,
            leptonCollection=lambda event: event.leadingLeptons[0],
            metInput=met_variable[year],
            globalOptions=globalOptions,
            outputName="EventObservables_nominal"
        )
    )

    analyzerChain.append(
        XGBEvaluation(
            systName="nominal",
            jetCollection=lambda event: event.selectedJets_nominal
        )
    )

    analyzerChain.append(
        PileupWeight(
            outputName="puweight",
            globalOptions=globalOptions
        )
    )

storeVariables = [
    [lambda tree: tree.branch("PV_npvs", "I"), lambda tree,
     event: tree.fillBranch("PV_npvs", event.PV_npvs)],
    [lambda tree: tree.branch("PV_npvsGood", "I"), lambda tree,
     event: tree.fillBranch("PV_npvsGood", event.PV_npvsGood)],
    [lambda tree: tree.branch("fixedGridRhoFastjetAll", "F"), lambda tree,
     event: tree.fillBranch("fixedGridRhoFastjetAll",
                            event.fixedGridRhoFastjetAll)],
]


if not globalOptions["isData"]:
    storeVariables.append([lambda tree: tree.branch("genweight", "F"),
                           lambda tree,
                           event: tree.fillBranch("genweight",
                           event.Generator_weight)])

    if args.isSignal:
        for coupling in range(1,68):
            storeVariables.append([
                lambda tree, coupling=coupling: tree.branch('LHEWeights_coupling_%i'%coupling,'F'),
                lambda tree, event, coupling=coupling: tree.fillBranch('LHEWeights_coupling_%i'%coupling,getattr(event,"LHEWeights_coupling_%i"%coupling)),
            ])

analyzerChain.append(EventInfo(storeVariables=storeVariables))
taggerTypes = ['EventCategorization', 'TaggerEvaluationProfiled']

if testMode:
    for ianalyzer, analyzer in enumerate(analyzerChain):
        if type(analyzer).__name__ == "PileupWeight":
            analyzerChain.pop(ianalyzer)   

if args.noTagger:
    analyzerChainNew = analyzerChain
    for ianalyzer, analyzer in enumerate(analyzerChain):
        if type(analyzer).__name__ not in taggerTypes:
            analyzerChainNew.append(analyzer)
    analyzerChain = analyzerChainNew

p = PostProcessor(
    args.output[0],
    [args.inputFiles],
    cut="(nJet>0)&&((nElectron+nMuon)>0)",
    modules=analyzerChain,
    maxEvents=-1,
    friend=True
)

p.run()
