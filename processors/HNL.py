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
parser.add_argument('--year', dest='year',
                    action='store', type=int, default=2016)
parser.add_argument('--input', dest='inputFiles', action='append', default=[])
parser.add_argument('--noTagger', dest='noTagger', action='store_true', default=False)
parser.add_argument('--profile', action='store_true', default=False)
parser.add_argument('--bdt', action='store_true', default=False)
parser.add_argument('output', nargs=1)

args = parser.parse_args()

testMode = args.testMode
isSignal = False

print "isData:",args.isData
print "inputs:",len(args.inputFiles)
print "Running tagger", not args.noTagger


for inputFile in args.inputFiles:
    if "dirac" in inputFile or "majorana" in inputFile or "LLPGun" in inputFile: 
        isSignal = True
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
print "isSignal:",isSignal

globalOptions = {
    "isData": args.isData,
    "isSignal": isSignal,
    "year": year
}

isMC = not args.isData

minMuonPt = {2016: 26., 2017: 29., 2018: 26.}
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

muonIso = MuonSelection.NONE if args.bdt else MuonSelection.TIGHT
eleId = "noIso_WP90" if args.bdt else "Iso_WP90"
leptonSelection = [
    EventSkim(selection=lambda event: event.nTrigObj > 0),
    MuonSelection(
        outputName="tightMuons",
        storeKinematics=['pt', 'eta', 'dxy', 'dxyErr', 'dz',
                         'dzErr', 'phi', 'pfRelIso04_all', 'charge'],
        storeWeights=True,
        muonMinPt=minMuonPt[globalOptions["year"]],
        muonMaxDxy=0.01,
        muonMaxDz=0.05,
        muonID=MuonSelection.TIGHT,
        muonIso=muonIso,
        globalOptions=globalOptions
    ),
    ElectronSelection(
        outputName="tightElectrons",
        storeKinematics=['pt', 'eta', 'dxy', 'dxyErr', 'dz',
                         'dzErr', 'phi','pfRelIso03_all', 'charge'],
        electronMinPt=minElectronPt[globalOptions["year"]],
        electronID=eleId,
        storeWeights=True,
        electronIPCuts=True,
        globalOptions=globalOptions
    ),
    EventSkim(selection=lambda event: event.ntightMuons + event.ntightElectrons > 0),
    SingleMuonTriggerSelection(
        inputCollection=lambda event: event.tightMuons,
        outputName="IsoMuTrigger",
        storeWeights=True,
        globalOptions=globalOptions
    ),
    SingleElectronTriggerSelection(
        inputCollection=lambda event: event.tightElectrons,
        outputName="IsoElectronTrigger",
        storeWeights=False,
        globalOptions=globalOptions
    ),
    MuonSelection(
        inputCollection=lambda event: event.tightMuons_unselected,
        outputName="looseMuons",
        storeKinematics=['pt', 'eta', 'dxy', 'dxyErr', 'dz',
                         'dzErr', 'phi', 'pfRelIso04_all',
                         'looseId', 'mediumId', 'tightId', 'charge'],
        muonMinPt=3.,
        muonID=MuonSelection.LOOSE,
        muonIso=MuonSelection.NONE,
        globalOptions=globalOptions
    ),
    MuonSelection(
        inputCollection=lambda event: event.looseMuons,
        outputName="looseIsoMuons",
        storeKinematics=['pt', 'eta', 'dxy', 'dxyErr', 'dz',
                         'dzErr', 'phi', 'pfRelIso04_all',
                         'looseId', 'mediumId', 'tightId', 'charge'],
        muonMinPt=3.,
        storeWeights=False,
        muonID=MuonSelection.LOOSE,
        muonIso=MuonSelection.TIGHT,
        globalOptions=globalOptions
    ),
    ElectronSelection(
        inputCollection=lambda event: event.tightElectrons_unselected,
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
    ElectronSelection(
        inputCollection=lambda event: event.looseElectrons,
        outputName="looseIsoElectrons",
        storeKinematics=['pt', 'eta', 'dxy', 'dxyErr', 'dz',
                         'dzErr', 'phi','pfRelIso03_all', 'charge',
                         'mvaFall17V2noIso_WP80', 'mvaFall17V2noIso_WP90', 'mvaFall17V2noIso_WPL',
                         'mvaFall17V2Iso_WP80', 'mvaFall17V2Iso_WP90', 'mvaFall17V2Iso_WPL',
                         'cutBased'],
        electronMinPt=5.,
        storeWeights=False,
        electronID="CustomIso",
        globalOptions=globalOptions
    ),
    LeptonCollecting(
        tightMuonsCollection=lambda event:event.tightMuons,
        tightElectronsCollection=lambda event:event.tightElectrons,
        looseMuonCollection=lambda event:event.looseMuons,
        looseElectronCollection=lambda event:event.looseElectrons,
        outputName = "Leptons"
    ),
    #EventSkim(selection=lambda event: event.nsubleadingLeptons > 0)
]

analyzerChain = []

analyzerChain.extend(leptonSelection)

if not args.bdt:
    analyzerChain.append(EventSkim(selection=lambda event: event.IsoMuTrigger_flag + event.IsoElectronTrigger_flag > 0))
    analyzerChain.append(EventSkim(selection=lambda event: event.isTriggered))


analyzerChain.append(
    InvariantSystem(
        inputCollection= lambda event:
            sorted(event.tightMuons+event.looseMuons+event.tightElectrons+event.looseElectrons,key=lambda x: -x.pt)[:2],
        outputName="dilepton"
    )
)

featureDictFile = "${CMSSW_BASE}/src/PhysicsTools/NanoAODTools/data/nn/201117/feature_dict.py"
modelPath = {
    2016: "${CMSSW_BASE}/src/PhysicsTools/NanoAODTools/data/nn/201117/weightMixed2016_ExtNominalNetwork_origSV_DA_20_lr001_201117.pb",
    2017: "${CMSSW_BASE}/src/PhysicsTools/NanoAODTools/data/nn/201117/weightMixed2017_ExtNominalNetwork_origSV_DA_20_lr001_201117.pb",
    2018: "${CMSSW_BASE}/src/PhysicsTools/NanoAODTools/data/nn/201117/weightMixed2018_ExtNominalNetwork_origSV_DA_20_lr001_201117.pb"
}

BDTmodelPath = {
    2016: "${CMSSW_BASE}/src/PhysicsTools/NanoAODTools/data/bdt/201117/nominal/bdt_2016.model",
    2017: "${CMSSW_BASE}/src/PhysicsTools/NanoAODTools/data/bdt/201117/nominal/bdt_2017.model",
    2018: "${CMSSW_BASE}/src/PhysicsTools/NanoAODTools/data/bdt/201117/nominal/bdt_2018.model"
}

BDTmodelPathUncorr = {
    2016: "${CMSSW_BASE}/src/PhysicsTools/NanoAODTools/data/bdt/201117/uncorrelated/bdt_2016.model",
    2017: "${CMSSW_BASE}/src/PhysicsTools/NanoAODTools/data/bdt/201117/uncorrelated/bdt_2017.model",
    2018: "${CMSSW_BASE}/src/PhysicsTools/NanoAODTools/data/bdt/201117/uncorrelated/bdt_2018.model"
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
            jetKeys = ['pt', 'eta', 'phi' , 'jetId', 'nConstituents', 'rawFactor'],
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
                leptonCollectionDRCleaning=lambda event: event.tightMuons+event.tightElectrons+event.looseIsoMuons+event.looseIsoElectrons,
                leptonCollectionP4Subraction=lambda event:event.looseMuons+event.looseElectrons,
                jetMinPt=15.,
                jetMaxEta=2.4,
                jetId=JetSelection.TIGHT,
                outputName="selectedJets_"+systName,
                globalOptions=globalOptions
            )
        )

        analyzerChain.append(
            JetSelection(
                inputCollection=jetCollection,
                leptonCollectionDRCleaning=lambda event: event.tightMuons+event.tightElectrons+event.looseIsoMuons+event.looseIsoElectrons,
                jetMinPt=30.,
                jetMinEta=2.4,
                jetMaxEta=5.,
                jetId=JetSelection.TIGHT,
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
                jetId=JetSelection.TIGHT,
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
        EventSkim(selection=lambda event: \
            getattr(event, "nselectedJets_nominal") > 0 or
            getattr(event, "nselectedJets_jesTotalUp") > 0 or
            getattr(event, "nselectedJets_jesTotalDown") > 0 or
            getattr(event, "nselectedJets_jerUp") > 0 or
            getattr(event, "nselectedJets_jerDown") > 0
        )
    )

    analyzerChain.append(
        EventSkim(selection=lambda event: \
            getattr(event, "nselectedJets_nominal") < 5 or
            getattr(event, "nselectedJets_jesTotalUp") < 5 or
            getattr(event, "nselectedJets_jesTotalDown") < 5 or
            getattr(event, "nselectedJets_jerUp") < 5 or
            getattr(event, "nselectedJets_jerDown") < 5
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
                metInput=metObject,
                globalOptions=globalOptions,
                outputName="EventObservables_"+systName
            )
        )

        analyzerChain.append(
           SimplifiedEventCategorization(
                maxDeltaR=1.3,
                looseLeptons=lambda event: event.subleadingLeptons,
                jetsCollection=jetCollection,
                outputName="category_"+systName,
                globalOptions=globalOptions
           )
        )

        analyzerChain.append(
            MassReconstruction(
                globalOptions=globalOptions,
                outputName="category_"+systName,
                tightLeptons=lambda event: event.leadingLeptons,
                looseLeptons=lambda event: event.subleadingLeptons,
                jets=jetCollection,
            )
        )


    analyzerChain.append(
        XGBEvaluation(
            systematics=["nominal", "jerUp", "jerDown", "jesTotalUp", "jesTotalDown", "unclEnUp", "unclEnDown"],
            jetCollections=[
            lambda event: event.selectedJets_nominal,
            lambda event: event.selectedJets_jesTotalUp,
            lambda event: event.selectedJets_jesTotalDown,
            lambda event: event.selectedJets_jerUp,
            lambda event: event.selectedJets_jerDown,
            lambda event: event.selectedJets_nominal,
            lambda event: event.selectedJets_nominal
        ],
            modelPath=BDTmodelPath[year],
            inputFeatures="${CMSSW_BASE}/src/PhysicsTools/NanoAODTools/data/bdt/201117/nominal/bdt_inputs.txt",
            outputName="bdt_score"
        )
    )

    analyzerChain.append(
        XGBEvaluation(
            systematics=["nominal", "jerUp", "jerDown", "jesTotalUp", "jesTotalDown", "unclEnUp", "unclEnDown"],
            jetCollections=[
            lambda event: event.selectedJets_nominal,
            lambda event: event.selectedJets_jesTotalUp,
            lambda event: event.selectedJets_jesTotalDown,
            lambda event: event.selectedJets_jerUp,
            lambda event: event.selectedJets_jerDown,
            lambda event: event.selectedJets_nominal,
            lambda event: event.selectedJets_nominal
        ],
            modelPath=BDTmodelPathUncorr[year],
            inputFeatures="${CMSSW_BASE}/src/PhysicsTools/NanoAODTools/data/bdt/201117/uncorrelated/bdt_inputs.txt",
            outputName="bdt_score_uncorr"
        )
    )

    analyzerChain.append(
        TaggerEvaluationProfiled(
            modelPath=modelPath[year],
            featureDictFile=featureDictFile,
            inputCollections=[
                lambda event: event.category_nominal_jets,
                lambda event: event.category_jesTotalUp_jets,
                lambda event: event.category_jesTotalDown_jets,
                lambda event: event.category_jerUp_jets,
                lambda event: event.category_jerDown_jets
            ],
            taggerName="llpdnnx",
            profiledLabelDict = {
                'LLP_Q': ['LLP_Q','LLP_QTAU_H','LLP_QTAU_3H'],
                'LLP_QE': [ 'LLP_QE'],
                'LLP_QMU': [ 'LLP_QMU']
            },
            globalOptions=globalOptions,
            evalValues = np.linspace(-1.9,1.9,5*4),
        )
    )

    for systName, jetCollection in [
        ("nominal", lambda event: event.category_nominal_jets),
        ("jerUp", lambda event: event.category_jesTotalUp_jets),
        ("jerDown", lambda event: event.category_jesTotalDown_jets),
        ("jesTotalUp", lambda event: event.category_jesTotalUp_jets),
        ("jesTotalDown", lambda event: event.category_jesTotalDown_jets),
    ]:
        analyzerChain.append(
            JetTaggerProfiledResult(
            inputCollection=jetCollection,
            taggerName="llpdnnx",
            outputName="category_"+systName+"_jets",
            profiledLabels = ['LLP_Q','LLP_QE','LLP_QMU'],
            globalOptions=globalOptions
            )
        )

else:
    analyzerChain.append(
        JetSelection(
            inputCollection=lambda event: Collection(event, "Jet"),
            leptonCollectionDRCleaning=lambda event: event.tightMuons+event.tightElectrons+event.looseIsoMuons+event.looseIsoElectrons,
            leptonCollectionP4Subraction=lambda event:event.looseMuons+event.looseElectrons,
            jetMinPt=15.,
            jetMaxEta=2.4,
            jetId=JetSelection.TIGHT,
            outputName="selectedJets_nominal",
            globalOptions=globalOptions
        )
    )


    analyzerChain.append(
        JetSelection(
            inputCollection=lambda event: Collection(event, "Jet"),
            leptonCollectionDRCleaning=lambda event: event.tightMuons+event.tightElectrons+event.looseIsoMuons+event.looseIsoElectrons,
            jetMinPt=30.,
            jetMinEta=2.4,
            jetMaxEta=5.,
            jetId=JetSelection.TIGHT,
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
            jetId=JetSelection.TIGHT,
            storeKinematics=[],
            outputName="selectedL1PreFiringJets_nominal",
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
            metInput=met_variable[year],
            globalOptions=globalOptions,
            outputName="EventObservables_nominal"
        )
    )


    analyzerChain.append(
        EventSkim(
            selection=lambda event: event.nselectedJets_nominal > 0 and event.nselectedJets_nominal < 5
        )
    )

    analyzerChain.append(
       SimplifiedEventCategorization(
            maxDeltaR=1.3,
            looseLeptons=lambda event: event.subleadingLeptons,
            jetsCollection=lambda event: event.selectedJets_nominal,
            outputName="category_nominal",
            globalOptions=globalOptions
       )
    )

    analyzerChain.append(
        MassReconstruction(
            globalOptions=globalOptions,
            outputName="category_nominal",
            tightLeptons=lambda event: event.leadingLeptons,
            looseLeptons=lambda event: event.subleadingLeptons,
            jets=lambda event: event.selectedJets_nominal,
        )
    )

    analyzerChain.append(
        XGBEvaluation(
            systematics=["nominal"],
            jetCollections=[lambda event: event.selectedJets_nominal],
            modelPath=BDTmodelPath[year],
            inputFeatures="${CMSSW_BASE}/src/PhysicsTools/NanoAODTools/data/bdt/201117/nominal/bdt_inputs.txt",
            outputName="bdt_score"
        )
    )

    analyzerChain.append(
        XGBEvaluation(
            systematics=["nominal"],
            jetCollections=[lambda event: event.selectedJets_nominal],
            modelPath=BDTmodelPath[year],
            inputFeatures="${CMSSW_BASE}/src/PhysicsTools/NanoAODTools/data/bdt/201117/uncorrelated/bdt_inputs.txt",
            outputName="bdt_score_uncorr"
        )
    )

    analyzerChain.append(
        TaggerEvaluationProfiled(
            modelPath=modelPath[year],
            featureDictFile=featureDictFile,
            inputCollections=[
                lambda event: event.category_nominal_jets,
            ],
            taggerName="llpdnnx",
            profiledLabelDict = {
                'LLP_Q': ['LLP_Q','LLP_QTAU_H','LLP_QTAU_3H'],
                'LLP_QE': [ 'LLP_QE'],
                'LLP_QMU': [ 'LLP_QMU']
            },
            globalOptions=globalOptions,
            evalValues = np.linspace(-1.9,1.9,5*4),
        )
    )

    analyzerChain.append(
        JetTaggerProfiledResult(
        inputCollection= lambda event: event.category_nominal_jets,
        taggerName="llpdnnx",
        outputName="category_nominal_jets",
        profiledLabels = ['LLP_Q','LLP_QE','LLP_QMU'],
        globalOptions=globalOptions
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

    if isSignal:
        for coupling in range(1,68):
            storeVariables.append([
                lambda tree, coupling=coupling: tree.branch('LHEWeights_coupling_%i'%coupling,'F'),
                lambda tree, event, coupling=coupling: tree.fillBranch('LHEWeights_coupling_%i'%coupling,getattr(event,"LHEWeights_coupling_%i"%coupling)),
            ])


if testMode:
    for ianalyzer, analyzer in enumerate(analyzerChain):
        if type(analyzer).__name__ == "PileupWeight":
            analyzerChain.pop(ianalyzer)

taggerTypes = ['TaggerEvaluationProfiled', 'XGBEvaluation']
if args.noTagger:
    analyzerChain = ([module for module in analyzerChain if type(module).__name__ not in taggerTypes])

analyzerChain.append(EventInfo(storeVariables=storeVariables))

p = PostProcessor(
    args.output[0],
    [args.inputFiles],
    cut="(nJet>0)&&((nElectron+nMuon)>0)",
    modules=analyzerChain,
    maxEvents=-1,
    friend=True
)

p.run()
