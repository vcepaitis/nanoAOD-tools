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

parser.add_argument('--isData', dest='isData',
                    action='store_true', default=False)
parser.add_argument('--isSignal', dest='isSignal',
                    action='store_true', default=False)
parser.add_argument('--year', dest='year',
                    action='store', type=int, default=-1)
parser.add_argument('--noiso', dest='noiso',
                    action='store_true', default=False)
parser.add_argument('--notrigger', dest='notrigger',
                    action='store_true', default=False)  
parser.add_argument('--notagger', dest='notagger',
                    action='store_true', default=False)  
parser.add_argument('--nobdt', dest='nobdt',
                    action='store_true', default=False) 
parser.add_argument('--leptons', dest='leptons', type=int, default=2, choices=[1,2])                     
parser.add_argument('--input', dest='inputFiles', action='append', default=[])
parser.add_argument('output', nargs=1)

args = parser.parse_args()

print "isData:",args.isData
print "inputs:",len(args.inputFiles)
isSignal = False

for inputFile in args.inputFiles:
    if args.isSignal or "dirac" in inputFile or "majorana" in inputFile or "LLPGun" in inputFile: 
        isSignal = True
    if args.year<0 and ("-2016" in inputFile or "Run2016" in inputFile):
        year = 2016
    elif args.year<0 and ("-2017" in inputFile or "Run2017" in inputFile):
        year = 2017
    elif args.year<0 and ("-2018" in inputFile or "Run2018" in inputFile):
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
print "isSignal:",isSignal
print "apply lepton iso: ","True (default)" if not args.noiso else "False"
print "apply trigger selection: ","True (default)" if not args.notrigger else "False" 
print "run BDT: ","True (default)" if not args.nobdt else "False" 
print "run tagger: ","True (default)" if not args.notagger else "False" 
print "channel: ","single lepton" if args.leptons==1 else "dilepton"
print "output directory:", args.output[0]

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


leptonSelection = [
    MuonSelection(
        outputName="tightMuons",
        storeKinematics=[],
        storeWeights=True,
        muonMinPt=minMuonPt[globalOptions["year"]],
        muonMaxDxy=0.01,
        muonMaxDz=0.05,
        muonID=MuonSelection.TIGHT,
        muonIso=MuonSelection.NONE if args.noiso else MuonSelection.TIGHT,
        globalOptions=globalOptions
    ),
    ElectronSelection(
        outputName="tightElectrons",
        storeKinematics=[],
        electronMinPt=minElectronPt[globalOptions["year"]],
        electronID="noIso_WP90" if args.noiso else "Iso_WP90",
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
        storeKinematics=[],
        muonMinPt=3.,
        muonID=MuonSelection.LOOSE,
        muonIso=MuonSelection.NONE,
        globalOptions=globalOptions
    ),
    MuonSelection(
        inputCollection=lambda event: event.looseMuons,
        outputName="looseIsoMuons",
        storeKinematics=[],
        muonMinPt=3.,
        storeWeights=False,
        muonID=MuonSelection.LOOSE,
        muonIso=MuonSelection.TIGHT,
        globalOptions=globalOptions
    ),
    ElectronSelection(
        inputCollection=lambda event: event.tightElectrons_unselected,
        outputName="looseElectrons",
        storeKinematics=[],

        electronMinPt=5.,
        electronID="Custom",
        globalOptions=globalOptions
    ),
    ElectronSelection(
        inputCollection=lambda event: event.looseElectrons,
        outputName="looseIsoElectrons",
        storeKinematics=[],
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
        outputName = "Leptons",
        storeLeadingKinematics=["pt", "eta", "phi", "charge/I", "isMuon/I", "isElectron/I", "relIso"],
        storeSubleadingKinematics=["pt", "eta", "phi", "charge/I", "isMuon/I", "isElectron/I", "relIso", "dxy", "dz", 'dxysig', 'dzsig']
    ),
    
]

analyzerChain = []

analyzerChain.extend(leptonSelection)

if not args.notrigger:
    analyzerChain.extend([
        EventSkim(selection=lambda event: (event.IsoMuTrigger_flag + event.IsoElectronTrigger_flag) > 0),
        EventSkim(selection=lambda event: event.leadingLeptons[0].isTriggerMatched>0),
    ])
    
if args.leptons==1:
    analyzerChain.append(
        EventSkim(selection=lambda event: event.nsubleadingLeptons==0)
    )
else:
    analyzerChain.extend([
        EventSkim(selection=lambda event: event.nsubleadingLeptons==1),
        InvariantSystem(
            inputCollection= lambda event: [event.leadingLeptons[0],event.subleadingLeptons[0]],
            outputName="dilepton"
        )
    ])




featureDictFile = "${CMSSW_BASE}/src/PhysicsTools/NanoAODTools/data/nn/201117/feature_dict.py"
taggerModelPath = {
    2016: "${CMSSW_BASE}/src/PhysicsTools/NanoAODTools/data/nn/201117/weightMixed2016_ExtNominalNetwork_allflavour_origSV_DA_30_wasserstein4_lr001_201117.pb",
    2017: "${CMSSW_BASE}/src/PhysicsTools/NanoAODTools/data/nn/201117/weightMixed2017_ExtNominalNetwork_allflavour_origSV_DA_30_wasserstein4_lr001_201117.pb",
    2018: "${CMSSW_BASE}/src/PhysicsTools/NanoAODTools/data/nn/201117/weightMixed2018_ExtNominalNetwork_allflavour_origSV_DA_30_wasserstein4_lr001_201117.pb"
}

BDT2lmodelPath = {
    2016: "${CMSSW_BASE}/src/PhysicsTools/NanoAODTools/data/bdt/201117/nominal/bdt_2016.model",
    2017: "${CMSSW_BASE}/src/PhysicsTools/NanoAODTools/data/bdt/201117/nominal/bdt_2017.model",
    2018: "${CMSSW_BASE}/src/PhysicsTools/NanoAODTools/data/bdt/201117/nominal/bdt_2018.model"
}

BDT2lmodelPathUncorr = {
    2016: "${CMSSW_BASE}/src/PhysicsTools/NanoAODTools/data/bdt/201117/uncorrelated/bdt_2l_2016.model",
    2017: "${CMSSW_BASE}/src/PhysicsTools/NanoAODTools/data/bdt/201117/uncorrelated/bdt_2l_2017.model",
    2018: "${CMSSW_BASE}/src/PhysicsTools/NanoAODTools/data/bdt/201117/uncorrelated/bdt_2l_2018.model"
}

BDT1lmodelPathUncorr = {
    2016: "${CMSSW_BASE}/src/PhysicsTools/NanoAODTools/data/bdt/201117/uncorrelated/bdt_1l_lr0.100_min0.1000_depths4_bagging0.75_2016.model",
    2017: "${CMSSW_BASE}/src/PhysicsTools/NanoAODTools/data/bdt/201117/uncorrelated/bdt_1l_lr0.100_min0.1000_depths4_bagging0.75_2017.model",
    2018: "${CMSSW_BASE}/src/PhysicsTools/NanoAODTools/data/bdt/201117/uncorrelated/bdt_1l_lr0.100_min0.1000_depths4_bagging0.75_2018.model",
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


def jetSelectionSequence(jetDict):
    sequence = []
    for systName,jetCollection in jetDict.items():
        sequence.extend([
            JetSelection(
                inputCollection=jetCollection,
                leptonCollectionDRCleaning=lambda event: event.tightMuons+event.tightElectrons ,
                leptonCollectionP4Subraction=lambda event:event.looseMuons+event.looseElectrons,
                jetMinPt=15.,
                jetMaxEta=2.4,
                globalFeatures = [],
                storeKinematics=['pt', 'eta', 'phi', 'minDeltaRSubtraction', 'ptLepton', 'ptOriginal', 'ptSubtracted', 'rawFactor', 'ptRaw'],
                jetId=JetSelection.TIGHT,
                outputName="selectedJets_"+systName,
                globalOptions=globalOptions
            ),
            JetSelection(
                inputCollection=jetCollection,
                leptonCollectionDRCleaning=lambda event: event.tightMuons+event.tightElectrons,
                leptonCollectionP4Subraction=lambda event:event.looseMuons+event.looseElectrons,
                jetMinPt=30.,
                jetMaxEta=2.4,
                globalFeatures = [],
                storeKinematics=[],
                jetId=JetSelection.TIGHT,
                outputName="selectedPt30Jets_"+systName,
                globalOptions=globalOptions
            ),
            JetSelection(
                inputCollection=jetCollection,
                leptonCollectionDRCleaning=lambda event: event.tightMuons+event.tightElectrons+event.looseIsoMuons+event.looseIsoElectrons,
                jetMinPt=30.,
                jetMinEta=2.4,
                jetMaxEta=5.,
                jetId=JetSelection.TIGHT,
                storeKinematics=[],
                globalFeatures=[],
                outputName="selectedFwdJets_"+systName,
                globalOptions=globalOptions
            ),
            JetSelection(
                inputCollection=jetCollection,
                jetMinPt=100.,
                jetMinEta=2.25,
                jetMaxEta=3.0,
                jetId=JetSelection.TIGHT,
                storeKinematics=[],
                globalFeatures=[],
                outputName="selectedL1PreFiringJets_"+systName,
                globalOptions=globalOptions
            )
        ])
        
    
        if isMC:
            sequence.append(
                JetTruthFlags(
                    inputCollection=lambda event, systName=systName: getattr(event, "selectedJets_"+systName),
                    originVariables = ['displacement_xy'],
                    outputName="selectedJets_"+systName,
                    globalOptions=globalOptions
                )
            )
    
    systNames = jetDict.keys()
    sequence.append(
        EventSkim(selection=lambda event, systNames=systNames: 
            any([getattr(event, "nselectedJets_"+systName) > 0 for systName in systNames])
        )
    )
            
    return sequence
    
    
def eventReconstructionSequence(jetMetDict):
    sequence = []
    for systName,(jetCollection,metObject) in jetMetDict.items():
        sequence.extend([
            EventObservables(
                lepton1Object=lambda event: event.leadingLeptons[0],
                lepton2Object=None if args.leptons==1 else (lambda event: event.subleadingLeptons[0]),
                jetCollection=jetCollection,
                metInput=metObject,
                globalOptions=globalOptions,
                outputName=systName
            ),

            HNLReconstruction(
                lepton1Object=lambda event: event.leadingLeptons[0],
                lepton2Object=None if args.leptons==1 else (lambda event: event.subleadingLeptons[0]),
                jetCollection=jetCollection,
                globalOptions=globalOptions,
                outputName=systName,
            )
        ])
        
    return sequence
    
    
def taggerSequence(jetDict, modelFile, taggerName):
    sequence = []
    if args.notagger:
        return []
    sequence.append(
        TaggerEvaluationProfiled(
            modelPath=taggerModelPath,
            featureDictFile=featureDictFile,
            inputCollections=jetDict.values(),
            taggerName=taggerName,
            profiledLabelDict = {
                'LLP_Q': ['LLP_Q','LLP_QTAU_H','LLP_QTAU_3H'],
                'LLP_QE': [ 'LLP_QE'],
                'LLP_QMU': [ 'LLP_QMU'],
                'LLP_QTAU_3H': ['LLP_QTAU_3H']
            },
            globalOptions=globalOptions,
            evalValues = np.linspace(-1.9,1.9,5*4),
        )
    )
   
    for systName,jetCollection in jetDict.items():
        sequence.append(
            JetTaggerProfiledResult(
                inputCollection=jetCollection,
                taggerName=taggerName,
                outputName="hnlJet_"+systName,
                profiledLabels = ['LLP_Q','LLP_QE','LLP_QMU','LLP_QTAU_3H'],
                kinds = ['single','ratio'],
                globalOptions={"isData": False}
            )
        )
    return sequence

    

def bdtSequence(systematics):
    sequence = []
    if args.nobdt:
        return []
    if args.leptons==1:
        sequence.append(
            XGBEvaluation(
                systematics=systematics,
                modelPath=BDT1lmodelPathUncorr[year],
                inputFeatures="${CMSSW_BASE}/src/PhysicsTools/NanoAODTools/data/bdt/201117/uncorrelated/bdt_1l_inputs.py",
                outputName="bdt_score"
            )
        )
    else:
        sequence.append(
            XGBEvaluation(
                systematics=systematics,
                modelPath=BDT2lmodelPathUncorr[year],
                inputFeatures="${CMSSW_BASE}/src/PhysicsTools/NanoAODTools/data/bdt/201117/uncorrelated/bdt_2l_inputs.py",
                outputName="bdt_score"
            )
        )
    return sequence


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
            jetKeys = ['jetId', 'nConstituents', 'rawFactor'],
        )
    )


    analyzerChain.extend(
        jetSelectionSequence({
            "nominal": lambda event: event.jets_nominal,
            "jerUp": lambda event: event.jets_jerUp,
            "jerDown": lambda event: event.jets_jerDown,
            "jesTotalUp": lambda event: event.jets_jesTotalUp,
            "jesTotalDown": lambda event: event.jets_jesTotalDown,
        })
    )
        
        
    analyzerChain.extend(
        eventReconstructionSequence({
            "nominal": (lambda event: event.selectedJets_nominal, lambda event: event.met_nominal),
            "jerUp": (lambda event: event.selectedJets_jerUp, lambda event: event.met_jerUp),
            "jerDown": (lambda event: event.selectedJets_jerDown, lambda event: event.met_jerDown),
            "jesTotalUp": (lambda event: event.selectedJets_jesTotalUp, lambda event: event.met_jesTotalUp),
            "jesTotalDown": (lambda event: event.selectedJets_jesTotalDown, lambda event: event.met_jesTotalDown),
            "unclEnUp": (lambda event: event.selectedJets_nominal, lambda event: event.met_unclEnUp),
            "unclEnDown": (lambda event: event.selectedJets_nominal, lambda event: event.met_unclEnDown),
        })
    )
    
    
    analyzerChain.extend(
        taggerSequence({
            "nominal": lambda event: event.hnlJets_nominal,
            "jerUp": lambda event: event.hnlJets_jerUp,
            "jerDown": lambda event: event.hnlJets_jerDown,
            "jesTotalUp": lambda event: event.hnlJets_jesTotalUp,
            "jesTotalDown": lambda event: event.hnlJets_jesTotalDown,
            "unclEnUp": lambda event: event.hnlJets_nominal,
            "unclEnDown": lambda event: event.hnlJets_nominal,
        },
        modelPath=taggerModelPath,
        modelFile=modelPath[year],
        taggerName='llpdnnx'
    ))
    
    analyzerChain.extend(
        bdtSequence([
            "nominal", 
            "jerUp", "jerDown", 
            "jesTotalUp", "jesTotalDown", 
            "unclEnUp", "unclEnDown"
        ])
    )
    
    analyzerChain.append(
        PileupWeight(
            outputName="puweight",
            #processName = "TTToSemiLeptonic_TuneCP5_PSweights_13TeV-powheg-pythia8-2016", #override pu profile
            globalOptions=globalOptions
        )
    )


else:
    analyzerChain.extend(
        jetSelectionSequence({
            "nominal": lambda event: Collection(event,"Jet")
        })
    )
    
    analyzerChain.extend(
        eventReconstructionSequence({
            "nominal": (lambda event: event.selectedJets_nominal, met_variable[year]),
        })
    )
    
    analyzerChain.extend(
        taggerSequence({
            "nominal": lambda event: event.hnlJets_nominal,
        },
        modelPath=taggerModelPath,
        modelFile=modelPath[year],
        taggerName='llpdnnx'
    ))
    
    analyzerChain.extend(
        bdtSequence([
            "nominal",
        ])
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

analyzerChain.append(EventInfo(storeVariables=storeVariables))

p = PostProcessor(
    args.output[0],
    [args.inputFiles],
    cut="((nElectron+nMuon)>0)",
    modules=analyzerChain,
    maxEvents=-1,
    friend=True
)

p.run()

