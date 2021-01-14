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
parser.add_argument('--inv', dest='invertLeptons', action='store_true', default=False)
parser.add_argument('output', nargs=1)

args = parser.parse_args()

testMode = args.testMode
print "isData:",args.isData
print "isSignal:",args.isSignal
print "inputs:",len(args.inputFiles)

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
    "year": year,
    "isSignal":args.isSignal
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


leptonSelection = [
    #EventSkim(selection=lambda event: event.nTrigObj > 0),
    MuonSelection(
        outputName="tightMuon",
        storeKinematics=[],
        storeWeights=True,
        muonMinPt=minMuonPt[globalOptions["year"]],
        triggerMatch=True,
        muonID=MuonSelection.TIGHT,
        muonIso=MuonSelection.INV if args.invertLeptons else MuonSelection.TIGHT,
        selectLeadingOnly=True,
        globalOptions=globalOptions
    ),
    ElectronSelection(
        outputName="tightElectron",
        storeKinematics=[],
        electronMinPt=minElectronPt[globalOptions["year"]],
        electronID="Inv" if args.invertLeptons else "Iso_WP90",
        storeWeights=True,
        triggerMatch=True,
        selectLeadingOnly=True,
        globalOptions=globalOptions
    ),
    #EventSkim(selection=lambda event: (event.ntightMuon + event.ntightElectron) > 0),
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
    #EventSkim(selection=lambda event: (event.IsoMuTrigger_flag + event.IsoElectronTrigger_flag) > 0),
    MuonSelection(
        inputCollection=lambda event: event.tightMuon_unselected,
        outputName="looseMuons",
        storeKinematics=[],
        muonMinPt=5.,
        muonID=MuonSelection.LOOSE,
        muonIso=MuonSelection.NONE,
        globalOptions=globalOptions
    ),
    ElectronSelection(
        inputCollection=lambda event: event.tightElectron_unselected,
        outputName="looseElectrons",
        storeKinematics=[],

        electronMinPt=5.,
        electronID="Custom",
        globalOptions=globalOptions
    ),

    #EventSkim(selection=lambda event: (event.ntightMuon + event.ntightElectron + event.nlooseMuons + event.nlooseElectrons ) <= 2),
    LeptonCollecting(
        tightMuonCollection=lambda event:event.tightMuon,
        tightElectronCollection=lambda event:event.tightElectron,
        looseMuonCollection=lambda event:event.looseMuons,
        looseElectronCollection=lambda event:event.looseElectrons
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

featureDictFile = "${CMSSW_BASE}/src/PhysicsTools/NanoAODTools/data/nn/201117/feature_dict.py"
taggers = [
    {"name":"nominal_ref", "modelfile":"${CMSSW_BASE}/src/PhysicsTools/NanoAODTools/data/nn/201117/weightMixed%i_NominalNetwork_ref_201117.pb"%year},
]

if year==2016:
    taggers.extend([
        {"name":"nominal", "modelfile":"${CMSSW_BASE}/src/PhysicsTools/NanoAODTools/data/nn/201117/weightMixed2016_NominalNetwork_201117.pb"},
        #{"name":"attention", "modelfile":"${CMSSW_BASE}/src/PhysicsTools/NanoAODTools/data/nn/201117/weightMixed2016_AttentionNetwork_201117.pb"},
        #{"name":"deepset", "modelfile":"${CMSSW_BASE}/src/PhysicsTools/NanoAODTools/data/nn/201117/weightMixed2016_DeepSetNetwork_201117.pb"},
        #{"name":"onlyglobal", "modelfile":"${CMSSW_BASE}/src/PhysicsTools/NanoAODTools/data/nn/201117/weightMixed2016_NominalNetworkOnlyGlobal_201117.pb"},
        #{"name":"p4attention", "modelfile":"${CMSSW_BASE}/src/PhysicsTools/NanoAODTools/data/nn/201117/weightMixed2016_P4AttentionNetwork_201117.pb"},
        #{"name":"p4attentionv2", "modelfile":"${CMSSW_BASE}/src/PhysicsTools/NanoAODTools/data/nn/201117/weightMixed2016_P4AttentionNetworkv2_201117.pb"},    
    ])

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

met_variable = {
    2016: lambda event: Object(event, "MET"),
    2017: lambda event: Object(event, "METFixEE2017"),
    2018: lambda event: Object(event, "MET")
}


analyzerChain.append(
     MetFilter(
        globalOptions=globalOptions,
        outputName="MET_filter"
     )
)


if not isMC:
    print "ERROR - performance can only be evaluated on MC"
    sys.exit(1)


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
        propagateJER = False, # need to fix, poor modelling
        jetKeys = ['pt', 'eta', 'phi' , 'jetId', 'nConstituents'],
    )
)

for systName, jetCollection in [
    ("nominal", lambda event: event.jets_nominal),
    #("jerUp", lambda event: event.jets_jerUp),
    #("jerDown", lambda event: event.jets_jerDown),
    #("jesTotalUp", lambda event: event.jets_jesTotalUp),
    #("jesTotalDown", lambda event: event.jets_jesTotalDown),
]:

    analyzerChain.append(
        JetSelection(
            inputCollection=jetCollection,
            leptonCollectionDRCleaning=lambda event: event.leadingLeptons,
            leptonCollectionP4Subraction=lambda event: event.subleadingLeptons,
            jetMinPt=15.,
            jetMaxEta=2.395,
            jetMinNConstituents=2,
            jetId=JetSelection.LOOSE,
            storeKinematics=[
                'pt', 'eta', 'phi',
                'nConstituents',
                'minDPhiClean','minDRClean'
            ],
            #dRCleaning = -1. if args.invertLeptons else 0.4,
            outputName="selectedJets_nominal",
            globalOptions=globalOptions
        )
    )


    analyzerChain.append(
        JetTruthFlags(
            inputCollection=lambda event, systName=systName: getattr(event, "selectedJets_"+systName),
            outputName="selectedJets_"+systName,
            originVariables=['displacement','displacement_xy','displacement_z'],
            globalVariables=["numberCpf","numberNpf","numberSv","numberMuon","numberElectron"],
            flags={
                'isE': ['isPrompt_E'],
                'isMU': ['isPrompt_MU'],
                'isTAU': ['isPrompt_TAU'],

                'isB': ['isB', 'isBB', 'isLeptonic_B'],
                'isC': ['isC', 'isCC', 'isLeptonic_C'],
                'isUDS': ['isS', 'isUD'],
                'isG': ['isG'],
                'isPU': ['isPU'],

                'isLLP_Q': ['isLLP_Q','isLLP_QQ','isLLP_RAD'],

                'isLLP_QE': ['isLLP_E','isLLP_QE','isLLP_QQE'],
                'isLLP_QMU': ['isLLP_MU','isLLP_QMU','isLLP_QQMU'],
                'isLLP_QTAU': ['isLLP_TAU','isLLP_QTAU','isLLP_QQTAU'],

            },
            globalOptions=globalOptions
        )


    )


analyzerChain.append(
    EventSkim(selection=lambda event: \
        getattr(event, "nselectedJets_nominal") > 0
        #or getattr(event, "nselectedJets_jesTotalUp") > 0 \
        #or getattr(event, "nselectedJets_jesTotalDown") > 0 \
        #or getattr(event, "nselectedJets_jerUp") > 0 \
        #or getattr(event, "nselectedJets_jerDown") > 0
    )
)



for tagger in taggers:
    analyzerChain.append(
        TaggerEvaluationProfiled(
            modelPath=tagger["modelfile"],
            featureDictFile=featureDictFile,
            inputCollections=[
                lambda event: event.selectedJets_nominal[:4],
                #lambda event: event.selectedJets_jesTotalUp[:4],
                #lambda event: event.selectedJets_jesTotalDown[:4],
                #lambda event: event.selectedJets_jerUp[:4],
                #lambda event: event.selectedJets_jerDown[:4],
            ],
            taggerName=tagger["name"],
            profiledLabelDict = {
                'E':['E'],
                'MU':['MU'],
                'TAU':['TAU'],
                'UDS':['UDS'],
                'G':['G'],
                'B':['B'],
                'C':['C'],
                'PU':['PU'],
                'LLP_Q': ['LLP_Q'],
                'LLP_QE': ['LLP_QE'],
                'LLP_QMU': ['LLP_QMU'],
                'LLP_QTAU_H': ['LLP_QTAU_H'],
                'LLP_QTAU_3H': ['LLP_QTAU_3H'],
                'TAUANY':['TAU','LLP_QTAU_H','LLP_QTAU_3H']
            },
            globalOptions=globalOptions,
            evalValues = np.linspace(-1.9,1.9,5*4)#np.linspace(-2,2,5*4+1),
            
            #np.linspace(-2,2,5*4+1) #np.linspace(-1.1,2.1,17),
        )
    )




for systName, jetCollection, metObject in [
    ("nominal", lambda event: event.selectedJets_nominal,
        lambda event: event.met_nominal),
    #("jerUp", lambda event: event.selectedJets_jerUp,
    #    lambda event: event.met_jerUp),
    #("jerDown", lambda event: event.selectedJets_jerDown,
    #    lambda event: event.met_jerDown),
    #("jesTotalUp", lambda event: event.selectedJets_jesTotalUp,
    #    lambda event: event.met_jesTotalUp),
    #("jesTotalDown", lambda event: event.selectedJets_jesTotalDown,
    #    lambda event: event.met_jesTotalDown),
    #("unclEnUp", lambda event: event.selectedJets_nominal,
    #    lambda event: event.met_unclEnUp),
    #("unclEnDown", lambda event: event.selectedJets_nominal,
    #    lambda event: event.met_unclEnDown),
]:
    for tagger in taggers:
        analyzerChain.append(
            JetTaggerProfiledResult(
                inputCollection=jetCollection,
                taggerName=tagger["name"],
                outputName="selectedJets_"+systName,
                profiledLabels = ['E','MU','TAU','UDS','G','PU','B','C','TAUANY','LLP_Q','LLP_QE','LLP_QMU','LLP_QTAU_H','LLP_QTAU_3H'],
                globalOptions=globalOptions
            )
        )
    '''
    analyzerChain.extend([
        WbosonReconstruction(
            leptonCollectionName='leadingLepton',
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
    '''

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


p = PostProcessor(
    args.output[0],
    [args.inputFiles],
    cut="1",#"((nElectron+nMuon)>0)",
    modules=analyzerChain,
    maxEvents=25000,
    friend=True
)

p.run()
