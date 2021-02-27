import os
import pickle
import pandas as pd
import numpy
from xgboost import XGBClassifier, Booster, DMatrix
from sklearn.preprocessing import LabelEncoder
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module

class XGBEvaluation(Module):
    def __init__(
        self,
        modelPath="${CMSSW_BASE}/src/PhysicsTools/NanoAODTools/data/bdt/201117/nominal/bdt_2016.model",
        inputFeatures="${CMSSW_BASE}/src/PhysicsTools/NanoAODTools/data/bdt/201117/nominal/bdt_inputs.txt",
        systematics=["nominal"],
        jetCollections= lambda event: event.selectedJets_nominal,
        leadingLeptonCollection= lambda event: event.leadingLeptons,
        subleadingLeptonCollection= lambda event: event.subleadingLeptons,
        outputName ="bdt_score"
    ):
        self.modelPath = modelPath
        self.systematics = systematics
        self.jetCollections = jetCollections
        self.leadingLeptonCollection = leadingLeptonCollection
        self.subleadingLeptonCollection = subleadingLeptonCollection
        self.outputName = outputName
        
        with open(os.path.expandvars(inputFeatures)) as f:
            self.array_bdt = sorted([line.rstrip() for line in f])
            print("BDT input features: ", self.array_bdt)

    def beginJob(self):

        pass
    def endJob(self):
        pass
    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.model = XGBClassifier()
        # self.model.load_model(os.path.expandvars(self.modelPath))
        self.booster = Booster(model_file=os.path.expandvars(self.modelPath))
        self.model._Booster = self.booster
        self.out = wrappedOutputTree
        for systematic in self.systematics:
            self.out.branch(self.outputName+"_"+systematic, "F")

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def analyze(self, event):
        
        leadingLeptons = self.leadingLeptonCollection(event)
        subleadingLeptons = self.subleadingLeptonCollection(event)

        i = 0
        for jetCollection, systematic in zip(self.jetCollections, self.systematics):
            jets = jetCollection(event)
            dict_list = {}
            dict_list["EventObservables_nominal_met"] = getattr(event, "EventObservables_"+systematic+"_met")
            dict_list["EventObservables_nominal_ht"] = getattr(event, "EventObservables_"+systematic+"_ht")
            dict_list["dilepton_mass"] = event.dilepton_mass
            dict_list["dilepton_charge"] = event.dilepton_charge
            dict_list["dilepton_deltaPhi"] = event.dilepton_deltaPhi
            dict_list["dilepton_deltaR"] = event.dilepton_deltaR
            dict_list["leadingLeptons_pt"] = leadingLeptons[0].pt
            dict_list["leadingLeptons_eta"] = leadingLeptons[0].eta
            dict_list["leadingLeptons_nominal_mtw"] = getattr(event, "leadingLeptons_"+systematic+"_mtw")
            dict_list["leadingLeptons_nominal_deltaPhi"] = getattr(event, "leadingLeptons_"+systematic+"_deltaPhi")
            dict_list["subleadingLeptons_pt"] = subleadingLeptons[0].pt if len(subleadingLeptons) > 0 else 0
            dict_list["subleadingLeptons_eta"] = subleadingLeptons[0].eta if len(subleadingLeptons) > 0 else 0
            dict_list["selectedJets_nominal_pt"] = jets[0].pt if len(jets) > 0 else 0
            dict_list["selectedJets_nominal_eta"] = jets[0].eta if len(jets) > 0 else 0
            dict_list["subselectedJets_nominal_pt"] = jets[1].pt if len(jets) > 1 else -1.
            dict_list["subselectedJets_nominal_eta"] = jets[1].eta if len(jets) > 1 else -1.
            dict_list["leadingLeptons_isElectron"] = leadingLeptons[0].isElectron
            dict_list["subleadingLeptons_isElectron"] = subleadingLeptons[0].isElectron if len(subleadingLeptons) > 0 else -1.
            dict_list["category_simplified_nominal_llpdnnx_m_llj"] = getattr(event, "category_"+systematic+"_m_llj")

            if i == 0:
                data = pd.DataFrame(data=dict_list, index=[0])
            else:
                _data = pd.DataFrame(data=dict_list, index=[0])
                data = data.append(_data, ignore_index=True)

            i += 1
        
        data = data[self.array_bdt]
        bdt_score = self.model.predict_proba(data)
        for i, systematic in enumerate(self.systematics):
            setattr(event, self.outputName+"_"+systematic, bdt_score[i, 1])
            self.out.fillBranch(self.outputName+"_"+systematic,bdt_score[i, 1])

        return True
