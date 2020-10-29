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
        modelPath="${CMSSW_BASE}/src/PhysicsTools/NanoAODTools/data/bdt/bdt.model",
        featurePath="${CMSSW_BASE}/src/PhysicsTools/NanoAODTools/data/bdt/bdt_inputs.txt",
        systematics=["nominal"],
        jetCollections= lambda event: event.selectedJets_nominal,
        leadingLeptonCollection= lambda event: event.leadingLeptons,
        subleadingLeptonCollection= lambda event: event.subleadingLeptons,
        outputName ="bdt_score"
    ):
        self.modelPath = modelPath
        self.featurePath = featurePath
        self.systematics = systematics
        self.jetCollections = jetCollections
        self.leadingLeptonCollection = leadingLeptonCollection
        self.subleadingLeptonCollection = subleadingLeptonCollection
        self.outputName = outputName

    def beginJob(self):
        pass
    def endJob(self):
        pass
    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        for systematic in self.systematics:
            self.out.branch(self.outputName+"_"+systematic, "F")

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def analyze(self, event):
        
        leadingLeptons = self.leadingLeptonCollection(event)
        subleadingLeptons = self.subleadingLeptonCollection(event)
        '''
        if len(jets) == 0 or len(leadingLeptons) == 0 or len(subleadingLeptons) == 0:
            setattr(event, "bdt_score", "-999.")
            self.out.fillBranch(self.outputName,-999.)
            return True
        '''

        '''
        with open(os.path.expandvars(self.featurePath)) as f:
            array_list = [line.rstrip() for line in f]
        '''
        i = 0
        for jetCollection, systematic in zip(self.jetCollections, self.systematics):
            jets = jetCollection(event)
            dict_list = {}
            dict_list["EventObservables_nominal_met"] = getattr(event, "EventObservables_"+systematic+"_met")
            dict_list["EventObservables_nominal_ht"] = getattr(event, "EventObservables_"+systematic+"_ht")
            dict_list["dilepton_mass"] = event.dilepton_mass
            dict_list["dilepton_deltaPhi"] = event.dilepton_deltaPhi
            dict_list["dilepton_deltaR"] = event.dilepton_deltaR
            dict_list["leadingLeptons_pt"] = leadingLeptons[0].pt
            dict_list["leadingLeptons_eta"] = leadingLeptons[0].eta
            dict_list["leadingLeptons_nominal_mtw"] = getattr(event, "leadingLeptons_"+systematic+"_mtw")
            dict_list["leadingLeptons_nominal_deltaPhi"] = getattr(event, "leadingLeptons_"+systematic+"_deltaPhi")
            dict_list["subleadingLeptons_pt"] = subleadingLeptons[0].pt if len(subleadingLeptons) > 0 else 0
            dict_list["subleadingLeptons_eta"] = subleadingLeptons[0].eta if len(subleadingLeptons) > 0 else 0
            dict_list["selectedJets_nominal_ptLeptonSubtracted"] = jets[0].pt if len(jets) > 0 else 0
            dict_list["selectedJets_nominal_eta"] = jets[0].eta if len(jets) > 0 else 0
            if i == 0:
                data = pd.DataFrame(data=dict_list, index=[0])
            else:
                _data = pd.DataFrame(data=dict_list, index=[0])
                data = data.append(_data, ignore_index=True)

            i += 1
        data = data.reindex(sorted(data.columns), axis=1)
        model = XGBClassifier()
        booster = Booster(model_file=os.path.expandvars(self.modelPath))
        model._Booster = booster
        bdt_score = model.predict_proba(data)
        for i, systematic in enumerate(self.systematics):
            setattr(event, "bdt_score", bdt_score[i, 1])
            self.out.fillBranch(self.outputName+"_"+systematic,bdt_score[i, 1])

        return True
