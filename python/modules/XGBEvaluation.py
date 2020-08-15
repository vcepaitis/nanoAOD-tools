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
        systName="nominal",
        jetCollection= lambda event: event.selectedJets_nominal,
        leadingLeptonCollection= lambda event: event.leadingLeptons,
        subleadingLeptonCollection= lambda event: event.subleadingLeptons,
        outputName ="bdt_score"
    ):
        self.modelPath = modelPath
        self.featurePath = featurePath
        self.systName = systName
        self.jetCollection = jetCollection
        self.leadingLeptonCollection = leadingLeptonCollection
        self.subleadingLeptonCollection = subleadingLeptonCollection
        self.outputName = outputName+"_"+systName

    def beginJob(self):
        pass
    def endJob(self):
        pass
    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        self.out.branch(self.outputName,"F")

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def analyze(self, event):
        
        jets = self.jetCollection(event)
        leadingLeptons = self.leadingLeptonCollection(event)
        subleadingLeptons = self.subleadingLeptonCollection(event)

        if len(jets) == 0 or len(leadingLeptons) == 0 or len(subleadingLeptons) == 0:
            setattr(event, "bdt_score", "-999.")
            self.out.fillBranch(self.outputName,-999.)
            return True

        '''
        with open(os.path.expandvars(self.featurePath)) as f:
            array_list = [line.rstrip() for line in f]
        '''

        dict_list = {}
        dict_list["EventObservables_nominal_met"] = getattr(event, "EventObservables_"+self.systName+"_met")
        dict_list["EventObservables_nominal_met_phi"] = getattr(event, "EventObservables_"+self.systName+"_met_phi")
        dict_list["EventObservables_nominal_ht"] = getattr(event, "EventObservables_"+self.systName+"_ht")
        dict_list["EventObservables_nominal_mT_met_lep"] = getattr(event, "EventObservables_"+self.systName+"_mT_met_lep")
        dict_list["selectedJets_nominal_pt"] = jets[0].pt
        dict_list["selectedJets_nominal_eta"] = jets[0].eta
        dict_list["nselectedJets_nominal"] = len(jets)
        dict_list["dilepton_mass"] = event.dilepton_mass
        dict_list["dilepton_deltaPhi"] = event.dilepton_deltaPhi
        dict_list["dilepton_deltaR"] = event.dilepton_deltaR
        dict_list["leadingLeptons_pt"] = leadingLeptons[0].pt
        dict_list["leadingLeptons_eta"] = leadingLeptons[0].eta
        dict_list["leadingLeptons_phi"] = leadingLeptons[0].phi
        dict_list["subleadingLeptons_pt"] = subleadingLeptons[0].pt
        dict_list["subleadingLeptons_eta"] = subleadingLeptons[0].eta


        '''
        for feature in array_list:
            feature = feature.replace("nominal", self.systName)
            value = getattr(event, feature)
            if isinstance(value, list):
                if len(value) > 0:
                    dict_list[feature] = value[0]
            else:
                dict_list[feature] = value
        '''

        data = pd.DataFrame(data=dict_list, index=[0])
        data = data.reindex(sorted(data.columns), axis=1)
        model = XGBClassifier()
        booster = Booster(model_file=os.path.expandvars(self.modelPath))
        model._Booster = booster
        bdt_score = model.predict_proba(data)
        setattr(event, "bdt_score", bdt_score[:, 1])
        self.out.fillBranch(self.outputName,bdt_score[:, 1])

        return True
