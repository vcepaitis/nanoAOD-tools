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
        modelPath,
    ):
        self.modelPath = modelPath

    def beginJob(self):
        pass
    def endJob(self):
        pass
    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def analyze(self, event):

        with open(os.path.expandvars('${CMSSW_BASE}/src/PhysicsTools/NanoAODTools/data/bdt/bdt_inputs.txt')) as f:
            array_list = [line.rstrip() for line in f]

        if event.nselectedJets_nominal == 0 or event.ntightMuon == 0 or event.nlooseMuons == 0:
            setattr(event, "bdt_score", "-999.")
            return True

        dict_list = {}
        for feature in array_list:
            value = getattr(event, feature)
            if isinstance(value, list):
                if len(value) > 0:
                    dict_list[feature] = value[0]
            else:
                dict_list[feature] = value
        data = pd.DataFrame(data=dict_list, index=[0])

        #print dict_list

        model = XGBClassifier()
        booster = Booster()
        #model._le = LabelEncoder().fit([1])
        booster.load_model(self.modelPath)
        booster.feature_names = sorted(array_list)
        model._Booster = booster
        bdt_score = model.predict_proba(data)
        setattr(event, "bdt_score", bdt_score[:, 1])
        return True
