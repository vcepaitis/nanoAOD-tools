import os
import pickle
print pickle.__version__
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
        outputName ="bdt_score"
    ):
        self.modelPath = modelPath
        self.featurePath = featurePath
        self.outputName = outputName

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

        with open(os.path.expandvars(self.featurePath)) as f:
            array_list = [line.rstrip() for line in f]

        if event.nselectedJets_nominal == 0 or event.leadingLepton == 0 or event.nsubleadingLepton == 0:
            setattr(event, "bdt_score", "-999.")
            self.out.fillBranch(self.outputName,-999.)
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
        data = data.reindex(sorted(data.columns), axis=1)
        model = XGBClassifier()
        booster = Booster(model_file=os.path.expandvars(self.modelPath))
        model._Booster = booster
        bdt_score = model.predict_proba(data)
        setattr(event, "bdt_score", bdt_score[:, 1])
        self.out.fillBranch(self.outputName,bdt_score[:, 1])

        return True
