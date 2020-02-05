import pickle
import pandas as pd
from xgboost import XGBClassifier, Booster, DMatrix
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
        data = pd.DataFrame(data={
                "lepJet_llpdnnx_-1_isLLP_QMU_QQMU" : event.lepJet_llpdnnx_1_isLLP_QMU_QQMU,
                "lepJet_llpdnnx_0_isLLP_QMU_QQMU" : event.lepJet_llpdnnx_0_isLLP_QMU_QQMU,
                "lepJet_llpdnnx_1_isLLP_QMU_QQMU" : event.lepJet_llpdnnx_1_isLLP_QMU_QQMU,
                "lepJet_llpdnnx_2_isLLP_QMU_QQMU" : event.lepJet_llpdnnx_2_isLLP_QMU_QQMU,
                "dimuon_mass" : event.dimuon_mass,
                "dimuon_deltaR" : event.dimuon_deltaR,
                "lepJet_pt" : event.lepJet_pt,
                "lepJet_eta" : event.lepJet_eta,
                "lepJet_deltaR" : event.lepJet_deltaR,
                "MET_pt" : event.MET_pt,
                "MET_phi" : event.MET_phi,
                "looseMuons_pt" : event.looseMuons_pt,
                "looseMuons_eta" : event.looseMuons_eta,
                "looseMuons_dxy" : event.looseMuons_dxy,
                "tightMuons_pt" : event.tightMuons_pt,
                "tightMuons_eta" : event.tightMuons_eta,
                "tightMuons_dxy" : event.tightMuons_dxy,
                },index=[0])
        model = XGBClassifier()
        booster = Booster()
        #model._le = LabelEncoder().fit([1])
        booster.load_model(self.modelPath)
        model._Booster = booster
        bdt_score = model.predict_proba(data)[:, 1]
        setattr(event, "bdt_score", bdt_score)
        return True
        
