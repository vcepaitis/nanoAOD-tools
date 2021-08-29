import math

def jetSyst(syst):
    if syst in ['unclEnUp','unclEnDown']:
        return 'nominal'
    else:
        return syst

features = [
    ('leadingLeptons_pt', lambda event, syst: event.leadingLeptons[0].pt),
    ('nominal_met', lambda event, syst: getattr(event,syst+"_met")), 
    ('nominal_mtw', lambda event, syst: getattr(event,syst+"_mtw")),
    ('nominal_dPhi_met_l1', lambda event, syst: math.fabs(getattr(event,syst+"_dPhi_met_l1"))),
    ('leadingLeptons_eta', lambda event, syst: math.fabs(event.leadingLeptons[0].eta)), 
    ('dilepton_charge', lambda event, syst: event.dilepton_charge), 
    ('leadingLeptons_isElectron', lambda event, syst: event.leadingLeptons[0].isElectron),
    ('subleadingLeptons_isElectron', lambda event, syst: event.subleadingLeptons[0].isElectron),
    ('nominal_eventShape_aplanarity', lambda event, syst: getattr(event,syst+"_eventShape_aplanarity")),
]

features = sorted(features, key=lambda feature_tuple: feature_tuple[0])
