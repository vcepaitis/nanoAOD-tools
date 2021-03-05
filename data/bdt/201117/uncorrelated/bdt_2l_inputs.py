import math

def jetSyst(syst):
    if syst in ['unclEnUp','unclEnDown']:
        return 'nominal'
    else:
        return syst
       
    

features = [
    ('met', lambda event, syst: getattr(event,syst+"_met")), 
    ('dilepton_charge', lambda event, syst: event.dilepton_charge), 
    ('l1_eta', lambda event, syst: math.fabs(event.leadingLeptons[0].eta)), 
    ('l1_isElectron', lambda event, syst: event.leadingLeptons[0].isElectron),
    ('met_l1_dphi', lambda event, syst: math.fabs(getattr(event,syst+"_met_l1_dphi"))),
    ('mtw', lambda event, syst: getattr(event,syst+"_mtw")),
    ('l1_pt', lambda event, syst: event.leadingLeptons[0].pt),
    ('j1_eta', lambda event, syst: -1. if len(getattr(event,"selectedJets_"+jetSyst(syst)))<1 else math.fabs(getattr(event,"selectedJets_"+jetSyst(syst))[0].eta)),
    ('l2_eta', lambda event, syst: math.fabs(event.subleadingLeptons[0].eta)),
    ('l2_isElectron', lambda event, syst: event.subleadingLeptons[0].isElectron),
    ('l2_pt', lambda event, syst: event.subleadingLeptons[0].pt),
    ('j2_eta', lambda event, syst: -1. if len(getattr(event,"selectedJets_"+jetSyst(syst)))<2 else math.fabs(getattr(event,"selectedJets_"+jetSyst(syst))[1].eta)),
]

