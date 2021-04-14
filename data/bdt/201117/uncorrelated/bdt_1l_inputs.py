import math

def jetSyst(syst):
    if syst in ['unclEnUp','unclEnDown']:
        return 'nominal'
    else:
        return syst
       
    

features = [
    ('mtw', lambda event, syst: getattr(event,syst+"_mtw")),
    ('dPhi_met_l1', lambda event, syst: math.fabs(getattr(event,syst+"_dPhi_met_l1"))),
    ('ht', lambda event, syst: getattr(event,syst+"_ht")),
    ('lj_dphi', lambda event, syst: getattr(event,syst+"_dPhi_l1j")),
    ('lj_deta', lambda event, syst: getattr(event,syst+"_dEta_l1j")),
    ('circularity', lambda event, syst: getattr(event,syst+"_eventShape_circularity")),
    ('isotropy', lambda event, syst: getattr(event,syst+"_eventShape_isotropy")),
    ('dPhi_mht_l', lambda event, syst: getattr(event,syst+"_dPhi_mht_l1")),
    ('ptR_mht_l', lambda event, syst: getattr(event,syst+"_ptR_mht_l1")),
    ('dPhi_mht_met', lambda event, syst: getattr(event,syst+"_dPhi_mht_met")),
]

