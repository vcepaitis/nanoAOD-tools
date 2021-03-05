import math

def jetSyst(syst):
    if syst in ['unclEnUp','unclEnDown']:
        return 'nominal'
    else:
        return syst
       
    

features = [
    ('mtw', lambda event, syst: getattr(event,syst+"_mtw")),
    ('met_l1_dphi', lambda event, syst: math.fabs(getattr(event,syst+"_met_l1_dphi"))),
    ('ht', lambda event, syst: getattr(event,syst+"_ht")),
    ('lj_dphi', lambda event, syst: getattr(event,syst+"_deltaPhi_l1j")),
    ('lj_deta', lambda event, syst: getattr(event,syst+"_deltaEta_l1j")),
    ('circularity', lambda event, syst: getattr(event,syst+"_eventShape_circularity")),
    ('isotropy', lambda event, syst: getattr(event,syst+"_eventShape_isotropy")),
    ('mht_l1_dphi', lambda event, syst: getattr(event,syst+"_mht_l1_dphi")),
    ('mht_l1_ptRatio', lambda event, syst: getattr(event,syst+"_mht_l1_ptRatio")),
    ('mht_met_dphi', lambda event, syst: getattr(event,syst+"_mht_met_dphi")),
]

