import uproot
import pandas as pd
import numpy as np
from xgboost import XGBClassifier
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score
from sklearn import metrics
from matplotlib import pyplot
import pickle

def get_df(path, branches, sig):
	iters = []
	for data in uproot.pandas.iterate(path, "Friends", branches=branches, flatten=False):
		#print(path)

		data = data[data.ntightMuon == 1]
		data = data[data.nlepJet_nominal == 1]
		data = data[data.nlooseMuons == 1]
		data = data[data.lepJet_nominal_deltaR < 0.4]
		data = data[data.dilepton_charge == -1]
		data = data[data.IsoMuTrigger_flag == True]

		if sig == True:
			data = data[data.LHEWeights_coupling_12 > 0]

		iters.append(data)
		# add cuts

	return pd.concat(iters)

'''
array_list = [
			  "dimuon_mass", "dimuon_deltaR",
			  "lepJet_pt", "lepJet_eta", "lepJet_deltaR",
			  "MET_pt", "MET_phi",
			  "looseMuons_pt", "looseMuons_eta", "looseMuons_dxy",
			  "tightMuons_pt", "tightMuons_eta", "tightMuons_dxy"
			 ]


array_list = [
			"EventObservables_nominal_met", "EventObservables_nominal_mht_NoMu",
			"EventObservables_nominal_mht", "EventObservables_nominal_ht",
			"EventObservables_nominal_MET_Mu_mT", "selectedJets_nominal_pt[0]",
			"selectedJets_nominal_eta[0]", "selectedJets_nominal_phi[0]",
			"nselectedJets_nominal", "dilepton_mass", "dilepton_pt", "dilepton_eta",
			"looseMuons_pt[0]", "looseMuons_eta[0]", "looseMuons_phi[0]",
			"tightMuon_pt[0]", "tightMuon_eta[0]", "tightMuon_phi[0]",
			"ntightMuon", "nlepJet_nominal", "nlooseMuons", "lepJet_nominal_deltaR",
			"dilepton_charge", "IsoMuTrigger_flag", "LHEWeights_coupling_12"
			]
'''
array_preselection = [
			"ntightMuon", "nlepJet_nominal", "nlooseMuons", "lepJet_nominal_deltaR",
			"dilepton_charge", "IsoMuTrigger_flag"
			]

array_bdt = [
			"EventObservables_nominal_met", "EventObservables_nominal_met_phi",
			"EventObservables_nominal_ht", "EventObservables_nominal_mT_met_Mu",
			"selectedJets_nominal_pt", "selectedJets_nominal_eta",
			"nselectedJets_nominal", "dilepton_mass", "dilepton_deltaPhi",
			"dilepton_deltaR", "tightMuon_pt", "tightMuon_eta", "tightMuon_phi",
			"looseMuons_pt", "looseMuons_eta"
			]

array_list = array_preselection + array_bdt

n_events=150000

sig_df = get_df("/vols/cms/jd918/LLP/CMSSW_10_2_18/src/nanoAOD_friends_200622/2016//HNL_*/nano_*_Friend.root", array_list + ["LHEWeights_coupling_12"], True)

wjets_df = get_df("/vols/cms/jd918/LLP/CMSSW_10_2_18/src/nanoAOD_friends_200622/2016//WToLNu_*/nano_*_Friend.root", array_list, False).sample(n=n_events)
tt_df = get_df("/vols/cms/jd918/LLP/CMSSW_10_2_18/src/nanoAOD_friends_200622/2016//TTToSemiLeptonic_*/nano_*_Friend.root", array_list, False).sample(n=n_events)
dyjets_df = get_df("/vols/cms/jd918/LLP/CMSSW_10_2_18/src/nanoAOD_friends_200622/2016//DYJetsToLL*amcatnlo*/nano_*1_Friend.root", array_list, False).sample(n=n_events)

'''
sig_df = get_df("/vols/cms/jd918/LLP/CMSSW_10_2_18/src/nanoAOD_friends_200622/2016//HNL_dirac_all_ctau1p0e01_massHNL4p5_Vall4p549e-03-2016/nano_1_Friend.root", array_list + ["LHEWeights_coupling_12"], True)
wjets_df = get_df("/vols/cms/jd918/LLP/CMSSW_10_2_18/src/nanoAOD_friends_200622/2016//WToLNu_0J_13TeV-amcatnloFXFX-pythia8-2016/nano_1_Friend.root", array_list, False)
tt_df = get_df("/vols/cms/jd918/LLP/CMSSW_10_2_18/src/nanoAOD_friends_200622/2016//TTToSemiLeptonic_*/nano_1_Friend.root", array_list, False)
dyjets_df = get_df("/vols/cms/jd918/LLP/CMSSW_10_2_18/src/nanoAOD_friends_200622/2016//DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8-2016/nano_1_Friend.root", array_list, False)
#bkg_df = pd.concat([wjets_df, dyjets_df, tt_df]).sample(n=n_events)
'''
bkg_df = pd.concat([wjets_df, dyjets_df, tt_df])
print(bkg_df.shape[0])

sig_label = np.ones(sig_df.shape[0])
bkg_label = np.zeros(bkg_df.shape[0])

print("sig_df")
print(list(sig_df))

print("bkg_df")
print(list(bkg_df))

df = pd.concat([sig_df[array_bdt], bkg_df[array_bdt]])
label = np.concatenate([sig_label, bkg_label])

feature_array = [
				'tightMuon_pt', 'tightMuon_eta', 'tightMuon_phi',
				'looseMuons_pt', 'looseMuons_eta',
				'selectedJets_nominal_pt', 'selectedJets_nominal_eta'
				]

for feat in feature_array:
	df[feat] = df[feat].map(lambda x: x[0])

df["EventObservables_nominal_met_phi-tightMuon_phi"] = df["EventObservables_nominal_met_phi"] - df["tightMuon_phi"]
df.drop("tightMuon_phi", axis=1, inplace=True)
print(df)


X_train, X_test, y_train, y_test = train_test_split(df, label, test_size=0.2, random_state=42, stratify=label)


# fit model no training data
model = XGBClassifier(objective='binary:logistic', learning_rate=0.1, max_depth=10, n_estimators=200, nthread=-1)
eval_set = [(X_train, y_train), (X_test, y_test)]
model.fit(X_train, y_train, eval_metric=["error", "logloss"], eval_set=eval_set, verbose=True)
# make predictions for test data
y_pred = model.predict(X_test)
predictions = [round(value) for value in y_pred]
# evaluate predictions
accuracy = accuracy_score(y_test, predictions)
print("Accuracy: %.2f%%" % (accuracy * 100.0))

# retrieve performance metrics
results = model.evals_result()
epochs = len(results['validation_0']['error'])
x_axis = range(0, epochs)

# plot log loss
fig, ax = pyplot.subplots()
ax.plot(x_axis, results['validation_0']['logloss'], label='Train')
ax.plot(x_axis, results['validation_1']['logloss'], label='Test')
ax.legend()
pyplot.ylabel('Log Loss')
pyplot.title('XGBoost Log Loss')
#pyplot.show()
pyplot.savefig('BDT_loss.pdf')
pyplot.savefig('BDT_loss.png')

# plot classification error
fig, ax = pyplot.subplots()
ax.plot(x_axis, results['validation_0']['error'], label='Train')
ax.plot(x_axis, results['validation_1']['error'], label='Test')
ax.legend()
pyplot.ylabel('Classification Error')
pyplot.title('XGBoost Classification Error')
#pyplot.show()
pyplot.savefig('BDT_error.pdf')
pyplot.savefig('BDT_error.png')

# take the second column because the classifier outputs scores for
# the 0 class as well
probs = model.predict_proba(X_test)[:, 1]

fig, axes = pyplot.subplots()
pyplot.hist(probs[y_test==0], label='background')
pyplot.hist(probs[y_test==1], label='signal')
axes.legend()
pyplot.ylabel('BDT output')
pyplot.xlabel('signal (1) vs background (0)')
pyplot.title('XGBoost separation')
pyplot.savefig('BDT_output.pdf')
pyplot.savefig('BDT_output.png')

# fpr means false-positive-rate
# tpr means true-positive-rate
fpr, tpr, _ = metrics.roc_curve(y_test, probs)

auc_score = metrics.auc(fpr, tpr)
print('AUC = {:.3f}'.format(auc_score))

fig, ax = pyplot.subplots()
ax.plot(tpr, fpr, label='AUC = {:.3f}'.format(auc_score))
ax.legend(loc='lower right')
pyplot.yscale('log')
pyplot.xlabel('Signal Efficiency')
pyplot.ylabel('Background Efficiency')
pyplot.title('XGBoost ROC curve')
#pyplot.show()
pyplot.savefig('BDT_roc.pdf')
pyplot.savefig('BDT_roc.png')

model._Booster.save_model("/vols/cms/jd918/LLP/CMSSW_10_2_18/src/PhysicsTools/NanoAODTools/data/bdt/bdt.model")
#pyplot.show()
