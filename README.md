# nanoAOD-tools
Tools for working with NanoAOD for the HNL legacy hadronic analysis

![build tests](https://travis-ci.org/LLPDNNX/nanoAOD-tools.svg?branch=HNL)

## Checkout instructions:

Using the particular CMSSW release which is shipped with all the required software (especially tensorflow C++ API v1.6). To be run from the imperial cluster lx0[23].hep.ph.ic.ac.uk:

```
cmsrel CMSSW_10_2_18
cd CMSSW_10_2_18/src
cmsenv
git clone -b HNL https://github.com/LLPDNNX/nanoAOD-tools.git PhysicsTools/NanoAODTools
scram b
```

You will likely need a proxy to access the files stored on the grid:

```
voms-proxy-init --voms cms
```

Example use: make a single ntuple with the HNL framework:

```
python PhysicsTools/NanoAODTools/processors/HNL.py --input root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mkomm/HNL/NANOX_201117/HNL_dirac_all_ctau1p0e00_massHNL10p0_Vall1p664e-03-2016/HNL_dirac_all_ctau1p0e00_massHNL10p0_Vall1p664e-03/HNL_dirac_all_ctau1p0e00_massHNL10p0_Vall1p664e-03-2016/201118_145719/0000/nano_9.root .
```

Running on data instead:

```
python PhysicsTools/NanoAODTools/processors/HNL.py --isData --input root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mkomm/HNL/NANOX_201117/SingleMuon_Run2016H/SingleMuon/SingleMuon_Run2016H/201121_114955/0001/nano_1055.root .
```

For more instructions refer to the twiki.

