function run_test()
{
    cd ~
    source /cvmfs/cms.cern.ch/cmsset_default.sh
    export SCRAM_ARCH=slc6_amd64_gcc700 || return 1
    scramv1 project CMSSW CMSSW_10_2_18 || return 1
    cd CMSSW_10_2_18/src || return 1
    eval `scram runtime -sh` || return 1
    rsync -r /scripts/* LLPReco  || return 1
    scram b || return 1
    mkdir -p PhysicsTools/NanoAODTools
    rsync -r --stats /scripts/ PhysicsTools/NanoAODTools/. || return 1
    scram b || return 1
    echo "--- Test HNL script ---"
    python PhysicsTools/NanoAODTools/processors/HNL.py --testMode --input=https://github.com/LLPDNNX/test-files/raw/master/nanoaod/Moriond17_aug2018_miniAODv3_HNL_nanoAODv2p1.root . || return 1
}

run_test
