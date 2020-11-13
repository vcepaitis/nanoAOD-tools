import json
import ROOT
import argparse
import os
import uproot
parser = argparse.ArgumentParser()

parser.add_argument('--input')
parser.add_argument('--output')
parser.add_argument('--file')

args = parser.parse_args()
input_file = args.input
output_path = args.output
resubmit_file = args.file


def check_root_file(f):
    _f = ROOT.TFile.Open(f)
    if not(_f):
        return 1
    try:
        status = _f.IsZombie()
    except:
        ReferenceError
        _f.Close()
        return 1
    else:
        if _f.TestBit(ROOT.TFile.kRecovered):
            _f.Close()
            return 1
        if status:
            _f.Close()
            return 1
        tree = _f.Get("Friends")
        if not tree:
            _f.Close()
            return 1
        _f.Close()
        return 0

def check_uproot_file(f):
    try: 
        _f = uproot.open(f)
    except ValueError:
        print("Broken file: {}".format(f))
        return 1
    else:
        if not _f:
            return 1
        tree = _f["Friends"]
        if not tree:
            return 1
        return 0

bad_files = []
with open(input_file) as f:
    for l in f:
        l = l.rstrip()
        proc_name = l.split("/")[6].replace('.txt', '')
        print(proc_name)
        output_path_proc = os.path.join(output_path, proc_name)
        with open(l) as subf:
            for root_file in subf:
                root_file_friend = root_file.rstrip().split("/")[-1].replace('.root', '_Friend.root')
                root_file_output_path = os.path.join(output_path_proc, root_file_friend)
                if check_root_file(root_file_output_path):
                    bad_files.append(root_file.rstrip())

with open('batch/resubmit/{}'.format(resubmit_file), 'w') as f:
    for bad_file in bad_files:
        f.write(bad_file+"\n")