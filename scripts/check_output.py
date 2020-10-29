import json
import ROOT
import argparse
import os
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
    try:
        status = _f.IsZombie()
    except:
        ReferenceError
        return 1
    else:
        _f.Close()
        if status:
            return 1
        return 0


bad_files = []
with open(input_file) as f:
    for l in f:
        l = l.rstrip()
        proc_name = l.split("/")[7].replace('.txt', '')
        output_path_proc = os.path.join(output_path, proc_name)
        with open(l) as subf:
            for root_file in subf:
                root_file_friend = root_file.rstrip().split("/")[-1].replace('.root', '_Friend.root')
                root_file_output_path = os.path.join(output_path_proc, root_file_friend)
                if check_root_file(root_file_output_path):
                    bad_files.append(root_file_output_path)

with open('batch/resubmit/{}'.format(resubmit_file), 'w') as f:
    for bad_file in bad_files:
        f.write(bad_file+"\n")