#include <iostream>
#include <string>
#include <fstream>
#include <dirent.h>
#include "TH1.h"
#include "TFile.h"
#include "TTree.h"
#include "json.hpp"

using json = nlohmann::json;

int main(int argc, char **argv){
    std::string path{argv[1]};
    std::string output_path{argv[2]};
    //std::map<std::string, TH1F> pileupHists;
    std::map<std::string, std::map<std::string, double>> processDict;
    //std::map<std::string, TH1F> mapOfHists;

    if (auto dir = opendir(path.c_str())) {
        while (auto f = readdir(dir)) {
            if (!f->d_name || f->d_name[0] == '.')
                continue; // Skip everything that starts with a dot
            std::string process_string = std::string(f->d_name);
            std::string file_path = path+"/"+process_string;


            std::string delimiter = ".";
            std::string process = process_string.substr(0, process_string.find(delimiter));

            if (process.find("HNL") == std::string::npos ) {
                continue;
            }

            std::cout << "Reading in folder: " << file_path << std::endl;


            //pileupHists.insert(std::make_pair(process, TH1F(process.c_str(), "", 101, 0, 100)));
            //pileupHists[process].SetDirectory(0);
            for (int weight_index=2; weight_index<68; weight_index++){
                std::map<std::string, double> dummy_map;
                std::string name_string = "LHEWeights_coupling_"+std::to_string(weight_index);
                dummy_map.insert(std::make_pair(name_string, 0.));
                processDict.insert(std::make_pair(process, dummy_map));
            }

            std::ifstream file(file_path);
            if (file.is_open()) {
                std::string line;
                while (std::getline(file, line)) {
                    std::cout << "Reading in file: " << line << std::endl;
                    // Open with root
                    TFile *rootFile = TFile::Open(line.c_str());
                    if (rootFile->IsZombie()) {
                        std::cout << "Error opening file" << std::endl;
                        continue;
                    }
                    TTree* tree = (TTree*)rootFile->Get("Events");
                    TH1F* h = new TH1F("pu","",101,0,100);
                    for (int weight_index=2; weight_index<68; weight_index++){
                        std::string weight_string = "genWeight*LHEWeights_coupling_"+std::to_string(weight_index);
                        std::string name_string = "LHEWeights_coupling_"+std::to_string(weight_index);
                        tree->Project(h->GetName(), "Pileup_nTrueInt", weight_string.c_str());
                        processDict[process][name_string] += h->Integral();
                    }
                    //pileupHists[process].Add(h);
                    delete tree;
                    rootFile->Close();
                }
            }
        }
        closedir(dir);
    }

    /*
    for (std::pair<std::string, std::pair<int, double>> x: processDict) {
        std::cout << x.first << " => " << x.second << '\n';
    }
    */

    std::string output_string = (output_path+"/pileup.root");

    /*
    TFile *rootFile = TFile::Open(output_string.c_str(), "RECREATE");
    for (std::pair<std::string, TH1F> x: pileupHists) {
        x.second.SetDirectory(rootFile);
        x.second.SetName(x.first.c_str());
        x.second.Write();
    }

    rootFile->Close();
    */
    json j_map(processDict);

    output_string = output_path+"/eventyields.json";
    std::ofstream o(output_string.c_str());
    o << j_map.dump(0) << std::endl;

    return 0;
}
