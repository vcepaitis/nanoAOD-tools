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
    std::map< std::string, TH1F > pileupHists;
    std::map< std::string, std::map<std::string, double> > processDict;

    if (auto dir = opendir(path.c_str())) {
        while (auto f = readdir(dir)) {
            if (!f->d_name || f->d_name[0] == '.')
                continue; // Skip everything that starts with a dot
            std::string process_string = std::string(f->d_name);
            std::string file_path = path+"/"+process_string;

            std::cout << "Reading in folder: " << file_path << std::endl;

            std::string delimiter = ".";
            std::string process = process_string.substr(0, process_string.find(delimiter));

            if (process.find("Single") != std::string::npos or process.find("EGamma") != std::string::npos) {
                continue;
            }

            std::pair<std::string, double> _yield("unweighted", 0.);
            std::pair<std::string, double> _yieldWeighted("weighted", 0.);
            std::pair<std::string, double> _yieldPositive("posFraction", 0.);
            std::map<std::string, double> processMap;
            processMap.insert(_yield);
            processMap.insert(_yieldWeighted);
            processMap.insert(_yieldPositive);

            double yield(0);
            double yieldWeighted(0);
            double yieldPositive(0);

            pileupHists.insert(std::make_pair(process, TH1F(process.c_str(), "", 101, 0, 100)));
            pileupHists[process].SetDirectory(0);

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
                    tree->Project(h->GetName(),"Pileup_nTrueInt", "genWeight");

                    TH1F* hPos = new TH1F("puPos","",101,0,100);
                    tree->Project(hPos->GetName(),"Pileup_nTrueInt", "(genWeight>0)");

                    yield += h->GetEntries();
                    yieldWeighted += h->Integral();
                    yieldPositive += hPos->Integral();

                    pileupHists[process].Add(h);
                    delete tree;
                    rootFile->Close();
                }
            }
            processMap["unweighted"] = yield;
            processMap["weighted"] = yieldWeighted;
            processMap["posFraction"] = yieldPositive/yield;
            processDict.insert(std::make_pair(process, processMap));
        }
        closedir(dir);
    }
    for (std::pair<std::string, std::map<std::string, double>> x: processDict) {
        std::string processName = x.first;
        std::map<std::string, double> processMap = x.second;
        double yield = processMap["unweighted"];
        double weighted = processMap["weighted"];
        double yieldPositive = processMap["posFraction"];
        //double yield = processMap["yield"];

        std::cout << x.first << " => unweighted: " << yield << ", weighted: " << weighted << ", pos fraction: "   << yieldPositive <<'\n';
    }  

    std::string output_string = (output_path+"/pileup.root");

    TFile *rootFile = TFile::Open(output_string.c_str(), "RECREATE");
    for (std::pair<std::string, TH1F> x: pileupHists) {
        x.second.SetDirectory(rootFile);
        x.second.SetName(x.first.c_str());
        x.second.Write();
    }  

    rootFile->Close();
    json j_map(processDict);
    output_string = output_path+"/eventyields.json";
    std::ofstream o(output_string.c_str());
    o << j_map.dump(0) << std::endl;

    return 0;
}
