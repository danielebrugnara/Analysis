#include "../../Lib/inc/MugastData.h"
#include <TTreeReader.h>
#include <TTreeReaderValue.h>


void CreateSimulationFile(std::map<std::string, double> files, const std::string& outputName="simulationsum.root"){
    std::cout << "Output name : " << outputName << std::endl;
    TFile* outFile = new TFile(outputName.c_str(), "recreate");
    TTree* outTree = new TTree("AnalyzedTree", "AnalyzedTree");
    MugastData outData(1,0);
    MugastData* inData{0};
    outTree->Branch("MugastData", &outData);

    for(const auto& it: files){
        TFile* inFile = new TFile(it.first.c_str(), "read");
        if (inFile == nullptr) throw std::runtime_error("file not found: "+it.first);
        TTree* inTree = (TTree*) inFile->Get("AnalyzedTree");
        if (inTree == nullptr) throw std::runtime_error("Tree  not found in: "+it.first);
        TTreeReader reader;
        reader.SetTree(inTree);
        TTreeReaderValue<MugastData> inData{reader, "MugastData"};

        for (int i=0; i<inTree->GetEntries()*it.second; ++i){
           reader.SetEntry(i);

            outData.~MugastData();
            new (&outData) MugastData(inData->multiplicity, 0);

            for (int j{0}; j<inData->multiplicity; ++j) {
                outData.MG[j] = inData->MG[j];
                outData.Pos[j] = inData->Pos[j];
                outData.EmissionDirection[j] = inData->EmissionDirection[j];
                outData.SI_X[j] = inData->SI_X[j];
                outData.SI_Y[j] = inData->SI_Y[j];
                outData.SI_E[j] = inData->SI_E[j];
                outData.E[j] = inData->E[j];
                outData.Ex[j] = inData->Ex[j];
                outData.Theta_CM[j] = inData->Theta_CM[j];
            }

           outTree->Fill();

        }
        
    }
    outFile->cd();
    outTree->Write();
    outFile->Write();
    outFile->Close();
}
