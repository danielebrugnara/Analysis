#include "Analysis.h"

//ClassImp(Analysis);

Analysis::Analysis(){
    std::string file_with_runs = "./Configs/Runs.txt";
    std::ifstream file(file_with_runs);
    if (!file.is_open()) throw std::runtime_error(std::string("Unable to read :")+std::string(file_with_runs));
    std::string line;
    while(std::getline(file, line)){
        std::ifstream file_check(line);
        if (line.at(0)=='#') continue;
        if(file_check.fail()) throw std::runtime_error(std::string("File not present :")+ line);
        file_names.push_back(line);
    }
}

Analysis::~Analysis(){
    delete data;
}

bool Analysis::LoadFiles(){
    data = new TChain("PhysicsTree");
    int n_files=0;
    for (const auto &file_name: file_names){
        n_files = data->Add(file_name.c_str());
    }
    if (n_files!=file_names.size()) throw std::runtime_error(std::string("Not all files were loaded!"));
    else return true;
}

bool Analysis::RunAnalysis(int n_events=-1, int start_event=0){
    bool success_loading = LoadFiles();
    if (!success_loading) return false;
    int analyzed_entries = 0;
    if (n_events!=-1) analyzed_entries = data->Process(&selector, "", n_events, start_event);
    else analyzed_entries = data->Process(&selector, "");
    std::cout << "A total of " << analyzed_entries <<" was processed\n";
    return true;
}