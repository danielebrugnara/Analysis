#include "RunSelector.h"

RunSelector::RunSelector(std::string file_name){

	TFile* file = new TFile(file_name.c_str());
    if (!file) throw std::runtime_error("File not valid\n");
    TTree* tree = (TTree*)file->Get("TreeMaster");
    if (!tree) throw std::runtime_error("Tree not valid\n");

    TVector2 Start_Stop_TS;
    TBranch* AddTS_branch = tree->GetBranch("AddTS");
    long unsigned AddTS{0};
    AddTS_branch->SetAddress(&AddTS);
    int nevents = tree->GetEntries();
    AddTS_branch->GetEntry(0);
    Start_Stop_TS.SetX(static_cast<double>(AddTS)*1.E-8);//In seconds

    AddTS_branch->GetEntry(nevents-1);
    Start_Stop_TS.SetY(static_cast<double>(AddTS)*1.E-8);//In seconds

	auto* selector = new Selector();
	tree->Process(selector, "option");
    std::string file_out_name (file->GetName());

    file->Close();

    file_out_name.insert(0, "spectra_");
    TIter iter (selector->GetOutputList());
    TObject *obj;
    auto* out_file = new TFile(file_out_name.c_str(), "recreate");
    this->file_name = file_out_name;
    while ((obj=iter())){
        obj->Write();
    }
    Start_Stop_TS.Write("start_stop");
    out_file->Close();
}   

std::string RunSelector::GetFileName() const {
    return file_name;
}