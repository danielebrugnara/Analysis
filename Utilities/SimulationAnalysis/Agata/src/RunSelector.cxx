#include "RunSelector.h"
#include "Readsimu.h"

RunSelector::RunSelector(std::string file_name):nevts(0){
    if (file_name.compare(file_name.size()-5, 5, ".root")){
    //is not a root file, generate it!
        bool smearing = true;
        int it_max = -1;
        std::cout << "Calling readsimu\n";
        nevts = readsimu(file_name, file_name+".root", smearing, it_max);
        file_name = file_name+".root";
    }
    std::cout << "nevts read by readimu: " << nevts << std::endl;
	auto* file = new TFile(file_name.c_str());
    if (!file) throw std::runtime_error("File not valid\n");
    auto* tree = (TTree*)file->Get("SimulatedAgata");
    if (!tree) throw std::runtime_error("Tree not valid\n");
	auto* selector = new Selector();
	tree->Process(selector, "option");
    file_out_name = file->GetName();
    file->Close();

    file_out_name.insert(0, "spectra_");
    TIter iter (selector->GetOutputList());
    TObject *obj;
    auto* out_file = new TFile(file_out_name.c_str(), "recreate");
    while ((obj=iter())){
        obj->Write();
    }
    out_file->Close();
    delete selector;
}   