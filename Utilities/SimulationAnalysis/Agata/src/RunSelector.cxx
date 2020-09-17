#include "RunSelector.h"

#include "Selector.h"
#include "Readsimu.h"


RunSelector::RunSelector(std::string file_name):nevts(0){
    bool run_readsimu = false;
    if (file_name.compare(file_name.size()-5, 5, ".root")) {
        run_readsimu = true;
    }

    if (run_readsimu){
    //is not a root file, generate it!
        bool smearing = true;
        int it_max = -1;
        std::cout << "Converting simulation\n";
        nevts = readsimu(file_name, file_name+".root", smearing, it_max);
        file_name = file_name+".root";
        std::cout << "\nnevts read by readimu: " << nevts << std::endl;
    }

	auto* file = new TFile(file_name.c_str());
    if (! file->IsOpen()) throw std::runtime_error("File not valid\n");
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
    TVector2 start_stop(0.,nevts);
    start_stop.Write("start_stop");
    out_file->Close();
    delete selector;
}   
