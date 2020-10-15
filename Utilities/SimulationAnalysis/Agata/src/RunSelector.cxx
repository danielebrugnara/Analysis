#include <utility>

#include "RunSelector.h"



RunSelector::RunSelector(TTree& tree, std::string  out_file_name):
                nevts(0),
                the_tree(tree),
                fits(),
                out_file_name(std::move(out_file_name)){}

void RunSelector::Run(const std::vector<double>& energies) {
    nevts = the_tree.GetEntries();
    auto *selector = new Selector();
    auto nentries = TVirtualTreePlayer::kMaxEntries;
    //int nentries = 10;
    auto firstentry = 0;
    the_tree.Process(selector, "option", nentries,firstentry);
    if (nentries<nevts) nevts = nentries;

    TIter iter (selector->GetOutputList());
    TObject *obj;
    auto* out_file = new TFile(out_file_name.c_str(), "recreate");

    //std::string spetrum_name = "addb_spec_DC";
    std::string spetrum_name = "addb_spec_pos_1";
    double eff = 0;
    while ((obj=iter())){
        obj->Write();
        if (spetrum_name==obj->GetName() && !energies.empty()){
            for(const auto& energy: energies) {
                std::vector<std::pair<double, double>> eners;
                eners.emplace_back(energy, 0);
                fits.emplace_back(*dynamic_cast<TH1D *>(obj), eners, false);
                fits.back().EnableCanvas();
                std::vector<Fitter::FitRes> results = fits.back().Fit();
                eff = results.at(0).integral.first / (double) nevts;
                std::cout   << "Efficiency for " << energy
                            <<" keV : " << eff << std::endl;
            }
        }
    }
    TVector2 start_stop(0.,nevts);
    start_stop.Write("start_stop");
    out_file->Close();
    delete selector;
}
