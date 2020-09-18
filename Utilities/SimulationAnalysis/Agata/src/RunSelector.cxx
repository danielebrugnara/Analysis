#include "RunSelector.h"

#include "Selector.h"


RunSelector::RunSelector(TTree& tree, const std::string& out_file_name):nevts(0){
    nevts = tree.GetEntries();
	auto *selector = new Selector();
	auto nentries = TVirtualTreePlayer::kMaxEntries;
	//int nentries = 10;
	auto firstentry = 0;
	tree.Process(selector, "option", nentries,firstentry);
	if (nentries<nevts) nevts = nentries;

    TIter iter (selector->GetOutputList());
    TObject *obj;
    auto* out_file = new TFile(out_file_name.c_str(), "recreate");
    while ((obj=iter())){
        obj->Write();
    }
    TVector2 start_stop(0.,nevts);
    start_stop.Write("start_stop");
    out_file->Close();
    delete selector;
}
