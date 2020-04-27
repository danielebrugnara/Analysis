#include "RunSelector.h"

RunSelector::RunSelector(std::string file_name){
	TFile* file = new TFile(file_name.c_str());
	TTree* tree = (TTree*) file->Get("SimulatedTree");
	std::cout << "Starting selector\n";
	Selector* selector = new Selector();
	tree->Process(selector);
}
