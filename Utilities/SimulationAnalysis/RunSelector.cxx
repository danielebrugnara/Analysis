#include "RunSelector.h"

RunSelector::RunSelector(std::string file_name){
	TFile* file = new TFile(file_name.c_str());
	TTree* tree = (TTree*) file->Get("SimulatedTree");
	Selector* selector = new Selector();
	std::cout << "Starting selector\n";
	tree->Process(selector, "option");
}
