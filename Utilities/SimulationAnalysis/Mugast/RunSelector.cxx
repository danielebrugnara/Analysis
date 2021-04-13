#include "RunSelector.h"

RunSelector::RunSelector(std::string file_name){
	TFile* file = new TFile(file_name.c_str());
  if (file == nullptr) std::cout << "File not found";
	TTree* tree = (TTree*) file->Get("PhysicsTree");
  if (tree == nullptr) std::cout << "Tree not found";
	Selector* selector = new Selector();
	std::cout << "Starting selector\n";
	tree->Process(selector, "option");
}
