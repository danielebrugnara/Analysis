#include "Identification.h"

Identification::Identification() {
}

Identification::~Identification() {
    for (const auto & it_cuts: cuts){
        delete it_cuts.second;
    }
}

void Identification::LoadCuts(std::string path) {
    std::ifstream input_file(path);
    if (!input_file) throw std::runtime_error("File "+path+" not present\n");
    input_file.close();
    TFile *cuts_file = new TFile(path.c_str(), "READ");
    if (!(cuts_file->IsOpen())) throw std::runtime_error("VAMOS file not opened\n");

    TIter contents(cuts_file->GetListOfKeys());
    TKey *key;
    TObject *obj;
#ifdef VERBOSE_DEBUG
    std::cout << "Starting to look for cuts:\n";
#endif
    while ((key = (TKey *)contents())) {
        obj = cuts_file->Get(key->GetName());
        if (obj->InheritsFrom("TCutG")) {
#ifdef VERBOSE_DEBUG
            std::cout << key->GetName() << "\n";
#endif
            TCutG *tmp = (TCutG *)obj;
            cuts[tmp->GetName()] = tmp;
        }
    }
}
