#include "Identification.h"

Identification::Identification() {}

Identification::~Identification() {}

void Identification::LoadCuts(std::string path) {
    TFile *cuts_file = new TFile(path.c_str(), "READ");
    if (!(cuts_file->IsOpen())) throw std::runtime_error("VAMOS file not opened\n");

    TIter contents(cuts_file->GetListOfKeys());
    TKey *key;
    TObject *obj;
    while ((key = (TKey *)contents())) {
        obj = cuts_file->Get(key->GetName());
        if (obj->InheritsFrom("TCutG")) {
            TCutG *tmp = (TCutG *)obj;
            cuts[tmp->GetName()] = tmp;
        }
    }
}
