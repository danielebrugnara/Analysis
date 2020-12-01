#include "Identification.h"

#ifdef VERBOSE_DEBUG
#  define DEBUG(x, y) std::cout << (x) << (y) << std::endl;
#else
#  define DEBUG(x, y) 
#endif

Identification::Identification()
= default;

Identification::~Identification()
{
    for (const auto &it_cuts : cuts)
    {
        delete it_cuts.second;
    }
}

void Identification::LoadCuts(const std::string& path)
{
    std::ifstream input_file(path);
    if (!input_file)
        throw std::runtime_error("File " + path + " not present\n");
    input_file.close();
    auto *cuts_file = new TFile(path.c_str(), "READ");
    if (!(cuts_file->IsOpen()))
        throw std::runtime_error("VAMOS file not opened\n");

    TIter contents(cuts_file->GetListOfKeys());
    TKey *key;
    TObject *obj;
    DEBUG("Starting to look for cuts:", "");
    while ((key = (TKey *)contents()))
    {
        obj = cuts_file->Get(key->GetName());
        if (obj->InheritsFrom("TCutG"))
        {
            DEBUG("Key : ", key->GetName());
            auto *tmp = (TCutG *)obj;
            cuts[tmp->GetName()] = tmp;
        }
    }
    DEBUG("Finished looking for cuts", "");
    }