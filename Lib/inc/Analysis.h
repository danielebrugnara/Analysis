#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <thread>

#include <TChain.h>

#include "Selector.h"


class Analysis{
    public:        
        Analysis();
        ~Analysis();

        bool RunAnalysis(int, int);
    private:
        bool LoadFiles();
        Selector selector;
        TChain * data;
        std::vector<std::string> file_names;
        //ClassDef(Analysis, 1);
};