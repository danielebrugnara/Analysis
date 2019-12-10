#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <thread>
#include <stack>
#include <mutex>

//Root headers
//#include <TChain.h>
#include <TFile.h>
#include <TTree.h>
#include <Selector.h>


class Analysis{
    public:        
        Analysis(int);
        ~Analysis();

        bool RunAnalysis(int, int);
        bool RunSelector(std::string);
    private:
        //bool LoadFiles();
        Selector selector;
        int n_threads;
        //TChain * data;
        std::stack<std::string> file_names;
        std::vector<std::thread> threads;
        std::mutex mtx;
        std::string GetRun();
        bool Job();
        //ClassDef(Analysis, 1);
};