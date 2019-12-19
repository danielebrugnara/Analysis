#ifndef __ANALYSIS_H__
#define __ANALYSIS_H__

#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <thread>
#include <stack>
#include <mutex>

#include <sys/wait.h>

//Root headers
//#include <TChain.h>
#include <TFile.h>
#include <TTree.h>
#include <Selector.h>


class Analysis{
    public:        
        Analysis(int);
        ~Analysis();

        bool RunAnalysis();
        bool RunSelector(std::string);
    private:
        //bool LoadFiles();
        int n_threads;
        //std::map<std::string, Selector> selector;
        //TChain * data;
        std::stack<std::string> file_names;
        std::vector<std::thread> threads;
        std::mutex mtx;
        std::string GetRun();
        bool Job();
        //ClassDef(Analysis, 1);
};

#endif