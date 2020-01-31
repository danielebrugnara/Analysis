#ifndef __ANALYSIS_H__
#define __ANALYSIS_H__

#include <fstream>
#include <iostream>
#include <mutex>
#include <stack>
#include <string>
#include <thread>
#include <vector>

#include <sys/wait.h>

//Root headers
#include <Selector.h>
#include <TFile.h>
#include <TTree.h>

class Analysis {
   public:
    Analysis(int);
    ~Analysis();

    bool RunAnalysis();
    bool RunSelector(std::string);

   private:
    int n_threads;
    std::stack<std::string> file_names;
    std::stack<std::string> processed_files;
    std::vector<std::thread> threads;
    bool generate_TW_ice_interpolation;
    std::mutex mtx;
    std::mutex mtx_data;
    std::string GetRun();
    bool Job();

    struct Data{
        std::vector<std::pair<double, double>> TW_vs_ice;
    } data;
    //ClassDef(Analysis, 1);
};

#endif
