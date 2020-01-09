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
    bool ReadRunsFromFile();
    int n_threads;
    std::stack<std::string> file_names;
    std::stack<std::string> processed_files;
    std::vector<std::thread> threads;
    std::mutex mtx;
    std::string GetRun();
    bool Job();
    //ClassDef(Analysis, 1);
};

#endif
