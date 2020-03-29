#ifndef __ANALYSIS_H__
#define __ANALYSIS_H__

#include <sys/wait.h>

#include <fstream>
#include <iostream>
#include <mutex>
#include <stack>
#include <string>
#include <thread>
#include <vector>

//Root headers
#include <Selector.h>
#include <TFile.h>
#include <TTree.h>

class Analysis {
   private:
    static constexpr int data_size = 1000;
    struct Data_partial {
        double TW_vs_ice[data_size][2];
        Data_partial(std::vector<std::pair<double, double>> const & TW_vs_ice){
            if (TW_vs_ice.size()> data_size)
                throw std::runtime_error("data_size in analysis class too small to contain data");
            long unsigned int ii=0;
            for (;ii < TW_vs_ice.size(); ++ii){
                this->TW_vs_ice[ii][0]=TW_vs_ice[ii].first;
                this->TW_vs_ice[ii][1]=TW_vs_ice[ii].second;
            }
            for (;ii<data_size; ++ii){
                this->TW_vs_ice[ii][0] = 0;
                this->TW_vs_ice[ii][1] = 0;
            }
        };
        Data_partial(){};
    };
    struct Data {
        std::vector<std::pair<double, double>> TW_vs_ice;
        //Data(std::vector<std::pair<double, double>> const & TW_vs_ice):
        //        TW_vs_ice(TW_vs_ice){};
        Data(){};
    } data;

   public:
    Analysis(int);
    ~Analysis();

    bool RunAnalysis();
    Data_partial * RunSelector(std::string);

   private:
    int n_threads;
    std::stack<std::string> file_names;
    std::stack<std::string> processed_files;
    std::vector<std::thread> threads;
    bool generate_TW_ice_interpolation;
    std::mutex mtx_job_que;
    std::mutex mtx_data;
    std::string GetRun();
    bool Job();

    void UpdateData(Data_partial & );
    //ClassDef(Analysis, 1);
};

#endif
