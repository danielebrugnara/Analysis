#include "Analysis.h"

//ClassImp(Analysis);

Analysis::Analysis(int n_threads) : n_threads(n_threads) {
    ReadRunsFromFile();
}

bool Analysis::ReadRunsFromFile(){
    std::string file_with_runs = "./Configs/Runs.txt";
    std::ifstream file(file_with_runs);
    if (!file.is_open()) throw std::runtime_error(std::string("Unable to read :") + std::string(file_with_runs));
    std::string line;
    while (std::getline(file, line)) {
        std::ifstream file_check(line);
        if (line.at(0) == '#') continue;
        if (file_check.fail()) throw std::runtime_error(std::string("File not present :") + line);
        file_names.push(line);
    }

}

Analysis::~Analysis() {
}

bool Analysis::RunAnalysis() {
    system("rm -rf ./Out");
    system("mkdir ./Out");
    for (int ii = 0; ii < n_threads; ++ii) {
        threads.push_back(std::thread(&Analysis::Job, this));
    }
    for (int ii = 0; ii < n_threads; ++ii) {
        threads.at(ii).join();
    }

    std::cout << "Joined threads, adding output files :\n";

    system("hadd -f ./Out/sum.root ./Out/*");

    return true;
}

bool Analysis::Job() {
    try {
        while (1) {
            std::string current_run = GetRun();
            std::cout << "Run : " << current_run << " assigned to thread\n";
            pid_t pid;
            int status;
            //Forking process to avoid ROOT threading problems
            if ((pid = fork()) == 0) {  //Child process
                RunSelector(current_run);
                //std::cout << "Finished analysis  of : "<< current_run <<"\n";
                exit(1);
            } else {            //Parent process
                wait(&status);  //Wait child process to complete instructions
                //std::cout << "Parent finished waiting for child process  : " << current_run <<"\n";
            }
        }
    } catch (std::runtime_error &e) {
        std::cout << e.what();
    }
    return false;
}

bool Analysis::RunSelector(std::string run) {
    TFile *file = new TFile(run.c_str());
    if (!file) throw std::runtime_error("Run : " + run + " Not opened\n");
    TTree *tree = (TTree *)file->Get("PhysicsTree");
    std::cout << "Starting selector for run : " << run << "\n";
    Selector *selector = new Selector();
    tree->Process(selector, ("analyzed_" + run.substr(run.find_last_of("/") + 1)).c_str());
    delete selector;
    return true;
}

std::string Analysis::GetRun() {
    std::string run;
    mtx.lock();
    if (file_names.empty()) {
        mtx.unlock();
        throw std::runtime_error("Finished Runs\n");
    }
    run = file_names.top();
    file_names.pop();
    mtx.unlock();
    return run;
}
