#include "Analysis.h"

//ClassImp(Analysis);

Analysis::Analysis(int n_threads) : n_threads(n_threads) {
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

    std::ifstream file_test("./Configs/Interpolations/TW_Ice_Thickness.root");
    generate_TW_ice_interpolation = file_test ? false : true;
}

Analysis::~Analysis() {
}

bool Analysis::RunAnalysis() {
    //Removing previous data
    system("rm -rf ./Out");
    system("mkdir ./Out");

    //Assigning to threads
    for (int ii = 0; ii < n_threads; ++ii) {
        threads.push_back(std::thread(&Analysis::Job, this));
    }

    //Waiting for threads to finish
    for (int ii = 0; ii < n_threads; ++ii) {
        threads.at(ii).join();
    }

    std::cout << "Joined threads, adding output files :\n";

    //Adding all data in unique file
    system("hadd -f ./Out/sum.root ./Out/*");
    std::ofstream output_file("./Configs/Interpolations/TW_Ice_Thickness.txt");

    if (generate_TW_ice_interpolation) {
        std::sort(data.TW_vs_ice.begin(), data.TW_vs_ice.end());
        TFile *root_file = new TFile("./Configs/Interpolations/TW_Ice_Thickness.root",
                                     "recreate");

        double X[data.TW_vs_ice.size()];
        double Y[data.TW_vs_ice.size()];
        for (int ii = 0; data.TW_vs_ice.size(); ++ii) {
            X[ii] = data.TW_vs_ice[ii].first;
            Y[ii] = data.TW_vs_ice[ii].second;
            output_file << X[ii] << "\t" << Y[ii] << "\n";
        }
        output_file.close();
        TGraph *gr = new TGraph(data.TW_vs_ice.size(),X, Y);
        TSpline3 *spl = new TSpline3("TW_vs_ice_thickness", gr);
        spl->Write();
        root_file->Write();
        root_file->Close();
    }

    return true;
}

//This job will be assigned to a thread
bool Analysis::Job() {
    try {
        while (1) {  //Repeat while there are jobs
            std::string current_run = GetRun();
            std::cout << "Run : " << current_run << " assigned to thread\n";
            pid_t pid;
            int status;
            //Forking process to avoid ROOT threading problems
            if ((pid = fork()) == 0) {  //Child process
                RunSelector(current_run);
                exit(1);
            } else {            //Parent process
                wait(&status);  //Wait child process to complete instructions
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
    std::cout << "---------->Starting selector for run : " << run << "<----------\n";
    Selector *selector = new Selector();
    tree->Process(selector, ("analyzed_" + run.substr(run.find_last_of("/") + 1)).c_str());
    if (generate_TW_ice_interpolation) {
        std::vector<std::pair<double, double>> tmp_vec = selector->GetTWvsIce(); 
        mtx_data.lock();
        data.TW_vs_ice.insert(data.TW_vs_ice.end(),
                              tmp_vec.begin(),
                              tmp_vec.end());
        mtx_data.unlock();
    }
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
