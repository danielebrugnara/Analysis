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
    //Creating needed files
    system("cd Configs/Cuts/; root -l -b -q E_TOF_cuts.cxx; cd -;");
    system("cd Configs/Cuts/; root -l -b -q MQ_Q_cuts.cxx; cd -;");
    system("cd Configs/Interpolations/; root -l -b -q TW_Brho_M46_Z18.C; cd -;");


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

    if (generate_TW_ice_interpolation) {
        std::sort(data.TW_vs_ice.begin(), data.TW_vs_ice.end());
        TFile *root_file = new TFile("./Configs/Interpolations/TW_Ice_Thickness.root",
                                     "recreate");

        double X[data.TW_vs_ice.size()];
        double Y[data.TW_vs_ice.size()];
        for (long unsigned int ii = 0; ii < data.TW_vs_ice.size(); ++ii) {
            X[ii] = data.TW_vs_ice[ii].first;
            Y[ii] = data.TW_vs_ice[ii].second;
        }
        TGraph *gr = new TGraph(data.TW_vs_ice.size(), X, Y);
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

            int p[2];

            if (pipe(p) < 0)
                throw std::runtime_error("Pipe error");

            //Forking process to avoid ROOT threading problems
            Data_partial *partial_data = nullptr;

            switch (pid = fork()) {
                case -1:
                    throw std::runtime_error("Fork error");
                    break;
                case 0:  //Child process
                    close(p[0]);
                    if (generate_TW_ice_interpolation) {
                        partial_data = RunSelector(current_run);
                        FILE *out = fdopen(p[1], "w");
                        fwrite(partial_data, sizeof(*partial_data), 1, out);
                        fclose(out);
                    } else {
                        RunSelector(current_run);
                    }
                    exit(1);
                    break;
                default:
                    close(p[1]);
                    if (generate_TW_ice_interpolation) {
                        partial_data = new Data_partial();
                        FILE *in = fdopen(p[0], "r");
                        fread(partial_data, sizeof(*partial_data), 1, in);
                        fclose(in);
                    }
                    wait(&status);  //Wait child process to complete instructions
            }
            if (generate_TW_ice_interpolation) {
                UpdateData(*partial_data);
                if (partial_data) delete partial_data;
            }
        }
    } catch (std::runtime_error &e) {
        std::cout << e.what();
    }
    return false;
}

Analysis::Data_partial *Analysis::RunSelector(std::string run) {
    TFile *file = new TFile(run.c_str());
    if (!file) throw std::runtime_error("Run : " + run + " Not opened\n");
    TTree *tree = (TTree *)file->Get("PhysicsTree");
    std::cout << "---------->Starting selector for run : " << run << "<----------\n";
    Selector *selector = new Selector();
    tree->Process(selector, ("analyzed_" + run.substr(run.find_last_of("/") + 1)).c_str(), 2000);
    if (generate_TW_ice_interpolation) {
        Data_partial *partial_data = new Data_partial(selector->GetTWvsIce());
        delete selector;
        return partial_data;
    } else {
        delete selector;
        Data_partial *dummy_data = nullptr;
        return dummy_data;
    }
}

std::string Analysis::GetRun() {
    std::string run;
    mtx_job_que.lock();
    if (file_names.empty()) {
        mtx_job_que.unlock();
        throw std::runtime_error("Finished Runs\n");
    }
    run = file_names.top();
    file_names.pop();
    mtx_job_que.unlock();
    return run;
}

void Analysis::UpdateData(Data_partial &partial_data) {
    data.TW_vs_ice.reserve(data.TW_vs_ice.size() + data_size);
    int ii{0};
    mtx_data.lock();
    while (partial_data.TW_vs_ice[ii][0] != 0 && partial_data.TW_vs_ice[ii][1] != 0) {
    //while (partial_data.TW_vs_ice[ii][0] != 0 ) {
        if (ii > data_size)
            throw std::runtime_error("Something wrong converting data in analysis class");

        data.TW_vs_ice.emplace_back(partial_data.TW_vs_ice[ii][0],
                                    partial_data.TW_vs_ice[ii][1]);
        std::cout << ii << std::endl;
        ++ii;
    }
    mtx_data.unlock();
}