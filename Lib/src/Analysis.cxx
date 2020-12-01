#include "Analysis.h"

//ClassImp(Analysis);

Analysis::Analysis(int n_threads) : n_threads(n_threads)
{
    std::string file_with_runs = "./Configs/Runs.txt";
    std::ifstream file(file_with_runs);
    if (!file.is_open())
        throw std::runtime_error(std::string("Unable to read :") + std::string(file_with_runs));
    std::string line;

    while (std::getline(file, line))
    {
        std::ifstream file_check(line);
        if (line.at(0) == '#')
            continue;
        if (file_check.fail())
            throw std::runtime_error(std::string("File not present :") + line);
        file_names.push(line);
    }

    std::ifstream file_test("./Configs/Interpolations/TW_Ice_Thickness.root");
    generate_TW_ice_interpolation = !file_test;
    
    threads_pid.resize(n_threads, 0);
}

Analysis::~Analysis()
= default;

bool Analysis::RunAnalysis()
{
    //Creating needed files
    system("cd Configs/Cuts/; root -l -b -q E_TOF_cuts.cxx; cd -;");
    system("cd Configs/Cuts/; root -l -b -q DE_E_cuts.cxx; cd -;");
    system("cd Configs/Cuts/; root -l -b -q MQ_Q_cuts.cxx; cd -;");
    system("cd Configs/Interpolations/; root -l -b -q TW_Brho_M46_Z18.C; cd -;");

    //Removing previous data
    system("rm -rf ./Out");
    system("mkdir ./Out");

    if (! NecessaryFilesPresent()) {
        n_threads = 1;
        std::cout << "n_threads moved to 1, as GraphsEnabled.txt file is missing\n";    
    }
    //Assigning to threads
    std::vector<int> exit_status;
    for (int ii = 0; ii < n_threads; ++ii)
    {
        exit_status.push_back(0);
        threads.emplace_back(&Analysis::Job, this, ii, &exit_status.back());
    }

    //Waiting for threads to finish
    for (int ii = 0; ii < n_threads; ++ii)
    {
        threads.at(ii).join();
    }

    for (const auto & stat: exit_status){
        if (stat != 0)
            throw stat;
    }

    std::cout << "Joined threads, adding output files :\n";

    //Adding all data in unique file
    system("hadd -f ./Out/sum.root ./Out/*");

    if (generate_TW_ice_interpolation)
    {
        std::sort(data.TW_vs_ice.begin(), data.TW_vs_ice.end());
        auto *root_file = new TFile("./Configs/Interpolations/TW_Ice_Thickness.root",
                                     "recreate");

        double X[data.TW_vs_ice.size()];
        double Y[data.TW_vs_ice.size()];
        for (long unsigned int ii = 0; ii < data.TW_vs_ice.size(); ++ii)
        {
            X[ii] = data.TW_vs_ice[ii].first;
            Y[ii] = data.TW_vs_ice[ii].second;
        }
        auto *gr = new TGraph(data.TW_vs_ice.size(), X, Y);
        auto *spl = new TSpline3("TW_vs_ice_thickness", gr);
        spl->Write();
        root_file->Write();
        root_file->Close();
    }

    return true;
}

//This job will be assigned to a thread
bool Analysis::Job(const int ind, int* exit_status)
{
    try
    {
        while (1)
        { //Repeat while there are jobs
            std::string current_run = GetRun();
            std::cout << "Run : " << current_run << " assigned to thread\n";

            pid_t pid;
            int status;

            int p[2];

            if (pipe(p) < 0)
                throw std::runtime_error("Pipe error");

            //Forking process to avoid ROOT threading problems
            Data_partial *partial_data = nullptr;

            mtx_fork.lock();
            bool fork_enabled = true;
            if (n_threads ==1)
                fork_enabled = false;

            if (fork_enabled)
                pid = fork();
            else //no need to fork!
                pid = 0;

            switch (pid){
            case -1:
                throw std::runtime_error("Fork error");
                break;
            case 0: //Child process
                close(p[0]);
                if (generate_TW_ice_interpolation){
                    partial_data = RunSelector(current_run);
                    FILE *out = fdopen(p[1], "w");
                    fwrite(partial_data, sizeof(*partial_data), 1, out);
                    fclose(out);
                }else{
                    RunSelector(current_run);
                }
                if (fork_enabled)
                    exit(0);
                //break;
            default:
                threads_pid[ind] = pid;
                mtx_fork.unlock();
                close(p[1]);

                if (generate_TW_ice_interpolation){
                    partial_data = new Data_partial();
                    FILE *in = fdopen(p[0], "r");
                    fread(partial_data, sizeof(*partial_data), 1, in);
                    fclose(in);
                }
                if (fork_enabled)
                    wait(&status); //Wait child process to complete instructions
                threads_pid[ind] = 0;
                if (WEXITSTATUS(status) == 1){
                    *exit_status=1;
                    return false;
                }
                    
            }
            if (generate_TW_ice_interpolation){
                UpdateData(*partial_data);
                if (partial_data)
                    delete partial_data;
            }
        }
    }
    catch (std::runtime_error &e){
        std::cout << e.what();
    }
    catch (int &e){
        std::cout << "Exception : " << e << std::endl;
    }
    return false;
}

Analysis::Data_partial *Analysis::RunSelector(const std::string& run) const
{
    auto *file = new TFile(run.c_str());
    if (!file->IsOpen())
        throw std::runtime_error("Run : " + run + " Not opened\n");
    auto *tree = (TTree *)file->Get("PhysicsTree");
    if (tree == nullptr)
        throw std::runtime_error("Tree not found in file\n");
    std::cout << "---------->Starting selector for run : " << run << "<----------\n";
    auto *selector = new Selector();

    tree->Process(selector, ("analyzed_" + run.substr(run.find_last_of('/') + 1)).c_str());
    if (generate_TW_ice_interpolation)
    {
        auto *partial_data = new Data_partial(selector->GetTWvsIce());
        delete selector;
        return partial_data;
    }
    else
    {
        delete selector;
        Data_partial *dummy_data = nullptr;
        return dummy_data;
    }
}

std::string Analysis::GetRun()
{
    std::string run;
    mtx_job_que.lock();
    if (file_names.empty())
    {
        mtx_job_que.unlock();
        throw std::runtime_error("Finished Runs\n");
    }
    run = file_names.top();
    file_names.pop();
    mtx_job_que.unlock();
    return run;
}

void Analysis::UpdateData(Data_partial &partial_data)
{
    data.TW_vs_ice.reserve(data.TW_vs_ice.size() + data_size);
    int ii{0};
    mtx_data.lock();
    while (partial_data.TW_vs_ice[ii][0] != 0 && partial_data.TW_vs_ice[ii][1] != 0)
    {
        //while (partial_data.TW_vs_ice[ii][0] != 0 ) {
        if (ii > data_size)
            throw std::runtime_error("Something wrong converting data in analysis class");

        data.TW_vs_ice.emplace_back(partial_data.TW_vs_ice[ii][0],
                                    partial_data.TW_vs_ice[ii][1]);
        ++ii;
    }
    mtx_data.unlock();
}

bool Analysis::NecessaryFilesPresent(){
    std::ifstream file_tmp("Configs/GraphsEnabled.txt");
    if (file_tmp.good()) return true;
    else return false;
}