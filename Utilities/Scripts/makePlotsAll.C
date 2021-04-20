#include <sys/wait.h>
#include <unistd.h>
#include <dirent.h>

#include <iostream>
#include <mutex>
#include <thread>
#include <vector>

std::vector<std::string> jobs;
std::mutex mtx;
std::mutex mtx_fork;

void job(const std::string&);
void jobAna(const std::string&);
void getJob();
void findJob();

void makePlotsAll() {
    int n_threads = 10;
    std::vector<std::thread> threads;
    //for (int i=15;i<32;++i){
    //    jobs.emplace_back(i);
    //}
    findJob();
    for (int i = 0; i < n_threads; ++i) {
        threads.emplace_back(&getJob);
    }

    for (int i = 0; i < n_threads; ++i) {
        std::cout << "try to join " << i << std::endl;
        threads.at(i).join();
        std::cout << "joined " << i << std::endl;
    }
    //system("hadd outFinished.root out_*.root");
}

int main(){
    makePlotsAll();
}

void getJob() {
    while (1) {
        mtx.lock();
        if (jobs.empty()) {
            mtx.unlock();
            break;
        }
        auto val = jobs.back();
        jobs.pop_back();
        mtx.unlock();

        pid_t pid;
        int status;

        mtx_fork.lock();

        pid = fork();
        switch (pid) {
            case -1:
                throw std::runtime_error("Fork error");
                break;
            case 0:
                job(val);
                jobAna(val);
                exit(0);
            default:
                mtx_fork.unlock();
                waitpid(pid, &status, 0);
        }
    }
}

void job(const std::string& val) {
    std::string command;
    command += "npsimulation -E ";  //./Reaction/46Ar3Hed47K_0keV_s12_franco.reaction
    command += "./Reaction/";  //./Reaction/46Ar3Hed47K_0keV_s12_franco.reaction
    command += val;
    command += " -D ./Detector/MUGAST_cryotarget_franco.detector -O ";
    command += val.substr(val.find_last_of("/") + 1, val.find_last_of(".") - val.find_last_of("/") - 1);
    command += ".root "; 
    command += " -B 46Ar_single.mac";

    std::cout << "Launching job : " << command << std::endl;
    //system(command.c_str());
}

void jobAna(const std::string& val) {
    std::string command;
    command += "npanalysis -T ./simout/";
    command += val.substr(val.find_last_of("/") + 1, val.find_last_of(".") - val.find_last_of("/") - 1);
    command += ".root SimulatedTree"; 
    command += " -D ./Detector/MUGAST_cryotarget_franco.detector -O ";
    command += val.substr(val.find_last_of("/") + 1, val.find_last_of(".") - val.find_last_of("/") - 1);
    command += ".root ";  

    std::cout << "Launching job : " << command << std::endl;
    system(command.c_str());

}

void findJob() {
    DIR* dir;
    struct dirent* ent;
    if ((dir = opendir("./Reaction/")) != NULL) {
        /* print all the files and directories within directory */
        while ((ent = readdir(dir)) != NULL) {
            std::string tmpname(ent->d_name);
            std::string necessaryBegin = "46Ar3Hed47K_";
            if(tmpname.compare(0, necessaryBegin.size(), necessaryBegin) != 0 ) continue;
            printf("%s\n", ent->d_name);
            jobs.emplace_back(tmpname);
        }
        closedir(dir);
    } else {
        /* could not open directory */
        perror("");
        return;
    }
}
