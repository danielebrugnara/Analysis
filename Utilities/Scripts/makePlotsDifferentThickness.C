#include <sys/wait.h>
#include <unistd.h>
#include <dirent.h>

#include <iostream>
#include <mutex>
#include <thread>
#include <vector>
#include <string>

struct Job{
  std::string detector;
  std::string reaction;
};

std::vector<Job> jobs;
std::mutex mtx;
std::mutex mtx_fork;

void jobSimu(const Job&);
void jobAna(const Job&);
void getJob();
void findJob();

void makePlotsAll() {
    int n_threads = 10;
    std::vector<std::thread> threads;
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
                jobSimu(val);
                jobAna(val);
                exit(0);
            default:
                mtx_fork.unlock();
                waitpid(pid, &status, 0);
        }
    }
}

void jobSimu(const Job& val) {
    std::string command;
    command += "npsimulation -E ";  //./Reaction/46Ar3Hed47K_0keV_s12_franco.reaction
    command += val.reaction;
    command += " -D ./Detector/";
    command += val.detector;                    //./Detector/MUGAST_cryotarget_franco.detector
    command += " -O ";
    command += val.reaction.substr(val.reaction.find_last_of("/") + 1, val.reaction.find_last_of(".") - val.reaction.find_last_of("/") - 1);
    command += "_";
    command += val.detector.substr(val.detector.find_last_of("_") + 1 , val.detector.find_last_of(".")- val.detector.find_last_of("_") - 1 );
    command += ".root "; 
    command += " -B 46Ar_single.mac";

    std::cout << "Launching job : " << command << std::endl;
    system(command.c_str());
}

void jobAna(const Job& val) {
    std::string command;
    command += "npanalysis -T ./simout/";
    command += val.reaction.substr(val.reaction.find_last_of("/") + 1, val.reaction.find_last_of(".") - val.reaction.find_last_of("/") - 1);
    command += "_";
    command += val.detector.substr(val.detector.find_last_of("_") + 1 , val.detector.find_last_of(".")- val.detector.find_last_of("_") - 1 );
    command += ".root SimulatedTree"; 
    command += " -D ./Detector/";
    command += val.detector; //./Detector/MUGAST_cryotarget_franco.detector -O
    command += " -O ";
    command += val.reaction.substr(val.reaction.find_last_of("/") + 1, val.reaction.find_last_of(".") - val.reaction.find_last_of("/") - 1);
    command += "_";
    command += val.detector.substr(val.detector.find_last_of("_") + 1 , val.detector.find_last_of(".")- val.detector.find_last_of("_") - 1 );
    command += ".root ";  

    std::cout << "Launching job : " << command << std::endl;
    //system(command.c_str());

}

void findJob() {

    std::vector<std::string> tmpDetector;
    DIR* dir;
    struct dirent* ent;
    if ((dir = opendir("./Detector/")) != NULL) {
        /* print all the files and directories within directory */
        while ((ent = readdir(dir)) != NULL) {
            std::string tmpname(ent->d_name);
            std::string necessaryBegin = "MUGAST_cryotarget_";
            if(tmpname.compare(0, necessaryBegin.size(), necessaryBegin) != 0 ) continue;
            printf("%s\n", ent->d_name);
            tmpDetector.emplace_back(tmpname);
        }
        closedir(dir);
    } else {
        /* could not open directory */
        perror("");
        return;
    }

    std::vector<std::string> tmpReaction{
      "./Reaction/46Ar3Hed47K_0keV_s12_0_0.reaction",
      "./Reaction/46Ar3Hed47K_360keV_d32_0_0.reaction",
      "./Reaction/46Ar3Hed47K_2020keV_f72_0_0.reaction",
      "./Reaction/46Ar3Hed47K_0keV_flat_0_0.reaction",
      "./Reaction/46Ar3Hed47K_360keV_flat_0_0.reaction",
      "./Reaction/46Ar3Hed47K_2020keV_flat_0_0.reaction"
    };

    for (const auto& itDetector: tmpDetector){
      for (const auto& itReaction: tmpReaction){
        Job tmpJob;
        tmpJob.reaction = itReaction;
        tmpJob.detector = itDetector;

        jobs.push_back(tmpJob);
      }
    }
}
