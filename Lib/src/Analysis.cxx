#include "Analysis.h"

//ClassImp(Analysis);

Analysis::Analysis(int n_threads):
                    n_threads(n_threads)
{
    std::string file_with_runs = "./Configs/Runs.txt";
    std::ifstream file(file_with_runs);
    if (!file.is_open()) throw std::runtime_error(std::string("Unable to read :")+std::string(file_with_runs));
    std::string line;
    while(std::getline(file, line)){
        std::ifstream file_check(line);
        if (line.at(0)=='#') continue;
        if(file_check.fail()) throw std::runtime_error(std::string("File not present :")+ line);
        file_names.push(line);
    }
}

Analysis::~Analysis(){
}


bool Analysis::RunAnalysis(){
    for (int ii=0; ii<n_threads; ++ii){
        //mtx.lock();//added because ROOT is not thread safe
        threads.push_back(std::thread(&Analysis::Job, this));
        //mtx.unlock();//added because root is not thread safe
    }

    for (int ii=0; ii<n_threads; ++ii){
        threads.at(ii).join();
    }

    remove("./Out/sum.root");
    system("hadd -f ./Out/sum.root ./Out/*");

    return true;
}

bool Analysis::Job(){
    try{
        while(1){
            std::string current_run = GetRun();
            std::cout << "Run : " << current_run <<" assigned to thread\n";
            pid_t pid;
            int status;
            if ((pid=fork()) == 0){ //Child process
                RunSelector(current_run);
                _exit(1);
            }
            while((pid=wait(&status))>0) {	/* waitpid(-1, &status, 0) */
	    	    if(WIFEXITED(status)) {
	    	    	printf("%d returns %d\n",
	    	  	       pid, WEXITSTATUS(status));
	    	} else if(WIFSIGNALED(status)) {
	    		printf("%d terminated by %d\n",
	    		       pid, WTERMSIG(status));
	    	}
	}
            
        }
    }catch(std::runtime_error & e){
        std::cout << e.what();
    }
    return false;
}

bool Analysis::RunSelector(std::string run){
    TFile *file = new TFile(run.c_str());
    if (!file) throw std::runtime_error("Run : "+run+" Not opened\n"); 
    TTree *tree = (TTree *)file->Get("PhysicsTree");
    std::cout << "Starting selector for run : " << run << "\n";
    tree->Process(&selector, ("analyzed_"+run.substr(run.find_last_of("/")+1)).c_str());
    return true;
}

std::string Analysis::GetRun(){
    std::string run;
    mtx.lock();
    if (file_names.empty()) throw std::runtime_error("Finished Runs\n");
    run = file_names.top();
    file_names.pop();
    mtx.unlock();
    return run;
}
