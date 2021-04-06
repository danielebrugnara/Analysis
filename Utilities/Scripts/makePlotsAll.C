#include <sys/wait.h>


std::vector<int> jobs;
std::vector<std::thread> threads;
std::mutex mtx;
std::mutex mtx_fork;

void run(int);
void run2(int);
void run3(int);
void getJob();


void makePlotsAll(){
    int n_threads = 10;
    for (int i=15;i<32;++i){
        jobs.emplace_back(i);
    }
    for (int i=0;i<n_threads;++i){
        threads.emplace_back(&getJob);
    }

    for (int i=0;i<n_threads;++i){
        std::cout << "try to join " << i << std::endl;
        threads.at(i).join();
        std::cout << "joined " << i << std::endl;
    }
    system("hadd outFinished.root out_*.root");
}


void getJob(){
    while (1){
        mtx.lock();
        if (jobs.empty()) break;
        int val = jobs.back();
        jobs.pop_back();
        mtx.unlock();

        pid_t pid;
        int status;

        mtx_fork.lock();

        pid = fork();
        switch(pid){
            case -1:
                throw std::runtime_error("Fork error");
                break;
            case 0:
                //run(val);
                //run2(val);
                run3(val);
                exit(0);
            default:
                mtx_fork.unlock();
                waitpid(pid, &status, 0);
        }
    }

}




void run(int number){

    TFile f("sum.root", "read");
    TTree *tree = (TTree*)f.Get("AnalyzedTree");
    TFile* f_out = new TFile(Form("out_%i.root", number), "recreate");
    std::set<int> mg = {1, 3, 4, 5, 7, 11};
    int i=number;
    for (const auto& it: mg){
        std::string histo_name = Form("h_corrected_%i_%i", i,it);
        auto* histo = new TH1D(histo_name.c_str(),histo_name.c_str(),150, -15, 15);
        std::cout << "h_corrected" << i << std::endl;
        tree->Draw(Form("MugastData.Ex_corrected[0][%i]>>%s", i, histo_name.c_str()), Form("MugastData.M ==2 && MugastData.Z == 1 && VamosData.id_Z == 19 && VamosData.id_M == 47 && MugastData.MG == %i",it));
        histo->Write();

    }
    for (const auto& it: mg){
        std::string histo_name = Form("h_uncentered_%i_%i", i,it);
        auto* histo = new TH1D(histo_name.c_str(), histo_name.c_str(),150, -15, 15);
        std::cout << "h_uncentered" << i << std::endl;
        tree->Draw(Form("MugastData.Ex_uncentered[0][%i]>>%s", i, histo_name.c_str()), Form("MugastData.M ==2 && MugastData.Z == 1 && VamosData.id_Z == 19 && VamosData.id_M == 47 && MugastData.MG == %i", it));
        histo->Write();

    }
    f_out->Write();
    f_out->Close();
}

void run2(int number){

    TFile f("sum.root", "read");
    TTree *tree = (TTree*)f.Get("AnalyzedTree");
    TFile* f_out = new TFile(Form("out2_%i.root", number), "recreate");
    std::set<int> mg = {1, 3, 4, 5, 7, 11};
    int i=number;
    for (const auto& it: mg){
        std::string histo_name = Form("h_corrected_%i_%i", i,it);
        auto* histo = new TH2D(histo_name.c_str(),histo_name.c_str(),1000, 0, TMath::Pi(), 1000, 0, 15);
        std::cout << "h_corrected" << i << std::endl;
        tree->Draw(Form("MugastData.E_corrected[0][%i]:MugastData.EmissionDirection_corrected[%i].Theta()>>%s", i, i, histo_name.c_str()), Form("MugastData.M ==2 && MugastData.Z == 1 && VamosData.id_Z == 19 && VamosData.id_M == 47 && MugastData.MG == %i",it));
        histo->Write();

    }
    for (const auto& it: mg){
        std::string histo_name = Form("h_uncentered_%i_%i", i,it);
        auto* histo = new TH2D(histo_name.c_str(),histo_name.c_str(),1000, 0, TMath::Pi(), 1000, 0, 15);
        std::cout << "h_uncentered" << i << std::endl;
        tree->Draw(Form("MugastData.E_uncentered[0][%i]:MugastData.EmissionDirection_uncentered[%i].Theta()>>%s", i, i, histo_name.c_str()), Form("MugastData.M ==2 && MugastData.Z == 1 && VamosData.id_Z == 19 && VamosData.id_M == 47 && MugastData.MG == %i", it));
        histo->Write();

    }
    f_out->Write();
    f_out->Close();
}

void run3(int number){

    TFile f("sum.root", "read");
    TTree *tree = (TTree*)f.Get("AnalyzedTree");
    if (tree == nullptr) std::cout << "error\n";
    TFile* f_out = new TFile(Form("out3_%i.root", number), "recreate");
    std::set<int> mg = {1, 3, 4, 5, 7, 11};
    int i=number;
    for (const auto& it: mg){
        std::string histo_name = Form("h_corrected_%i_%i", i,it);
        auto* histo = new TH2D(histo_name.c_str(),histo_name.c_str(),1000, 0, TMath::Pi(), 1000, -15, 15);
        std::cout << "h_corrected" << i << std::endl;
        tree->Draw(Form("MugastData.Ex_corrected[0][%i]:MugastData.EmissionDirection_corrected[%i].Theta()>>%s", i, i, histo_name.c_str()), Form("MugastData.M ==2 && MugastData.Z == 1 && VamosData.id_Z == 19 && VamosData.id_M == 47 && MugastData.MG == %i",it));
        histo->Write();

    }
    for (const auto& it: mg){
        std::string histo_name = Form("h_uncentered_%i_%i", i,it);
        auto* histo = new TH2D(histo_name.c_str(),histo_name.c_str(),1000, 0, TMath::Pi(), 1000, -15, 15);
        std::cout << "h_uncentered" << i << std::endl;
        tree->Draw(Form("MugastData.Ex_uncentered[0][%i]:MugastData.EmissionDirection_uncentered[%i].Theta()>>%s", i, i, histo_name.c_str()), Form("MugastData.M ==2 && MugastData.Z == 1 && VamosData.id_Z == 19 && VamosData.id_M == 47 && MugastData.MG == %i", it));
        histo->Write();

    }
    f_out->Write();
    f_out->Close();
}
