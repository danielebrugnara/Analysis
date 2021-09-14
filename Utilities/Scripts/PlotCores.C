#include "tktToRoot.C"

#include <TH1D.h>
#include <TCanvas.h>
#include <THStack.h>


#include <unistd.h>
#include <dirent.h>
#include <string>
#include <map>

std::vector<std::string> getFolders(const std::string& pos);


void PlotCores(){
    std::string folder{"./GammaSpectra/"};
    std::vector<std::string> folders = getFolders(folder);
    std::map<std::string, TH1D*> histogramsBefore;
    std::map<std::string, TH1D*> histogramsAfter;

    int color{1};
    auto* stackBefore = new THStack();
    auto* stackAfter = new THStack();
    for(const auto&it: folders){
        try {
            histogramsBefore[it] = tktToRoot(folder + it + "/Post__5-40-16384-UI__Ener.spec",
                                       std::vector<unsigned int>{4, 39},
                                       1);
            histogramsBefore[it]->SetName(it.c_str());
            histogramsAfter[it] = tktToRoot(folder + it + "/Post__5-40-16384-UI__Ener.spec",
                                             std::vector<unsigned int>{4, 39},
                                             1);
            histogramsAfter[it]->SetName(it.c_str());

            if (it != "06C") {
                stackBefore->Add(histogramsBefore[it]);
                stackAfter->Add(histogramsAfter[it]);
                color +=2;
                histogramsBefore[it]->SetLineColorAlpha(color, 0.2);
                histogramsAfter[it]->SetLineColorAlpha(color, 0.1);
                histogramsBefore[it]->SetName(it.c_str());
                histogramsAfter[it]->SetName(it.c_str());
                histogramsBefore[it]->SetTitle(it.c_str());
                histogramsAfter[it]->SetTitle(it.c_str());
            }
        } catch (const std::exception& e) {
            continue;
        }
        //histogramsBefore[it]->Draw("same");
    }

    auto* cv = new TCanvas();
    cv->Divide(2,1);
    cv->cd(1)->SetLogy();
    stackBefore->Draw("nostack");
    stackBefore->GetXaxis()->SetTitle("Energy [keV]");
    stackBefore->GetYaxis()->SetTitle("Counts");
    stackBefore->GetXaxis()->SetLabelSize(0.05);
    stackBefore->GetXaxis()->SetTitleSize(0.05);
    stackBefore->GetYaxis()->SetLabelSize(0.05);
    stackBefore->GetYaxis()->SetTitleSize(0.05);
    stackBefore->GetYaxis()->SetTitleOffset(0.9);
    cv->cd(2);
    stackAfter->Draw("nostack");
    stackAfter->GetXaxis()->SetRangeUser(100, 400);
    stackAfter->GetXaxis()->SetTitle("Energy [keV]");
    stackAfter->GetYaxis()->SetTitle("Counts");
    stackAfter->GetXaxis()->SetLabelSize(0.05);
    stackAfter->GetXaxis()->SetTitleSize(0.05);
    stackAfter->GetYaxis()->SetLabelSize(0.05);
    stackAfter->GetYaxis()->SetTitleSize(0.05);
    stackAfter->GetYaxis()->SetTitleOffset(0.9);
}

std::vector<std::string> getFolders(const std::string& pos) {
    DIR *dir;
    struct dirent *ent;
    std::vector<std::string> files;
    if ((dir = opendir(pos.c_str())) != NULL) {
    /* print all the files and directories within directory */
        while ((ent = readdir(dir)) != NULL) {
            std::string tmpname(ent->d_name);
            if (tmpname.size()!=3) continue;

            printf("%s\n", ent->d_name);
            files.emplace_back(tmpname);
        }
        closedir(dir);
    } else {
/* could not open directory */
        perror("");
        throw std::runtime_error("error\n");
    }
    return files;
}

