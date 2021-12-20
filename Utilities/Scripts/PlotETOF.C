#include "TTree.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TH2D.h"
#include "TLegend.h"

void PlotETOF(){
    std::vector<int> MG{1,3,4,5,7,11};
    std::map<int, double> minLim;
    std::map<int, double> maxLim;

    minLim[1] = -2.5;
    maxLim[1] = -9;
    minLim[3] = -2.5;
    maxLim[3] = -13;
    minLim[4] = -4;
    maxLim[4] = -14;
    minLim[5] = -4;
    maxLim[5] = -13;
    minLim[7] = -5;
    maxLim[7] = -20;
    minLim[11] = -2;
    maxLim[11] = -7;
    

    std::vector<std::string> gates{"E_TOF_m1_z1_MG", "E_TOF_m2_z1_MG", "E_TOF_m4_z2_MG"};
    TFile* f = new TFile("../../build/Out/sum.root", "read");
    TFile* fcuts = new TFile("../../build/Configs/Cuts/MUGAST.root", "read");
    TTree* t = (TTree*) f->Get("AnalyzedTree");
    TTree* td = (TTree*) f->Get("AnalyzedTreeDeuterons");
    for(const auto&it: MG){
        //bool useProtons{false};
        bool useProtons{true};
        std::string name = std::to_string(it);
        TCanvas* cv = new TCanvas(name.c_str(), name.c_str());
        TH2D* h   = nullptr;
        TH2D* hk  = nullptr;
        TH2D* hkd = nullptr;
        if (useProtons){
            h = new TH2D(name.c_str(), name.c_str(), 1000, 0, 30, 1000, -30, 30);
            hk = new TH2D((name + "k").c_str(), (name + "k").c_str(), 1000, 0, 30, 1000, -30, 30);
            hkd = new TH2D((name + "kd").c_str(), (name + "kd").c_str(), 1000, 0, 30, 1000, -30, 30);
            t ->Draw(("MugastData.T_proton:MugastData.SI_E>>"+name).c_str(), Form("MugastData.MG == %i", it), "");
            t ->Draw(("MugastData.T_proton:MugastData.SI_E>>"+(name+"k")).c_str(), Form("MugastData.MG == %i && VamosData.id_Z == 19 && VamosData.id_M == 47", it), "same");
            td->Draw(("MugastData.T_proton:MugastData.SI_E>>"+(name+"kd")).c_str(), Form("MugastData.MG == %i", it), "same");
            //t->Draw(("MugastData.T_proton:MugastData.SI_E>>"+(name+"kd")).c_str(), Form("VamosData.id_Z == 19 && VamosData.id_M == 47 && MugastData.T_proton<%f && MugastData.T_proton>%f && MugastData.MG == %i",minLim[it], maxLim[it], it), "same");


            //Plot gates
            for(const auto& gate:gates){
                std::string gateName = gate+std::to_string(it);    
                std::cout << "drawing gate " << gateName << std::endl;
                ((TCutG*)fcuts->Get(gateName.c_str()))->Draw("same");
            }
        }else {
            h = new TH2D(name.c_str(), name.c_str(), 1000, 0, 30, 1000, 250, 400);
            hk = new TH2D((name + "k").c_str(), (name + "k").c_str(), 1000, 0, 30, 1000, 250, 400);
            hkd = new TH2D((name + "kd").c_str(), (name + "kd").c_str(), 1000, 0, 30, 1000, 250, 400);
            t ->Draw(("MugastData.T:MugastData.SI_E>>"+name).c_str(), Form("MugastData.MG == %i", it), "");
            t ->Draw(("MugastData.T:MugastData.SI_E>>"+(name+"k")).c_str(), Form("MugastData.MG == %i && VamosData.id_Z == 19 && VamosData.id_M == 47", it), "same");
            td->Draw(("MugastData.T:MugastData.SI_E>>"+(name+"kd")).c_str(), Form("MugastData.MG == %i", it), "same");
            //t->Draw(("MugastData.T:MugastData.SI_E>>"+(name+"kd")).c_str(), Form("VamosData.id_Z == 19 && VamosData.id_M == 47 && MugastData.T_proton<%f && MugastData.T_proton>%f && MugastData.MG == %i",minLim[it], maxLim[it], it), "same");
        }

        h->SetTitle(";Energy [MeV]; TOF [ns]");
        h->GetXaxis()->SetTitleOffset(0.9);
        h->GetYaxis()->SetTitleOffset(0.9);
        h->GetXaxis()->SetTitleSize(0.05);
        h->GetYaxis()->SetTitleSize(0.05);
        h->GetXaxis()->SetLabelSize(0.05);
        h->GetYaxis()->SetLabelSize(0.05);

        auto* legend = new TLegend(0.4, 0.1, 0.9, 0.4);
        legend->AddEntry(h, "No selection", "pf");
        legend->AddEntry(hk, "^{47}K in VAMOS", "pf");
        legend->AddEntry(hkd, "^{47}K in vamos and ^{2}H in MUGAST", "pf");
        legend->Draw();

        h->SetLineColor(kBlack);
        hk->SetMarkerStyle(4);
        hk->SetMarkerSize(0.7);
        hk->SetMarkerColor(kBlue+1);
        hk->SetLineColor(kBlue+1);
        hkd->SetMarkerStyle(8);
        hkd->SetMarkerSize(0.7);
        hkd->SetMarkerColor(kRed+1);
        hkd->SetLineColor(kRed+1);
        cv->Update();
    }
}
