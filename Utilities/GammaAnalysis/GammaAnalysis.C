#include <TFile.h>
#include <TH1D.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TLegend.h>

#include <fstream>


void GammaAnalysis(){

    std::ifstream datafile("data.root");
    TH1D* h{nullptr};
    TH1D* hbkg{nullptr};
    if(!datafile.is_open()) {

        auto *file = TFile::Open("../../DataAnalyzed/sum.root", "read");
        auto *tree = static_cast<TTree *>(file->Get("AnalyzedTree"));

        int nbins{100};
        double min{50};
        double max{2000};

        h = new TH1D("h", ";Energy [keV]; Counts", nbins, min, max);
        tree->Draw("AgataData.EDC>>h", "AgataData.in_coincidence && VamosData.id_M == 47 && VamosData.id_Z == 19", "");

        hbkg = new TH1D("hbkg", "hbkg", nbins, min, max);
        tree->Draw("AgataData.EDC>>hbkg", "!AgataData.in_coincidence", "");

        auto* datafile = TFile::Open("data.root", "recreate");
        h->Write();
        hbkg->Write();
        datafile->Write();
        datafile->Close();

    }else{
        auto* datafile = TFile::Open("data.root", "read");
        h = static_cast<TH1D *>(datafile->Get("h"));
        hbkg = static_cast<TH1D *>(datafile->Get("hbkg"));
    }


    h->Add(hbkg, -0.03);

    h->SetMarkerColor(kBlack);
    h->SetMarkerStyle(8);
    h->SetLineColor(kBlack);



    auto* fileSimu = TFile::Open("../../DataAnalyzed/gammaspec_simulated_2020.root", "read");
    auto* simu = static_cast<TH1D*>(fileSimu->Get("addb_spec_DC"));
    simu->Scale(0.0011);
    simu->Rebin(16);
    simu->SetLineColor(kRed+1);

    gStyle->SetOptStat(0);
    auto* cv = new TCanvas();
    //cv->SetLogy();
    h->GetXaxis()->SetLabelSize(0.05);
    h->GetXaxis()->SetTitleSize(0.05);
    h->GetXaxis()->SetTitleOffset(0.95);
    h->GetXaxis()->SetRangeUser(50, 2000);

    h->GetYaxis()->SetLabelSize(0.05);
    h->GetYaxis()->SetTitleSize(0.05);
    h->GetYaxis()->SetTitleOffset(0.7);
    h->GetYaxis()->SetRangeUser(-5, 80);

    h->Draw("E1");
    simu->Draw("histo, same");

    auto* legend = new TLegend(0.7, 0.7, 0.9, 0.9);
    legend->AddEntry(h, "Simulation");
    legend->AddEntry(simu, "Data");

    legend->Draw();

    cv->SaveAs("gamma.pdf");
}
