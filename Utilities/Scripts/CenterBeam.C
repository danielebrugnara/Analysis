#include <vector>
#include <TH1D.h>
#include <TH2D.h>

void CenterBeam(TH2* h1, TH2* h2, TH2* h3, TH2*h4){

//    TFile f("center_beam.root");
    std::vector<std::pair<TH2*, TH1D*>> histos;
    TH2* tmps[4] = {h1, h2, h3, h4};

    for (int i=0; i< 4; ++i){
//        auto* tmp = (TH2D*)f.Get(Form("h%i", i+1));
        auto* tmp = tmps[i];  
        auto* tmp2 = (TH1D*)tmp->ProjectionY((std::string(tmp->GetName())+"_normalized").c_str());
        histos.emplace_back(tmp, tmp2);

    }

    for (const auto& histo: histos){

        for ( int i=0; i<histo.first->GetNbinsY(); ++i){
            int nbins = 0;
            for ( int j=0; j<histo.first->GetNbinsX(); ++j){
                int content = histo.first->GetBinContent(j,i);
                if (content > 0) ++nbins;
            }
            if (histo.second->GetBinContent(i)>0){
                histo.second->SetBinContent(i, histo.second->GetBinContent(i)/nbins);
            }
        }
        new TCanvas;
        histo.second->Draw();
        histo.second->Rebin(8);
    }
    new TCanvas;
    histos[0].second->Draw();
    histos[1].second->SetLineColor(kRed);
    histos[1].second->Draw("same");
    histos[2].second->SetLineColor(kBlue);
    histos[2].second->Draw("same");
    histos[3].second->SetLineColor(kGreen);
    histos[3].second->Draw("same");
}
