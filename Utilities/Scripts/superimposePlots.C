#include <NPReaction.h>

void superimposePlots(){
    TFile* f = new TFile("outFinal.root", "read");
    TFile* f_out = new TFile("out3.root", "recreate");
    f_out->cd();
    std::set<int> mg = {1, 3, 4, 5, 7, 11};
    std::map<double, int> exs = {{0., kRed}, {2.02, kBlack}};
    std::map<int, int> colors = {{1,1}, {3,2}, {4,4}, {5,9}, {7,8}, {11, 35}};
    std::set<std::string> sp = {"h_corrected", "h_uncentered"};


    for (const auto& s: sp){
        for (int i=0;i<32;++i){
            std::cout << Form("%s_%i_%i", s.c_str(), i, 1) << std::endl;
            auto* n =  Form("%s_%i", s.c_str(), i);
            auto* hs = new THStack(n, n);
            std::vector<TF1*> fits;
            
            for (const auto& it: mg){
                auto*h = (TH1D*)f->Get(Form("%s_%i_%i", s.c_str(), i, it));
                TF1* gr = nullptr;
                if(it != 11){
                    gr = new TF1(Form("fit_%s_%i_%i", s.c_str(), i, it),
                                 "gaus(0)+gaus(3)",
                                 -2, 5);
                    //means
                    gr->FixParameter(1, 0);
                    gr->FixParameter(4, 2.02);

                    //ampliudes
                    gr->SetParameter(0, 10);
                    gr->SetParLimits(0, 0, 50);

                    gr->SetParameter(3, 10);
                    gr->SetParLimits(3, 0, 50);

                    //sigmas
                    gr->SetParameter(2, 1.5);
                    gr->SetParLimits(2, 0.3, 3);

                    gr->SetParameter(5, 1.5);
                    gr->SetParLimits(5, 0.3, 3);
                }else{
                    gr = new TF1(Form("fit_%s_%i_%i", s.c_str(), i, it),
                                 "gaus(0)",
                                 -2, 2);
                    //means
                    gr->SetParameter(1, 0);
                    gr->SetParLimits(1, 0, 0.4);

                    //ampliudes
                    gr->SetParameter(0, 25);
                    gr->SetParLimits(0, 0, 50);

                    //sigmas
                    gr->SetParameter(2, 1.5);
                    gr->SetParLimits(2, 0.3, 3);
                }
                h->SetLineColor(colors[it]);
                h->SetMarkerColor(colors[it]);
                h->SetDrawOption("histo");
                h->Rebin(2);
                h->SetLineWidth(3);
                //h->Scale(1./h->Integral());
                hs->Add(h);

                h->Fit(gr, "", ""); 
                gr->SetLineColor(colors[it]);
                fits.push_back(gr);
            }
            auto* c = new TCanvas(Form("%s_%i", s.c_str(), i));
            c->cd();
            //hs->Draw("histo nostack"); 
            hs->Draw("histo"); 
            for (const auto& it2: fits){
                it2->Draw("same");
            }
            //hs->Draw("histo"); 
            c->Write();
            c->Draw();
            //c->WaitPrimitive();
        }
    }
    
    f = new TFile("outFinal2.root", "read");
    
    for (const auto& s: sp){
        for (int i=0;i<32;++i){
            for (const auto& it: mg){
                auto* c = new TCanvas(Form("%s_%i_%i", s.c_str(), i, it));
                auto*h = (TH1D*)f->Get(Form("%s_%i_%i", s.c_str(), i, it));
                h->SetMarkerStyle(8);
                h->SetMarkerSize(0.5);
                h->Draw();
            
                NPL::Reaction reac("46Ar(3He,d)47K@458MeV"); 
                for (const auto exx: exs){
                    reac.SetExcitation4(exx.first);
                    auto* gg = reac.GetKinematicLine3();
                    gg->SetLineColor(exx.second);
                    gg->SetMarkerSize(3);
                    for (int i=0; i<gg->GetN(); ++i){
                        double x,y; 
                        gg->GetPoint(i,x,y);
                        gg->SetPoint(i,x*TMath::DegToRad(), y);
                    }
                    gg->Draw("same L");
                }
                f_out->cd();
                c->Write();

            }
        }
    }
    f_out->Write();

}
