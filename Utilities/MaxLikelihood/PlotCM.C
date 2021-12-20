void PlotCM(){
    auto* f = new TFile("histofile.root", "read");
    auto* data = (TH1D*)f->Get("dataCMNormalized");
    auto* simu = (TH1D*)f->Get("simuCMNormalized");
    auto* l0 = (TGraph*)f->Get("gs12");
    auto* l2 = (TGraph*)f->Get("gd32");
    auto* l3 = (TGraph*)f->Get("gf72");
    auto* ls = (TGraph*)f->Get("gsum");

    auto* cv = new TCanvas();
    cv->SetLogy();
    double scalefac=2.5;
    data->Scale(scalefac);
    simu->Scale(scalefac);

    data->SetMarkerStyle(8);
    simu->SetMarkerStyle(8);
    data->SetMarkerSize(0.7);
    simu->SetMarkerSize(0.7);
    data->SetMarkerColorAlpha(kBlack, 0.8);
    data->SetLineColorAlpha(kBlack, 0.8);
    simu->SetMarkerColorAlpha(kRed+1, 0.8);
    simu->SetLineColorAlpha(kRed+1, 0.);
    data->SetTitle(";Center of Mass Angle [deg]; d#sigma/d#Omega [mb/sr]");
    simu->SetTitle(";Center of Mass Angle [deg]; d#sigma/d#Omega [mb/sr]");
    data->GetXaxis()->SetLabelSize(0.05);
    data->GetXaxis()->SetTitleSize(0.05);
    data->GetYaxis()->SetLabelSize(0.05);
    data->GetYaxis()->SetTitleSize(0.05);
    data->GetYaxis()->SetTitleOffset(0.9);
    
    
    l0->SetLineColorAlpha(kBlue, 0.8);
    l2->SetLineColorAlpha(kOrange+1, 0.8);
    l3->SetLineColorAlpha(kGreen+2, 0.8);
    l0->SetLineStyle(8);
    l2->SetLineStyle(8);
    l3->SetLineStyle(8);
    ls->SetLineColorAlpha(kRed+1, 0.7);
    l0->SetLineWidth(3);
    l2->SetLineWidth(3);
    l3->SetLineWidth(3);
    ls->SetLineWidth(6);

    data->Draw("e1");
    data->GetXaxis()->SetRangeUser(0, 23);
    data->GetYaxis()->SetRangeUser(0.5, 100);
    simu->Draw("same, p");
    l0->Draw("same, l");
    l2->Draw("same, l");
    l3->Draw("same, l");
    ls->Draw("same, l");
    auto* leg = new TLegend(0.5, 0.5, 0.9, 0.9);
    leg->AddEntry(data, "Data", "ep");
    leg->AddEntry(simu, "Simulation", "ep");
    leg->AddEntry(l0, "L=0", "l");
    leg->AddEntry(l2, "L=2", "l");
    leg->AddEntry(l3, "L=3", "l");
    leg->AddEntry(ls, "Convolution", "l");
    leg->Draw();
}
