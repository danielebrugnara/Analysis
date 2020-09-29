void CheckGeom(){
    TFile* file_dati = new TFile("./spectra_Tree_run99.root", "read");

    TFile* file_simu_2 = new TFile("./spectra_GammaEvents.0065.root", "read");
    TFile* file_simu = new TFile("./spectra_GammaEvents.0066.root", "read");

    auto* geom_dati =(TH3D*)file_dati->Get("geom");
    geom_dati->Scale(1./geom_dati->Integral("width"));
    auto* geom_simu =(TH3D*)file_simu->Get("geom");
    geom_simu->Scale(1./geom_simu->Integral("width"));
    auto* geom_simu_2 =(TH3D*)file_simu_2->Get("geom");
    geom_simu_2->Scale(1./geom_simu_2->Integral("width"));

    TH1D* projx_dati = geom_dati->ProjectionX("px_dati");
    TH1D* projx_simu = geom_simu->ProjectionX("px_simu");
    TH1D* projx_simu_2 = geom_simu_2->ProjectionX("px_simu_2");
    
    TH1D* projy_dati = geom_dati->ProjectionY("py_dati");
    TH1D* projy_simu = geom_simu->ProjectionY("py_simu");
    TH1D* projy_simu_2 = geom_simu_2->ProjectionY("py_simu_2");

    TH1D* projz_dati = geom_dati->ProjectionZ("pz_dati");
    TH1D* projz_simu = geom_simu->ProjectionZ("pz_simu");
    TH1D* projz_simu_2 = geom_simu_2->ProjectionZ("pz_simu_2");

    TH2D* projxy_dati = (TH2D*)geom_dati->Project3D("yx");
    projxy_dati->SetName("xy_dati");
    TH2D* projxy_simu = (TH2D*)geom_simu->Project3D("yx");
    projxy_simu->SetName("xy_simu");
    TH2D* projxy_simu_2 = (TH2D*)geom_simu_2->Project3D("yx");
    projxy_simu_2->SetName("xy_simu_2");

    TH2D* projxy_diff = (TH2D*)projxy_dati->Clone("projxy_diff");
    projxy_diff->Add(projxy_simu, -1);

    auto* cvx = new TCanvas("projx");
    projx_dati->SetMarkerColor(kRed);
    projx_dati->SetLineColor(kRed);
    projx_dati->Draw("histo");
    projx_simu->Draw("same histo");
    projx_simu_2->SetLineColor(kGreen);
    projx_simu_2->Draw("same histo");
    
    auto* cvy = new TCanvas("projy");
    projy_dati->SetMarkerColor(kRed);
    projy_dati->SetLineColor(kRed);
    projy_dati->Draw("histo");
    projy_simu->Draw("same histo");
    projy_simu_2->SetLineColor(kGreen);
    projy_simu_2->Draw("same histo");

    auto* cvz = new TCanvas("projz");
    projz_dati->SetMarkerColor(kRed);
    projz_dati->SetLineColor(kRed);
    projz_dati->Draw("histo");
    projz_simu->Draw("same histo");
    projz_simu_2->SetLineColor(kGreen);
    projz_simu_2->Draw("same histo");

    auto* cvxy_dati = new TCanvas("projxy_dati");
    projxy_dati->Draw("colz");

    auto* cvxy_simu = new TCanvas("projxy_simu");
    projxy_simu->Draw("colz");

    auto* cvxy_simu_2 = new TCanvas("projxy_simu_2");
    projxy_simu_2->Draw("colz");

    auto* cvxy_diff = new TCanvas("projxy_diff");
    projxy_diff->Draw("colz");
}
