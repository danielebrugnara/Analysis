void PlotIcVsIc(){
    TFile* fd = new TFile("~/Projects/Analysis/build/Out/sum.root", "read");
    TFile* fs = new TFile("~/Projects/Simulation/Vamos/build/SimuOutput.root", "read");
    TTree* td = (TTree*)fd->Get("AnalyzedTree");
    TTree* ts = (TTree*)fs->Get("ICDet");

    //s1
    TH2D* hd1 = new TH2D("hd1", ";Row_{0} [MeV];Row_{1} [MeV]", 1000, 0, 60, 1000, 0, 60);
    TH2D* hs1 = new TH2D("hs1", ";Row_{0} [MeV];Row_{1} [MeV]", 1000, 0, 60, 1000, 0, 60);

    auto* cv1 = new TCanvas();
    cv1->SetLogz();
    hd1->GetXaxis()->SetTitleSize(0.04);
    hd1->GetYaxis()->SetTitleSize(0.04);
    hd1->GetXaxis()->SetLabelSize(0.04);
    hd1->GetYaxis()->SetLabelSize(0.04);
    //hd1->GetXaxis()->SetLabelOffset(0.9);
    //hd1->GetYaxis()->SetLabelOffset(0.9);
    hs1->SetMarkerStyle(8);
    hs1->SetMarkerSize(0.7);
    hs1->SetMarkerColorAlpha(kRed, 0.03);
    td->Draw("VamosData.IC[1]*1.02-2:VamosData.IC[0]>>hd1", "", "col");
    ts->Draw("IC1:IC0>>hs1", "", "same");

    //s2
    TH2D* hd2 = new TH2D("hd2", ";Row_{1} [MeV];Row_{2} [MeV]", 1000, 0, 60, 1000, 0, 120);
    TH2D* hs2 = new TH2D("hs2", ";Row_{1} [MeV];Row_{2} [MeV]", 1000, 0, 60, 1000, 0, 120);

    auto* cv2 = new TCanvas();
    cv2->SetLogz();
    hd2->GetXaxis()->SetTitleSize(0.04);
    hd2->GetYaxis()->SetTitleSize(0.04);
    hd2->GetXaxis()->SetLabelSize(0.04);
    hd2->GetYaxis()->SetLabelSize(0.04);
    //hd2->GetXaxis()->SetLabelOffset(0.9);
    //hd2->GetYaxis()->SetLabelOffset(0.9);
    hs2->SetMarkerStyle(8);
    hs2->SetMarkerSize(0.7);
    hs2->SetMarkerColorAlpha(kRed, 0.03);
    td->Draw("VamosData.IC[2]*1.072-6.4:VamosData.IC[1]*1.02-2>>hd2", "", "col");
    ts->Draw("IC2:IC1>>hs2", "", "same");

    //s3
    TH2D* hd3 = new TH2D("hd3", ";Row_{2} [MeV];Row_{3} [MeV]", 1000, 0, 120, 1000, 0, 120);
    TH2D* hs3 = new TH2D("hs3", ";Row_{2} [MeV];Row_{3} [MeV]", 1000, 0, 120, 1000, 0, 120);

    auto* cv3 = new TCanvas();
    cv3->SetLogz();
    hd3->GetXaxis()->SetTitleSize(0.04);
    hd3->GetYaxis()->SetTitleSize(0.04);
    hd3->GetXaxis()->SetLabelSize(0.04);
    hd3->GetYaxis()->SetLabelSize(0.04);
    //hd3->GetXaxis()->SetLabelOffset(0.9);
    //hd3->GetYaxis()->SetLabelOffset(0.9);
    hs3->SetMarkerStyle(8);
    hs3->SetMarkerSize(0.7);
    hs3->SetMarkerColorAlpha(kRed, 0.03);
    td->Draw("VamosData.IC[3]*1.04-5:VamosData.IC[2]*1.072-6.4>>hd3", "", "col");
    ts->Draw("IC3:IC2>>hs3", "", "same");

    //s3
    TH2D* hd4 = new TH2D("hd4", ";Row_{3} [MeV];Row_{4} [MeV]", 1000, 0, 120, 1000, 0, 120);
    TH2D* hs4 = new TH2D("hs4", ";Row_{3} [MeV];Row_{4} [MeV]", 1000, 0, 120, 1000, 0, 120);

    auto* cv4 = new TCanvas();
    cv4->SetLogz();
    hd4->GetXaxis()->SetTitleSize(0.04);
    hd4->GetYaxis()->SetTitleSize(0.04);
    hd4->GetXaxis()->SetLabelSize(0.04);
    hd4->GetYaxis()->SetLabelSize(0.04);
    //hd4->GetXaxis()->SetLabelOffset(0.9);
    //hd4->GetYaxis()->SetLabelOffset(0.9);
    hs4->SetMarkerStyle(8);
    hs4->SetMarkerSize(0.7);
    hs4->SetMarkerColorAlpha(kRed, 0.03);
    td->Draw("VamosData.IC[4]*0.8-4.5:VamosData.IC[3]*1.04-5>>hd4", "", "col");
    ts->Draw("IC4:IC3>>hs4", "", "same");
}
