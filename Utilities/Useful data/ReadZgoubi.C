struct VamosFragment{
    double th;
    double ph;
    double delta;
    double xf;
    double thf;
    double yf;
    double phf;
    double L;
    int chState;
    VamosFragment():
        th(-10000),
        ph(-10000),
        delta(-10000),
        xf(-10000),
        thf(-10000),
        yf(-10000),
        phf(-10000),
        L(-10000),
        chState(-10000){};
};


void ReadZgoubi(){
    std::ifstream file_vamos_acceptance("Zgoubi_output_CRYOTARGET_phi.dat");
    std::string vamos_accpetance_line;
    VamosFragment vamosFragment;
    TFile* f = new TFile("map.root", "recreate");
    TTree* tr = new TTree("acc", "acc");
    tr->Branch("th",    &(vamosFragment.th   ), "th/D");
    tr->Branch("ph",    &(vamosFragment.ph   ), "ph/D");
    tr->Branch("delta", &(vamosFragment.delta), "delta/D");
    tr->Branch("xf",    &(vamosFragment.xf   ), "xf/D");
    tr->Branch("thf",   &(vamosFragment.thf  ), "thf/D");
    tr->Branch("yf",    &(vamosFragment.yf   ), "yf/D");
    tr->Branch("phf",   &(vamosFragment.phf  ), "phf/D");
    tr->Branch("L",     &(vamosFragment.L    ), "L/D");

    while(std::getline(file_vamos_acceptance, vamos_accpetance_line)){
        if (vamos_accpetance_line.compare(0, 1, "/") == 0) continue;
        sscanf(vamos_accpetance_line.c_str(), "%lf %lf %lf %lf %lf %lf %lf %lf",
                &vamosFragment.th,
                &vamosFragment.ph,
                &vamosFragment.delta,
                &vamosFragment.xf,
                &vamosFragment.thf,
                &vamosFragment.yf,
                &vamosFragment.phf,
                &vamosFragment.L
              );    
        

        vamosFragment.th    *= -1;
        vamosFragment.ph    *= -1;
        vamosFragment.delta *=  1;
        vamosFragment.xf    *=-10;
        vamosFragment.thf   *= -1;
        vamosFragment.yf    *=-10;
        vamosFragment.phf   *= -1;
        vamosFragment.L     *=  1;

        //printf("%lf %lf %lf %lf %lf %lf %lf %lf\n",
        //        vamosFragment.th,
        //        vamosFragment.ph,
        //        vamosFragment.delta,
        //        vamosFragment.xf,
        //        vamosFragment.thf,
        //        vamosFragment.yf,
        //        vamosFragment.phf,
        //        vamosFragment.L
        //      );    

        tr->Fill();
    }
    tr->Write();

    gStyle->SetOptStat(kFALSE);

    TH2D* h= new TH2D("h", ";Theta_{t} [mrad];Delta",41, -100, 100, 1401, 0.7, 1.4);
    auto* cv = new TCanvas();
    tr->Draw("delta:th>>h", "", "col");
    h->GetXaxis()->SetTitleSize(0.05);
    h->GetYaxis()->SetTitleSize(0.05);
    h->GetXaxis()->SetLabelSize(0.05);
    h->GetYaxis()->SetLabelSize(0.05);
    h->GetXaxis()->SetTitleOffset(0.9);
    h->GetYaxis()->SetTitleOffset(0.9);
    cv->SaveAs("vamos_theta_delta.pdf");

    TH2D* h1= new TH2D("h1", ";Phi_{t} [mrad];Delta",41, -100, 100, 1401, 0.7, 1.4);
    auto* cv1 = new TCanvas();
    tr->Draw("delta:ph>>h1", "", "col");
    h1->GetXaxis()->SetTitleSize(0.05);
    h1->GetYaxis()->SetTitleSize(0.05);
    h1->GetXaxis()->SetLabelSize(0.05);
    h1->GetYaxis()->SetLabelSize(0.05);
    h1->GetXaxis()->SetTitleOffset(0.9);
    h1->GetYaxis()->SetTitleOffset(0.9);
    cv1->SaveAs("vamos_phi_delta.pdf");


    TH2D* h2= new TH2D("h2", ";Theta_{t} [mrad];Phi_{t} [mrad]",41, -100, 100, 41, -100, 100);
    auto* cv2 = new TCanvas();
    tr->Draw("ph:th>>h2", "", "col");
    h2->GetXaxis()->SetTitleSize(0.05);
    h2->GetYaxis()->SetTitleSize(0.05);
    h2->GetXaxis()->SetLabelSize(0.05);
    h2->GetYaxis()->SetLabelSize(0.05);
    h2->GetXaxis()->SetTitleOffset(0.9);
    h2->GetYaxis()->SetTitleOffset(0.9);
    cv2->SaveAs("vamos_theta_phi.pdf");

    TH3D* h3= new TH3D("h3", ";Delta;Phi_{t} [mrad];Theta_{t} [mrad]",1400, 0.69, 1.41, 41, -100, 100, 41, -100, 100);
    h3->SetFillColor(kRed+3);
    auto* cv3 = new TCanvas();
    tr->Draw("th:ph:delta>>h3", "", "gliso");
    h3->GetXaxis()->SetTitleSize(0.05);
    h3->GetYaxis()->SetTitleSize(0.05);
    h3->GetXaxis()->SetLabelSize(0.05);
    h3->GetYaxis()->SetLabelSize(0.05);
    h3->GetXaxis()->SetTitleOffset(0.9);
    h3->GetYaxis()->SetTitleOffset(0.9);
    cv3->SaveAs("vamos_theta_phi_delta.pdf");

    f->Write();
    f->Close();
}
