void PlotEx(){
    gStyle->SetOptStat(kFALSE);
    gStyle->SetOptFit(1111);
    auto* f= new TFile("../../DataAnalyzed/sum.root", "read");
    auto* t= (TTree*)f->Get("AnalyzedTreeDeuterons");

    int nbins{100};
    double start{-10};
    double stop{10};
    auto* h= new TH1D("h", Form(";Excitation energy [MeV];Counts/%.2f MeV", (stop-start)/nbins), nbins, start, stop);
    h->SetMarkerStyle(8);
    h->SetMarkerSize(0.5);
    h->SetMarkerColorAlpha(kBlack, 0.7);
    h->SetLineColorAlpha(kBlack, 0.7);
    h->GetXaxis()->SetTitleSize(0.05);
    h->GetXaxis()->SetLabelSize(0.05);
    h->GetXaxis()->SetTitleOffset(0.9);

    h->GetYaxis()->SetTitleSize(0.05);
    h->GetYaxis()->SetLabelSize(0.05);
    h->GetYaxis()->SetTitleOffset(0.9);


    auto* fun = new TF1("f", "gaus(0)+gaus(3)", -3, 5);
    auto* g1 = new TF1("g1", "gaus(0)", -3, 5);
    auto* g2 = new TF1("g2", "gaus(0)", -3, 5);
    //mean
    fun->SetParameter(1, 0.25);
    fun->SetParameter(4, 2.02);
    //fun->SetParameter(1, 0);
    //fun->SetParameter(4, 2.02);
    fun->SetParLimits(1, 0, 0.5);
    fun->SetParLimits(4, 1.8, 2.6);
    //ampl
    fun->SetParameter(0, 30);
    fun->SetParameter(3, 20);
    //fun->SetParameter(0, 40);
    //fun->SetParameter(3, 20);
    fun->SetParLimits(0, 5, 50);
    fun->SetParLimits(3, 5, 50);
    //sigma
    //
    fun->SetParameter(2, 1.2);
    fun->SetParameter(5, 1.2);
    //fun->SetParameter(2, 1.2);
    //fun->SetParameter(5, 1.2);
    fun->SetParLimits(2, 0.5, 1.3);
    fun->SetParLimits(5, 0.5, 1.3);

    new TCanvas();
    t->Draw("MugastData.Ex>>h", "", "e");     
    auto res = h->Fit(fun, "RLS");

    for(int i{0};i<3;++i){
        g1->SetParameter(i, fun->GetParameter(i));
        g2->SetParameter(i, fun->GetParameter(i+3));
    }

    fun->Draw("same");
    g1->Draw("same");
    g2->Draw("same");
}
