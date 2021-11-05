void PlotEx_corrected(){
    gStyle->SetOptFit(11111111);
    auto* f = new TFile("../../build/Out/sum.root", "read");
    auto* t = (TTree*)f->Get("AnalyzedTreeDeuterons");

    auto* cv = new TCanvas();
    cv->Divide(6,6);

    //std::string conditions{"MugastData.MG == 11"};
    std::string conditions{""};
    for(int i{0}; i<32; ++i){
        cv->cd(i+1);
        std::cout << "spectrum nr : " << i << std::endl;
        auto* h = new TH1D(Form("h%i", i),Form("h%i", i), 50, -10, 6);
        t->Draw(Form("MugastData.Ex_corrected[0][%i]>>h%i", i, i),conditions.c_str() );
        //t->Draw(Form("MugastData.Ex_corrected[0][%i]:Time>>h%i", i, i));
        //t->Draw(Form("MugastData.Ex_corrected[0][%i]:CatsData.Pos.x()>>h%i", i, i));
        auto* fun = new TF1(Form("f%i", i), "gaus(0)+gaus(3)", -4, 6);
        //mean
        //fun->FixParameter(1, 0);
        //fun->FixParameter(4, 2.02);
        fun->SetParameter(1, 0);
        fun->SetParameter(4, 2.02);
        fun->SetParLimits(1, 0, 0.1);
        fun->SetParLimits(4, 1.8, 2.4);
        //ampl
        fun->SetParameter(0, 200);
        fun->SetParameter(3, 200);
        //sigma
        fun->SetParameter(2, 1.1);
        fun->SetParameter(5, 1.1);
        h->Fit(fun, "RL");
    }
    cv->cd(33);
    auto* h = new TH1D(Form("h"),Form("h"), 50, -10, 6);
    t->Draw(Form("MugastData.Ex>>h"), conditions.c_str());
    auto* fun = new TF1(Form("f"), "gaus(0)+gaus(3)", -4, 6);
    //mean
    //fun->FixParameter(1, 0);
    //fun->FixParameter(4, 2.02);
    fun->SetParameter(1, 0);
    fun->SetParameter(4, 2.02);
    fun->SetParLimits(1, 0, 0.1);
    fun->SetParLimits(4, 1.8, 2.4);
    //ampl
    fun->SetParameter(0, 200);
    fun->SetParameter(3, 200);
    //sigma
    fun->SetParameter(2, 1.1);
    fun->SetParameter(5, 1.1);
    h->Fit(fun, "RL");
}
