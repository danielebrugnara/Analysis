void PlotEx_corrected(){
    gStyle->SetOptFit(1111);
    auto* f = new TFile("../../DataAnalyzed/sum_100grid.root", "read");
    auto* t = (TTree*)f->Get("AnalyzedTreeDeuterons");
    auto* probGraph = new TGraph2D();
    bool useCorrected{true};

    auto* cv = new TCanvas();
    cv->Divide(10,10);
    double start{-1.1};
    double end{1.1};
    if(!useCorrected){
        start = -5;    
        end   =  5; 
    }
    double pace = (end - start)/(10-1);

    //std::string conditions{"MugastData.MG == 11"};
    std::pair<int,double> prob{std::make_pair(0, 0)};
    TH1D* bestHisto;
    std::string conditions{""};
    //std::string conditions{"(Time<313 || Time>315) && (Time<278 || Time>283) && (Time<315 || Time>318)"};
    for(int i{0}; i<100; ++i){
        cv->cd(i+1);
        std::cout << "spectrum nr : " << i << std::endl;
        TH1D* h;
        if(useCorrected){
            h = new TH1D(Form("h%i", i),Form("h%i", i), 70, -10, 6);
            t->Draw(Form("MugastData.Ex_corrected[0][%i]>>h%i", i, i),conditions.c_str() );

        }else{
            h = new TH1D(Form("h%i", i),Form("h%i", i), 100, -10, 6);
            t->Draw(Form("MugastData.Ex_uncentered[0][%i]>>h%i", i, i),conditions.c_str() );
        }
        //t->Draw(Form("MugastData.Ex_corrected[0][%i]:Time>>h%i", i, i));
        //t->Draw(Form("MugastData.Ex_corrected[0][%i]:CatsData.Pos.x()>>h%i", i, i));
        auto* fun = new TF1(Form("f%i", i), "gaus(0)+gaus(3)", -4, 6);
        //mean
        fun->FixParameter(1, 0);
        fun->FixParameter(4, 2.02);
        //fun->SetParameter(1, 0);
        //fun->SetParameter(4, 2.02);
        fun->SetParLimits(1, 0, 0.1);
        fun->SetParLimits(4, 1.8, 2.4);
        //ampl
        fun->SetParameter(0, 200);
        fun->SetParameter(3, 200);
        //sigma
        fun->SetParameter(2, 1.2);
        fun->SetParameter(5, 1.2);
        fun->SetParLimits(2, 1, 2);
        fun->SetParLimits(5, 1, 2);
        auto res = h->Fit(fun, "RLS");
        probGraph->AddPoint(start+(i%10)*pace, start+(i/10)*pace, fun->GetProb());
        if (fun->GetProb() > prob.second){
            prob.first = i;
            prob.second = fun->GetProb();
            bestHisto = h;
        }
        
    }
    std::cout << "best is : " << prob.first << std::endl;
    bestHisto->SetLineColor(kRed+2);
    bestHisto->SetLineWidth(4);
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


    new TCanvas();
    if(useCorrected){
        probGraph->SetTitle(";Coefficient X; Coefficient Y");
    }else{
        probGraph->SetTitle(";Position X [mm]; Position Y [mm]");
    }
    probGraph->GetXaxis()->SetLabelSize(0.05);
    probGraph->GetXaxis()->SetTitleSize(0.05);
    probGraph->GetXaxis()->SetTitleOffset(0.9);
    probGraph->GetYaxis()->SetLabelSize(0.05);
    probGraph->GetYaxis()->SetTitleSize(0.05);
    probGraph->GetYaxis()->SetTitleOffset(0.9);
    probGraph->SetLineColorAlpha(kBlack, 0.26);
    probGraph->SetLineStyle(5);
    probGraph->SetLineWidth(3);
    
    probGraph->Draw("colz");
    probGraph->Draw("cont3 same");
}
