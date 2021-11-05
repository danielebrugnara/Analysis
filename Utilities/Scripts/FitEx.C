#include <TFitResultPtr.h>
#include <TFitResult.h>

void FitEx(){
    double min{-10};
    double max{10};
    int nbins{70};

    TFile* f = new TFile("../../DataAnalyzed/sum.root", "read");
    TH1D* h = new TH1D("h", "h", nbins, min, max);
    
    TTree* t = (TTree*) f->Get("AnalyzedTreeDeuterons");

    t->Draw("MugastData.Ex>>h", "(Time<313 || Time>315) && (Time<278 || Time>283) && (Time<315 || Time>318)", "");

    double minfit{-2};
    double maxfit{5};
    TF1* fun = new TF1("fun", "gaus(0)+gaus(3)",minfit, maxfit);

    fun->SetParameter(0, 40);
    fun->SetParameter(1, 0);
    fun->SetParLimits(1, 0, 0.5);
    fun->SetParameter(2, 1.1);

    fun->SetParameter(3, 20);
    fun->SetParameter(4, 2);
    fun->SetParLimits(4, 1.8, 2.2);
    fun->SetParameter(5, 1.1);
 
    
    h->Sumw2();
    TFitResultPtr resptr = h->Fit(fun, "SER");
    TFitResult res = *resptr;


    TF1* fun0 = new TF1("fun0", "gaus(0)",minfit, maxfit);
    TF1* fun3 = new TF1("fun3", "gaus(0)",minfit, maxfit);


    fun0->FixParameter(0, fun->GetParameter(0));
    fun0->FixParameter(1, fun->GetParameter(1));
    fun0->FixParameter(2, fun->GetParameter(2));
    fun0->SetLineColor(kGreen);
    fun0->SetLineStyle(3);

    fun3->FixParameter(0, fun->GetParameter(0+3));
    fun3->FixParameter(1, fun->GetParameter(1+3));
    fun3->FixParameter(2, fun->GetParameter(2+3));
    fun3->SetLineColor(kBlue);
    fun3->SetLineStyle(3);

    h->Draw();
    fun0->Draw("same");
    fun3->Draw("same");

    new TCanvas();
    ((TGraph*) gMinuit->Contour(50, 0, 3))->Draw();

    double p0{fun->GetParameter(0)};
    double p3{fun->GetParameter(3)};
    std::cout << "e0 cov " << sqrt(res.CovMatrix(0,0)) << std::endl;

    double e{sqrt(pow(p3/pow(p0+p3, 2), 2)*res.CovMatrix(0,0)+pow(p0/pow(p0+p3, 2), 2)*res.CovMatrix(3,3)+((p0-p3)/pow(p0+p3, 3))*res.CovMatrix(0,3))};
    std::cout << "Percent 0+2= " << p0/(p0+p3) <<" pm " << e<< std::endl;
}
