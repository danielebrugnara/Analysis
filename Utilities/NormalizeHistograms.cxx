#include <TH1D.h>

TH1D* NormalizeHistograms(TH1D* a, TH1D* b){
    TH1D* tmp = new TH1D("tmp", "tmp", 
                            a->GetXaxis()->GetNbins(),
                            a->GetXaxis()->GetXmin(),
                            a->GetXaxis()->GetXmax());
    for(int ii=0; ii<a->GetXaxis()->GetNbins();++ii){
        if (b->GetBinContent(ii)==0) continue;
        double conta = a->GetBinContent(ii);
        double contb = b->GetBinContent(ii);
        tmp->SetBinContent(ii,conta/contb);
        tmp->SetBinError(ii, sqrt(pow(sqrt(conta)/contb,2)+pow(conta*sqrt(contb)/pow(contb,2),2)));
    }
    return tmp;
}

TH1D* GetDistribution(TH1D * a){
    TH1D* tmp = new TH1D("tmp", "tmp", 
                            a->GetXaxis()->GetNbins(),
                            a->GetXaxis()->GetXmin(),
                            a->GetXaxis()->GetXmax());
    for(int ii=0; ii<a->GetXaxis()->GetNbins();++ii){
        if (a->GetBinContent(ii)==0) continue;
        tmp->SetBinContent(ii,a->GetBinContent(ii)/sin(a->GetBinCenter(ii)));
        tmp->SetBinError(ii, a->GetBinError(ii)/sin(a->GetBinCenter(ii)));
    }
    return tmp;
}
