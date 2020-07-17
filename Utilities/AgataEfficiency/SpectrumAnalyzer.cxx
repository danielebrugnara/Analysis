#include "SpectrumAnalyzer.h"

SpectrumAnalyzer::SpectrumAnalyzer(const std::string & file_name):
        effgraph(),
        levelscheme("adoptedLevels152Sm.csv"),
        fit_interval(15.),
        proj_interval(30.),
        fit_counter(0){

    auto* file = new TFile(file_name.c_str());

    gg = *((TH2D*) file->Get("mgamma_gamma"));
    file->Close();

    gamma_gamma = levelscheme.GetGammaGamma();
//    for (const auto & it: gamma_gamma){
//        std::cout <<"~~" << *it.first <<" -> "<< *it.second<< std::endl;
//    }
    Analyze();
}

SpectrumAnalyzer::~SpectrumAnalyzer() = default;

void SpectrumAnalyzer::Analyze() {
    TFile out_file("eff_curve.root", "recreate");
    int nponts=0;
    double Br_threashold = 0.01;
    double Integral_threashold1 = 1000;
    double Integral_threashold2 = 10;
    int Entries_threashold = 1E3;
    TH1D projx  = *gg.ProjectionX();
    for (const auto & it: gamma_gamma){
        if (it.first->Br.first < Br_threashold || it.second->Br.first < Br_threashold)
            continue;
        std::cout <<"~~" << *it.first <<" -> "<< *it.second<< std::endl;
        double engamma1 = it.first->Egamma.first/UNITS::keV;
        double engamma2 = it.second->Egamma.first/UNITS::keV;
        double branching2 = it.second->Br.first;

        std::cout << "En1 : " << engamma1<< "\n";
        std::cout << "En2: " << engamma2<< "\n";

        TH1D projy = *gg.ProjectionY("_py",
                                    gg.GetXaxis()->FindBin(engamma1-proj_interval),
                                    gg.GetXaxis()->FindBin(engamma1-proj_interval));
        std::cout << "Entries projy : " << projy.GetEntries()<< "\n";

        if (projy.GetEntries()<Entries_threashold)
            continue;

        double integral1 = GetPeakIntegral(projx, engamma1);
        if (integral1 < Integral_threashold1)
            continue;
        double integral2 = GetPeakIntegral(projy, engamma2);
        if (integral2 < Integral_threashold2)
            continue;

        double eff = integral2/(integral1*branching2);
        effgraph.SetPoint(nponts++,engamma2, eff);

    }

    effgraph.SetName("effgraph");
    effgraph.Write();
    out_file.Write();
    out_file.Close();
}

double SpectrumAnalyzer::GetPeakIntegral(TH1D& spec, const double& engamma){
    UpdateErrors(spec);
    TH1D bkgx   = *((TH1D*)spec.ShowBackground(10));
    UpdateErrors(bkgx);
    TH1D subtrx = spec - bkgx;
    fit_counter++;
    TF1 fitfunc(Form("gaussian_%f_%i", engamma, fit_counter),
                [](const double* x, const double*par)
                {return par[0]*exp(-pow((x[0]-par[1])/(par[2]),2));},
                engamma-fit_interval,
                engamma+fit_interval,
                3,
                1);

    fitfunc.SetParameter(0, spec.GetBinContent(spec.GetXaxis()->FindBin(engamma)));
    fitfunc.SetParameter(1, engamma);
    fitfunc.SetParameter(2, 5);

    fitfunc.SetParLimits(0, fitfunc.GetParameter(0)/3, fitfunc.GetParameter(0)*3);
    fitfunc.SetParLimits(1, engamma-10, engamma+10);
    fitfunc.SetParLimits(2, 0.8, 6);

    fitfunc.FixParameter(1, engamma);

    TFitResultPtr fit_res_ptr(0);

    fit_res_ptr = subtrx.Fit( &fitfunc,
                             "S",
                             "",
                             fitfunc.GetXmin(),
                             fitfunc.GetXmax());
    if (fit_res_ptr>=0){
        fit_res_ptr->Print();
    }else{
        return -1;
    }

    if (fit_res_ptr->Chi2()>5*fit_res_ptr->Ndf())
        return -1;

    TFitResult fit_result = *fit_res_ptr;

    double integral =  fitfunc.GetParameter(0)*sqrt(UNITS::CONSTANTS::pi)*fitfunc.GetParameter(2);
    subtrx.SetName(Form("%s_%f_%i__%f",spec.GetName(), engamma, fit_counter, integral));
    subtrx.Write();

    return integral;
}

void SpectrumAnalyzer::UpdateErrors(TH1D & spec) {
    for (int ii=0; ii<spec.GetNbinsX(); ++ii){
        spec.SetBinError(ii, sqrt(spec.GetBinContent(ii)));
    }
}

