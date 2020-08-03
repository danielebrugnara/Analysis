#include "SpectrumAnalyzer.h"

SpectrumAnalyzer::SpectrumAnalyzer(const std::string & file_name):
        effgraph(),
        levelscheme("adoptedLevels152Sm.csv"),
        fit_interval(15.),
        proj_interval(30.),
        fit_counter(0){

    eu152_intensities = ReadIntensities("intensity152Eu.txt");

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
    TH1D* subtrprojx = Subtract(projx);
    UpdateErrors(*subtrprojx);
    //Plotter pltx(projx, "");
    int cnt{0};
    std::vector<int> meaningful_combinations = {
            1, 13, 19, 21, 26, 27, 50,
            66, 67, 68, 69, 70, 72, 79,
            81, 159
    };


    for (const auto & it: gamma_gamma){
        if (it.first->Br.first < Br_threashold || it.second->Br.first < Br_threashold)
            continue;
        if (std::find(meaningful_combinations.begin(), meaningful_combinations.end(), cnt++)!=meaningful_combinations.end())
            continue;

        std::cout <<"~~" << *it.first <<" -> "<< *it.second<< std::endl;
        double engamma1 = it.first->Egamma.first/UNITS::keV;
        double engamma2 = it.second->Egamma.first/UNITS::keV;
        double branching2 = it.second->Br.first;

        std::cout << "En1 : " << engamma1<< "\n";
        std::cout << "En2: " << engamma2<< "\n";

        double integral1 = GetPeakIntegral(*subtrprojx, engamma1);
        if (integral1 < Integral_threashold1){
            continue;
        }

        double gate_width = 5.;
        double noise_gate_width = 50.;
        double noise_gate_dist = 10.;

        TH1D* projy = ProjectAndSubtract(gg,
                                         engamma1-gate_width/2., engamma1+gate_width/2.,
                                         engamma1-noise_gate_dist-noise_gate_width, engamma1-noise_gate_dist,
                                         engamma1+noise_gate_dist, engamma1+noise_gate_dist+noise_gate_width);
        std::cout << "Entries projy : " << projy->GetEntries()<< "\n";

        projy = Subtract(*projy);

        if (projy->GetEntries()<Entries_threashold) {
            delete projy;
            continue;
        }
        double integral2 = GetPeakIntegral(*projy, engamma2);

        if (integral2 < Integral_threashold2){
            delete projy;
            continue;
        }


        double eff = integral2/(integral1*branching2);
        canvas->SetTitle(Form("comb%i", cnt));
        plotter.WriteOnCanvas(std::string("int1 : ")+
                              std::to_string(integral1)+
                              " int2 : "+
                              std::to_string(integral2)+
                              " eff: "+
                              std::to_string(eff),
                              46);

        effgraph.SetPoint(nponts++,engamma2, eff);
        delete projy;

    }

    effgraph.SetName("effgraph");
    effgraph.Write();
    out_file.Write();
    out_file.Close();
}

double SpectrumAnalyzer::GetPeakIntegral(TH1D& spec, const double& engamma){
    plotter.PlotOnCanvas<TH1>(spec, "histo");
    UpdateErrors(spec);
    fit_counter++;
    TF1 fitfunc(Form("gaussian_%f_%i", engamma, fit_counter),
                [](const double* x, const double*par)
                {return par[0]*exp(-pow((x[0]-par[1])/(par[2]),2));},
                engamma-fit_interval,
                engamma+fit_interval,
                3,
                1);
    fitfunc.SetNpx(400);

    fitfunc.SetParameter(0, spec.GetBinContent(spec.GetXaxis()->FindBin(engamma)));
    fitfunc.SetParameter(1, engamma);
    fitfunc.SetParameter(2, 4);

    fitfunc.SetParLimits(0, fitfunc.GetParameter(0)/3, fitfunc.GetParameter(0)*3);
    fitfunc.SetParLimits(1, engamma-10, engamma+10);
    fitfunc.SetParLimits(2, 0.8, 5);

    fitfunc.FixParameter(1, engamma);

    TFitResultPtr fit_res_ptr(0);

    fit_res_ptr = spec.Fit( &fitfunc,
                              "S",
                              "",
                              fitfunc.GetXmin(),
                              fitfunc.GetXmax());
    plotter.PlotOnCanvas<TF1>(fitfunc, "same");
    if (fit_res_ptr>=0){
        fit_res_ptr->Print();
    }else{
        plotter.WriteOnCanvas("FitNotSuccessful!!");
        return -1;
    }

//    if (fit_res_ptr->Chi2()>5*fit_res_ptr->Ndf()) {
//        plotter.WriteOnCanvas(std::string("high chi : ")+
//                                            std::to_string(fit_res_ptr->Chi2())+
//                                            " with ndf : "+
//                                            std::to_string(fit_res_ptr->Ndf()));
//        return -1;
//    }

    TFitResult fit_result = *fit_res_ptr;

    double integral =  fitfunc.GetParameter(0)*sqrt(UNITS::CONSTANTS::pi)*fitfunc.GetParameter(2);
    plotter.WriteOnCanvas(std::string("correct chi : ")+
                          std::to_string(fit_res_ptr->Chi2())+
                          " with ndf : "+
                          std::to_string(fit_res_ptr->Ndf())+
                          ", integral : "+
                          std::to_string(integral));

    // plotter.SetRange(engamma-100, engamma+100);

    //subtrx.SetName(Form("%s_%f_%i__%f",spec.GetName(), engamma, fit_counter, integral));
    //subtrx.Write();

    return integral;
}

void SpectrumAnalyzer::UpdateErrors(TH1D & spec) {
    for (int ii=0; ii<spec.GetNbinsX(); ++ii){
        spec.SetBinError(ii, sqrt(spec.GetBinContent(ii)));
    }
}

TH1D* SpectrumAnalyzer::ProjectAndSubtract(  const TH2D& mat,
                                             const double& e_min_sig  ,const double& e_max_sig,
                                             const double& e_min_left ,const double& e_max_left,
                                             const double& e_min_right,const double& e_max_right){

    // Get the bin width only for Y axis information
    float bin_width = mat.GetYaxis()->GetBinWidth(100);
    // To define the new histograms
    Int_t nb_bin_mat = mat.GetYaxis()->GetNbins();
    Float_t range_min = mat.GetYaxis()->GetXmin();
    Float_t range_max = mat.GetYaxis()->GetXmax();

    //TH2 *mat_gg = (TH2D*)afile->Get(Form("gamma_sig_ring%i_%i",j,k));
    TH1D* signal = new TH1D("signal"  ,"sig"     ,nb_bin_mat,range_min,range_max);
    TH1D noise    ("noise"   ,"noise"   ,nb_bin_mat,range_min,range_max);
    TH1D proj_sig ("proj_sig","proj_sig",nb_bin_mat,range_min,range_max);
    TH1D proj_nog ("proj_nog","proj_nog",nb_bin_mat,range_min,range_max);
    TH1D proj_nod ("proj_nod","proj_nod",nb_bin_mat,range_min,range_max);



    mat.ProjectionY("proj_sig",
                    mat.GetXaxis()->FindBin(e_min_sig),
                    mat.GetXaxis()->FindBin(e_max_sig));

    mat.ProjectionY("proj_nog",
                    mat.GetXaxis()->FindBin(e_min_left),
                    mat.GetXaxis()->FindBin(e_max_left));

    mat.ProjectionY("proj_nod",
                    mat.GetXaxis()->FindBin(e_min_right),
                    mat.GetXaxis()->FindBin(e_max_right));

//  std::cout << "Signal: " <<bin_sig << "\t" << bin_sig+nb_bin << "\n"
//	    << "NoG: "    <<bin_nog << "\t" << bin_nog+nb_bin <<"\n"
//	    << "NoD: "    <<bin_nod << "\t" << bin_nod+nb_bin << std::endl;
//

    Float_t egamma = (e_min_sig+e_max_sig)/2.;

    float nb_bin_sig   = mat.GetXaxis()->FindBin(e_max_sig)-mat.GetXaxis()->FindBin(e_min_sig);
    float nb_bin_noise = mat.GetXaxis()->FindBin(e_max_left)-mat.GetXaxis()->FindBin(e_min_left)
                         +mat.GetXaxis()->FindBin(e_max_right)-mat.GetXaxis()->FindBin(e_min_right);
    noise.Sumw2();
    noise.Add(&proj_nog,&proj_nod,1,1);

    signal->Sumw2();
    signal->Add(&proj_sig,1);
    signal->Add(&noise,-nb_bin_sig/nb_bin_noise);

//  signal->Draw("hist");
    signal->SetTitle(Form("Gate on %.2f [%.2f,%.2f]",egamma,e_min_sig,e_max_sig));
    signal->GetXaxis()->SetTitle("E_{#gamma} [keV]");
    signal->GetXaxis()->CenterTitle(false);
    signal->GetXaxis()->SetTitleSize(0.05);
    signal->GetXaxis()->SetTitleOffset(0.80);
    signal->GetYaxis()->SetTitle(Form("Counts / %.2f keV",bin_width));
    signal->GetYaxis()->SetTitleSize(0.05);
    signal->GetYaxis()->CenterTitle(false);
    proj_sig.Delete();
    proj_nog.Delete();
    proj_nod.Delete();

    return signal;
}

TH1D* SpectrumAnalyzer::Subtract(TH1D & spec) {
    UpdateErrors(spec);
    TH1D bkgx   = *((TH1D*)spec.ShowBackground(10, ""));
    UpdateErrors(bkgx);

    TH1D* subtrx = new TH1D(spec);
    //subtrx->Sumw2();
    subtrx->Add(&bkgx, -1);
    return subtrx;
}

SpectrumAnalyzer::IntensityData SpectrumAnalyzer::ReadIntensities(const std::string & file_name) {
    std::ifstream file(file_name);
    std::string line;
    int cnt{0};
    std::vector<std::pair<double, double>> data; //Energy, intensity
    double energy, intensity;
    while (std::getline(file, line))
    {
        ++cnt;
        if (cnt<3) continue;
        std::istringstream ss(line);
        std::string substr;
        getline(ss, substr, '\t');

        //std::size_t found = substr.find(' ');
        //if (found !=std::string::npos)
        //    substr = substr.substr(0, found);

        energy = std::stod(substr);
        getline(ss, substr, '\t');
        intensity = std::stod(substr);
        data.emplace_back(energy, intensity);

    }
    return data;
}

TGraph SpectrumAnalyzer::GenerateRelativeEffGraph() {
    TH1D projx  = *gg.ProjectionX();
    TH1D* subtrprojx = Subtract(projx);
    UpdateErrors(*subtrprojx);
    for (const auto& it_eff: eu152_intensities){
        std::cout << "hi";
    }
}

