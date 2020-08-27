#include "SpectrumAnalyzer.h"

SpectrumAnalyzer::SpectrumAnalyzer(const std::string & file_name, const bool& debug_canvas):
        debug_canvas(debug_canvas),
        effgraph(),
        relative_effgraph(),
        relative_integralgraph(),
        relative_intgraph(),
        levelscheme("adoptedLevels152Sm.csv"),
        fit_interval(15.),
        proj_interval(30.),
        fit_counter(0){

    eu152_intensities = ReadIntensities("intensity152Eu.txt");

    auto* file = new TFile(file_name.c_str());

    auto gg_ptr = dynamic_cast<TH2D*>(file->Get("mgamma_gamma"));
    auto hspec_ptr = dynamic_cast<TH1D*>(file->Get("hspec"));
    start_stop = *dynamic_cast<TVector2*>(file->Get("start_stop"));

    if (gg_ptr == nullptr)
        throw std::runtime_error("mgamma_gamma  not found\n");

    gg = *gg_ptr;
    hspec = *hspec_ptr;
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

    GenerateRelativeEffGraph();
    double acq_time = start_stop.Y()-start_stop.X();
    double source_activity = 2.2E4;
    scaled_relative_effgraph.SetName("scaled_relative_effgraph");
    scaled_relative_effgraph.SetTitle("scaled_relative_effgraph");
    scaled_relative_effgraph.SetMarkerColor(8);
    scaled_relative_effgraph.SetMarkerStyle(21);
    for (int ii=0; ii<relative_effgraph.GetN(); ++ii){
        scaled_relative_effgraph.SetPoint(ii,
                                          relative_effgraph.GetPointX(ii),
                                          relative_effgraph.GetPointY(ii)/(acq_time*source_activity));
        scaled_relative_effgraph.SetPointError(     ii,
                                                    0,
                                                    relative_effgraph.GetErrorY(ii)/(acq_time*source_activity));
    }


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
        if (debug_canvas)
            canvas->SetTitle(Form("comb%i", cnt));
        if (debug_canvas)
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
    relative_effgraph.Write();
    scaled_relative_effgraph.Write();
    relative_integralgraph.Write();
    relative_intgraph.Write();
    effgraph.Write();
    out_file.Write();
    out_file.Close();
}

double SpectrumAnalyzer::GetPeakIntegral(TH1D& spec, const double& engamma){
    if (debug_canvas)
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
    fitfunc.SetParameter(2, 1.8);

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
    if (debug_canvas)
        plotter.PlotOnCanvas<TF1>(fitfunc, "same");
    if (debug_canvas)
        plotter.SetRange(engamma-50, engamma+50);

    if (fit_res_ptr>=0){
        fit_res_ptr->Print();
    }else{
        if (debug_canvas)
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
    if (debug_canvas)
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

std::vector<std::pair<double, double>> SpectrumAnalyzer::GetPeaksIntegral(TH1D &spec, const std::vector<std::pair<double,double>> &energies) {

    spec.GetXaxis()->UnZoom();
    std:: cout << "------------> Fitting new energies : " <<  energies.front().first << "<----------------------------------------------------------\n";
    const double fit_interval = 8;
    std::vector<double> heights;
    const double MIN_VARIANCE{0.8};
    const double MAX_VARIANCE{5.0};
    //if (energies.front().first<1080) return std::vector<std::pair<double, double>>();

    for (const auto& it_en: energies) {
        std:: cout << "===== Single peak fit :" <<  it_en.first << "\n";
        TF1 fitfunc(Form("gaussian_%f_%i", it_en.first, fit_counter),
                    [](const double *x, const double *par) {
                        return par[0] * exp(-pow((x[0] - par[1]) / (par[2]), 2));
                    },
                    it_en.first - fit_interval,
                    it_en.first + fit_interval,
                    3,
                    1);
        fitfunc.SetNpx(400);

        fitfunc.SetParameter(0, std::max(0.,spec.GetBinContent(spec.GetXaxis()->FindBin(it_en.first))));
        fitfunc.SetParameter(1, it_en.first);
        fitfunc.SetParameter(2, 1.8);

        fitfunc.SetParLimits(0, std::max(0.,fitfunc.GetParameter(0) / 3.), std::max(fitfunc.GetParameter(0) * 3., 100.));
        fitfunc.SetParLimits(1, it_en.first - 10, it_en.first + 10);
        fitfunc.SetParLimits(2, MIN_VARIANCE, MAX_VARIANCE);

        fitfunc.FixParameter(1, it_en.first);

        TFitResultPtr fit_res_ptr(0);


        fit_res_ptr = spec.Fit(&fitfunc,
                               "S",
                               "",
                               fitfunc.GetXmin(),
                               fitfunc.GetXmax());

        heights.push_back(fitfunc.GetParameter(0));

    }

    int number_of_gaussians = energies.size();
    const int npars = 4;
    double lim_inf = energies.front().first - fit_interval;
    double lim_sup = energies.back().first + fit_interval;
//    //https://en.wikipedia.org/wiki/Exponentially_modified_Gaussian_distribution
//    TF1 fitfunc(Form("gaussians_%f_%i", energies.front().first, fit_counter),
//                [number_of_gaussians](const double *x, const double *par) {
//                    long double result{0};
//                    for(int ii=0; ii<number_of_gaussians; ++ii){
//                        double argument = 1./sqrt(2)*(par[2+ii*npars]/par[3+ii*npars]+(x[0]-par[1+ii*npars])/par[2+ii*npars]);
//                        result +=   par[0+ii*npars] *
//                                    exp(-1./2.*pow((x[0]-par[1+ii*npars])/par[2+ii*npars],2))*
//                                    0.5/par[3+ii*npars]*
//                                    exp(pow(argument,2))*erfc(argument);
//                    };
//                    return static_cast<double>(result*par[number_of_gaussians*npars]);
//                },
//                lim_inf,
//                lim_sup,
//                npars*number_of_gaussians+1,
//                1);
//    fitfunc.SetNpx(200);

    std::vector<TF1Convolution> tailed_peaks_conv;
    std::vector<TF1> tailed_peaks;
    std::vector<TF1> gaussian_peaks;
    std::vector<TF1> tails;
    tailed_peaks_conv.reserve(number_of_gaussians);
    tailed_peaks.reserve(number_of_gaussians);
    gaussian_peaks.reserve(number_of_gaussians);
    tails.reserve(number_of_gaussians);
    for(int ii=0; ii<number_of_gaussians; ++ii) {
        tails.emplace_back(Form("Tail_%i", ii),
                           [](const double*x, const double*par){
                                if(x[0]<0)
                                    return exp(-x[0]*par[0]);
                                else
                                    return 0.;
                            },
                            lim_inf,
                            lim_sup,
                            1,
                            1);
        gaussian_peaks.emplace_back(Form("Tail_%i", ii),
                                    [](const double*x, const double*par){
                                        return par[0]*exp(-pow((x[0]-par[1])/par[2],2));
                                    },
                                    lim_inf,
                                    lim_sup,
                                    3,
                                    1);
        tailed_peaks_conv.emplace_back( &gaussian_peaks.back(),
                                        &tails.back(),
                                        lim_inf,
                                        lim_sup,
                                        true);
        tailed_peaks_conv.back().SetNofPointsFFT(1000);

        tailed_peaks.emplace_back(Form("Tailed_peak_%i", ii),
                                  tailed_peaks_conv.back(),
                                  lim_sup,
                                  lim_inf,
                                  tailed_peaks_conv.back().GetNpar());
    }
    TF1 fitfunc(Form("gaussians_%f_%i", energies.front().first, fit_counter),
                [number_of_gaussians, &tailed_peaks](const double *x, const double *par) {
                    long double result{0};
                    for(int ii=0; ii<number_of_gaussians; ++ii){
                        result += tailed_peaks[ii].Eval(x[0]);
                    };
                    return static_cast<double>(result*par[number_of_gaussians*npars]);
                },
                lim_inf,
                lim_sup,
                npars*number_of_gaussians+1,
                1);


    for(int ii=0; ii<number_of_gaussians; ++ii) {
        fitfunc.SetParameter(1+npars*ii, energies[ii].first);
        fitfunc.SetParameter(2+npars*ii, sqrt(1.8));
        fitfunc.SetParameter(3+npars*ii, 0.7);

        //fitfunc.SetParLimits(0+3*ii, std::max(1.,heights[ii]*0.1), std::max(heights[ii]*10, 100.));
        fitfunc.SetParLimits(1+3*ii, energies[ii].first - 0.2, energies[ii].first + 0.2 );
        fitfunc.SetParLimits(2+npars*ii, sqrt(MIN_VARIANCE), sqrt(MAX_VARIANCE));
        fitfunc.SetParLimits(3+npars*ii, 0.05, 3*sqrt(MAX_VARIANCE));


        fitfunc.FixParameter(1+npars*ii, energies[ii].first);
        fitfunc.FixParameter(0+npars*ii, energies[ii].second);
    }
    int max_idx = 0;
    for (int ii=0; ii<heights.size(); ++ii){
        if (heights[ii]>heights[max_idx])
            max_idx = ii;
    }

    fitfunc.SetParameter(number_of_gaussians*npars, heights[max_idx]/energies[max_idx].second);
    fitfunc.SetParLimits(number_of_gaussians*npars,
                         std::max(0.,fitfunc.GetParameter(number_of_gaussians*npars)*0.1),
                         std::max(1.,fitfunc.GetParameter(number_of_gaussians*npars)*10.) );

    if (debug_canvas)
        plotter.PlotOnCanvas<TH1D>(spec, "histo");
    TFitResultPtr fit_res_ptr(0);
    std:: cout << "===== Multi fit 1:\n";
    //TVirtualFitter::SetDefaultFitter("Minuit2");
    //ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2");
    fit_res_ptr = spec.Fit( &fitfunc,
                            "SMR",
                            "",
                            fitfunc.GetXmin(),
                            fitfunc.GetXmax());

    //std:: cout << "===== Multi fit 2:\n";
   // fit_res_ptr = spec.Fit( &fitfunc,
   //                         "SMERI",
   //                         "",
   //                         fitfunc.GetXmin(),
   //                         fitfunc.GetXmax());

    TFitResult fit_result = *fit_res_ptr;
    auto cov = fit_result.GetCovarianceMatrix();
    auto aa  = fit_result.IsValid();

    std::vector<std::pair<double, double>> integrals;
    for(int ii=0; ii<number_of_gaussians; ++ii) {
        double integral =   fitfunc.GetParameter(0+npars*ii)*fitfunc.GetParameter(number_of_gaussians*npars);
        double variance = pow(fitfunc.GetParameter(2+npars*ii),2);
        double tail = fitfunc.GetParameter(3+npars*ii);
        if (abs(variance-MAX_VARIANCE)<0.02 || abs(variance-MIN_VARIANCE)<0.02)
            integral = 0;
        double error =  fitfunc.GetParameter(0+npars*ii)*fitfunc.GetParError(number_of_gaussians*npars);
        integrals.emplace_back(integral, error);
    }

    if (debug_canvas) {
        plotter.PlotOnCanvas<TF1>(fitfunc, "same");
        plotter.SetRange(fitfunc.GetXmin() , fitfunc.GetXmax() );
        std::string integral_string;
        std::string sigma_string;
        for(int ii=0;ii<integrals.size();++ii){
            integral_string+=std::string("en :")+
                             std::to_string(energies[ii].first)+
                             " int-> "+
                             std::to_string(integrals[ii].first)+"; ";
            sigma_string+=std::string("en :")+
                          std::to_string(energies[ii].first)+
                          " var-> "+
                          std::to_string(fitfunc.GetParameter(2+npars*ii))+"; ";
        }
        plotter.WriteOnCanvas(integral_string,
                              46);
        plotter.WriteOnCanvas(sigma_string,
                              30,
                              "tl");
    }

    return integrals;
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
    //UpdateErrors(spec);
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

        energy = std::stod(substr);
        getline(ss, substr, '\t');
        intensity = std::stod(substr);
        data.emplace_back(energy, intensity*1E-2);

    }
    return data;
}

void SpectrumAnalyzer::GenerateRelativeEffGraph() {
    //TH1D projx  = *gg.ProjectionX();
    TH1D projx = hspec;
    UpdateErrors(projx);
    TH1D* subtrprojx = Subtract(projx);
    //UpdateErrors(*subtrprojx);

    std::vector<std::vector<int>> seen_transitions;
    double near_peak_threash = 6.5;
    for (int ii=0; ii<eu152_intensities.size(); ++ii) {
        if(eu152_intensities[ii].second < 0.001) continue;
        std::vector<int> tmp_vec;
        tmp_vec.push_back(ii);
        while(ii<eu152_intensities.size() && abs(eu152_intensities[ii+1].first-eu152_intensities[ii].first)<near_peak_threash){
            tmp_vec.push_back(++ii);
        }
        seen_transitions.push_back(tmp_vec);
    }

    std::vector<double> X, Y, Xerr, Yerr, Y2, Y3, Y3err;
    //std::vector<double> blacklisted_energies;
    std::vector<double> blacklisted_energies = {
            //125.69, 443.96, 494.0, 493.508, 571.6, 586.264, 719.349,764.9,768.944,919.33, 926.317, 930.58, 937.05, 1457.643
            919.33
    };

    for (const auto& it_transition: seen_transitions){
        std::vector<std::pair<double, double>> energies;
        energies.reserve(it_transition.size());
        for (const auto & it: it_transition){
            energies.push_back(eu152_intensities[it]);
        }
        std::vector<std::pair<double, double>> counts = GetPeaksIntegral(*subtrprojx, energies);
        //counts /= it_eff.second;
        for (int ii=0; ii<counts.size(); ++ii){
            //if (counts[ii].first<2E3) continue;
            //if (counts[ii].second/energies[ii].second >1E4) continue;
            if (counts[ii].first == 0) continue;
            bool skip = false;
            for(const auto& it_blacklist: blacklisted_energies){
                if (abs(energies[ii].first-it_blacklist)<0.1)
                    skip = true;
            }
            if (skip) continue;
            //if (counts[ii].second/intensities[ii]>1E4) continue;
            X.push_back(energies[ii].first);
            Y.push_back(counts[ii].first/energies[ii].second);
            Y3.push_back(counts[ii].first);
            Yerr.push_back(counts[ii].second/energies[ii].second);
            Y3err.push_back(counts[ii].second);
            Xerr.push_back(0);
            Y2.push_back(energies[ii].second);
        }
        relative_effgraph = TGraphErrors(X.size(), &X[0], &Y[0], &Xerr[0], &Yerr[0]);
        relative_integralgraph = TGraphErrors(X.size(), &X[0], &Y3[0], &Xerr[0], &Y3err[0]);
        relative_intgraph = TGraph(X.size(), &X[0], &Y2[0]);

        relative_effgraph.SetName("relative_eff_graoh");
        relative_effgraph.SetTitle("Relative eff graoh");
        relative_effgraph.SetMarkerColor(4);
        relative_effgraph.SetMarkerStyle(21);

        relative_integralgraph.SetName("relative_integral_graoh");
        relative_integralgraph.SetTitle("Relative integral graoh");
        relative_integralgraph.SetMarkerColor(2);
        relative_integralgraph.SetMarkerStyle(21);
        //if (canvas)
        //    plotter.WriteOnCanvas(std::string("energy : ")+
        //                        std::to_string(it_eff.first)+
        //                        " rel.eff : "+
        //                        std::to_string(counts)+
        //                        " nr : "+
        //                        std::to_string(it_transition.front()),
        //                      46);
    }
}


