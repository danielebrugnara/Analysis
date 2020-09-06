#include "SpectrumAnalyzer.h"

SpectrumAnalyzer::SpectrumAnalyzer(const std::string & file_name, const bool& debug_canvas):
        file_name(file_name),
        debug_canvas(debug_canvas),
        simulation(false),
        effgraph(),
        sigmagraph(),
        sigmagaussgraph(),
        relative_effgraph(),
        relative_integralgraph(),
        relative_intgraph(),
        levelscheme("files/adoptedLevels152Sm.csv"),
        fit_interval(15.),
        proj_interval(30.),
        fit_counter(0){

    eu152_intensities = ReadIntensities("files/intensity152Eu.txt");

    auto* file = new TFile(file_name.c_str());
    if (!file->IsOpen())
        throw std::runtime_error("File not found\n");

    auto gg_ptr = dynamic_cast<TH2D*>(file->Get("mgamma_gamma"));
    auto hspec_ptr = dynamic_cast<TH1D*>(file->Get("hspec"));

    if (gg_ptr == nullptr || hspec_ptr == nullptr) {
        simulation = true;
        gg_ptr = dynamic_cast<TH2D *>(file->Get("addb_gg"));
        hspec_ptr = dynamic_cast<TH1D *>(file->Get("addb_spec"));
        auto* tmp =  dynamic_cast<TVector2 *>(file->Get("start_stop"));
        nevts = tmp->Y();
    }

    if (gg_ptr == nullptr || hspec_ptr == nullptr)
        throw std::runtime_error("gamma_gamma or spec not found\n");

    if (!simulation)
        start_stop = *dynamic_cast<TVector2*>(file->Get("start_stop"));
    else
        start_stop = TVector2();


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
    if (!simulation){
        double acq_time = start_stop.Y()-start_stop.X();
        double source_activity = 22296.;
        double ndays = 2359.;
        double t12 = 13.517 *365.2425; //t12 in years * ndays
        double tau = t12 /log(2); //t12 in years * ndays
        source_activity = source_activity * exp(-ndays/tau);
        nevts = acq_time*source_activity;
    }
    scaled_relative_effgraph.SetName("scaled_relative_effgraph");
    scaled_relative_effgraph.SetTitle("scaled_relative_effgraph");
    scaled_relative_effgraph.SetMarkerColor(8);
    scaled_relative_effgraph.SetMarkerStyle(21);
    for (int ii=0; ii<relative_effgraph.GetN(); ++ii){
        scaled_relative_effgraph.SetPoint(ii,
                                          relative_effgraph.GetPointX(ii),
                                          relative_effgraph.GetPointY(ii)/nevts);
        scaled_relative_effgraph.SetPointError(     ii,
                                                    0,
                                                    relative_effgraph.GetErrorY(ii)/nevts);
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

    sigmagraph.Write();
    sigmagaussgraph.Write();
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
    bkgx.Sumw2();
    //UpdateErrors(bkgx);
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
    TH1D spect = hspec;
    UpdateErrors(spect);
    TH1D* subtrspec = Subtract(spect);



    std::vector<std::vector<int>> seen_transitions;
    double near_peak_threash = 6.5;
    for (unsigned long int ii=0; ii<eu152_intensities.size(); ++ii) {
        if(eu152_intensities[ii].second < 0.001) continue;
        std::vector<int> tmp_vec;
        tmp_vec.push_back(ii);
        while(ii<eu152_intensities.size() && abs(eu152_intensities[ii+1].first-eu152_intensities[ii].first)<near_peak_threash){
            tmp_vec.push_back(++ii);
        }
        seen_transitions.push_back(tmp_vec);
    }

    std::vector<double> X, Xerr, Y1, Y1err, Y2, Y3, Y3err, Y4, Y4err, Y5, Y5err;
    std::vector<double> blacklisted_energies;
//    std::vector<double> blacklisted_energies = {
//            329.425,330.54,503.474,764.9,768.944,926.317,930.58,937.05,1249.94,1457.64
//    };
    int fit_idx = 0;

    bool read_pars_from_file = false;

    std::string pars_file_name = "files/fitparams_";
    pars_file_name += file_name.substr(0,file_name.find(".root"));
    pars_file_name += ".txt";
    if (!read_pars_from_file) {
        std::ofstream pars_file(pars_file_name);
        pars_file << "List of parameters" << std::endl;
        pars_file.close();
    }

    for (const auto& it_transition: seen_transitions){
        std::vector<std::pair<double, double>> energies;
        energies.reserve(it_transition.size());
        for (const auto & it: it_transition){
            energies.push_back(eu152_intensities[it]);
        }
        //Fitter fitter(*subtrspec, energies);
        Fitter fitter(spect, energies);

        if (read_pars_from_file)
            fitter.ReadParsFromFile(pars_file_name,fit_idx++);
        else
            fitter.WriteParsOnFile(pars_file_name, fit_idx++);

        auto results = fitter.Fit();



        if (debug_canvas){
            //canvas->WaitPrimitive();
            //auto fitref = fitter.GetFitRef();
            //auto specref = fitter.GetSpecRef();

            //plotter.PlotOnCanvas<TH1D>(specref, "histo");
            //plotter.PlotOnCanvas<TF1>(fitref, "same");
            //plotter.SetRange(fitref.GetXmin(),fitref.GetXmax());
        }

        //counts /= it_eff.second;
        for (unsigned long int ii=0; ii<results.size(); ++ii){
            //if (counts[ii].first<2E3) continue;
            //if (counts[ii].second/energies[ii].second >1E4) continue;
            if (results[ii].integral.first == 0) continue;
            //if (results[ii]..first > 6) continue;
            bool skip = false;
            for(const auto& it_blacklist: blacklisted_energies){
                if (abs(energies[ii].first-it_blacklist)<0.1)
                    skip = true;
            }
            if (skip) continue;
            //if (counts[ii].second/intensities[ii]>1E4) continue;

            //
            X.push_back(energies[ii].first);
            Xerr.push_back(0);

            //
            Y1.push_back(results[ii].integral.first/energies[ii].second);
            Y1err.push_back(results[ii].integral.second/energies[ii].second);

            //
            Y2.push_back(energies[ii].second);

            //
            Y3.push_back(results[ii].integral.first);
            Y3err.push_back(results[ii].integral.second);

            //
            Y4.push_back(results[ii].ampl.first);
            Y4err.push_back(results[ii].ampl.second);

            //
            Y5.push_back(results[ii].sigma_gauss.first);
            Y5err.push_back(results[ii].sigma_gauss.second);
        }

        relative_effgraph = TGraphErrors(X.size(), &X[0], &Y1[0], &Xerr[0], &Y1err[0]);
        relative_effgraph.SetName("relative_eff_graph");
        relative_effgraph.SetTitle("Relative eff graph");
        relative_effgraph.SetMarkerColor(4);
        relative_effgraph.SetMarkerStyle(21);

        relative_intgraph = TGraph(X.size(), &X[0], &Y2[0]);
        relative_intgraph.SetName("relative_int_graph");
        relative_intgraph.SetTitle("Relative int graph");
        relative_intgraph.SetMarkerColor(3);
        relative_intgraph.SetMarkerStyle(21);

        relative_integralgraph = TGraphErrors(X.size(), &X[0], &Y3[0], &Xerr[0], &Y3err[0]);
        relative_integralgraph.SetName("relative_integral_graoh");
        relative_integralgraph.SetTitle("Relative integral graoh");
        relative_integralgraph.SetMarkerColor(2);
        relative_integralgraph.SetMarkerStyle(21);

        sigmagraph = TGraphErrors(X.size(), &X[0], &Y4[0], &Xerr[0], &Y4err[0]);
        sigmagraph.SetName("sigma_graph");
        sigmagraph.SetTitle("sigma_graph");
        sigmagraph.SetMarkerColor(5);
        sigmagraph.SetMarkerStyle(21);

        sigmagaussgraph = TGraphErrors(X.size(), &X[0], &Y5[0], &Xerr[0], &Y5err[0]);
        sigmagaussgraph.SetName("sigmagauss_graph");
        sigmagaussgraph.SetTitle("sigmagauss_graph");
        sigmagaussgraph.SetMarkerColor(6);
        sigmagaussgraph.SetMarkerStyle(21);
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


