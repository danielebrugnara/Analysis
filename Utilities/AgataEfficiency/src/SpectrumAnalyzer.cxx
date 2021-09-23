#include "SpectrumAnalyzer.h"

SpectrumAnalyzer::SpectrumAnalyzer(const std::string & file_name, const bool& debug_canvas, const std::string& use_spectra, const bool& use_main_transitions):
        file_name(file_name),
        debug_canvas(debug_canvas),
        use_main_transitions(use_main_transitions),
        simulation(false),
        levelscheme("files/adoptedLevels152Sm_withalpha.csv"),
        fit_interval(15.),
        proj_interval(30.),
        fit_counter(0){

    eu152_intensities = ReadIntensities("files/intensity152Eu.txt");

    auto* file = new TFile(file_name.c_str());
    if (!file->IsOpen())
        throw std::runtime_error("File not found\n");

    TH2D* gg_ptr = nullptr;
    TH1D* hspec_ptr = nullptr;

    std::string spect_name;
    std::string gg_name;
    if (use_spectra=="addb") {
        gg_name = "data_addb_gg";
        spect_name = "data_addb_spec";
    }else if(use_spectra=="core"){
        gg_name = "data_core_gg";
        spect_name = "data_core_spec";
    }else if(use_spectra=="tracked"){
        gg_name = "data_tracked_gg";
        spect_name = "data_tracked_spec";
    }else{
        throw std::runtime_error("unknown specrtrum choice\n");
    }

    gg_ptr = dynamic_cast<TH2D*>(file->Get(gg_name.c_str()));
    hspec_ptr = dynamic_cast<TH1D*>(file->Get(spect_name.c_str()));
    auto* start_stop_ptr= dynamic_cast<TVector2*>(file->Get("start_stop"));
    if (start_stop_ptr == nullptr)
        throw std::runtime_error("No start_stop vector, do you need to run the selector?\n");
    start_stop = *start_stop_ptr;
    if (gg_ptr == nullptr || hspec_ptr == nullptr) {//probabily a simulation
        simulation = true;
        gg_name = gg_name.substr(5);
        spect_name = spect_name.substr(5);
        gg_ptr = dynamic_cast<TH2D *>(file->Get(gg_name.c_str()));
        hspec_ptr = dynamic_cast<TH1D *>(file->Get(spect_name.c_str()));
        auto* tmp =  dynamic_cast<TVector2 *>(file->Get("start_stop"));
        if (tmp != nullptr) {
            start_stop = *tmp;
            nevts = tmp->Y();
        }
    }

    if (gg_ptr == nullptr)
        throw std::runtime_error("gamma_gamma not found, do you need to run the selector?\n");
    if (hspec_ptr == nullptr)
        throw std::runtime_error("spectrum not found, do you need to run the selector?\n");

    gg = *gg_ptr;
    hspec = *hspec_ptr;

    for (int i=0; i<46; ++i){
        cry_spec.push_back(*dynamic_cast<TH1D*>(file->Get(Form("crystal_spectra_%i",i))));
    }

    file->Close();

    gamma_gamma = levelscheme.GetGammaGamma();
//    for (const auto & it: gamma_gamma){
//        std::cout <<"~~" << *it.first <<" -> "<< *it.second<< std::endl;
//    }
    Analyze();
}

SpectrumAnalyzer::~SpectrumAnalyzer() = default;

void SpectrumAnalyzer::Analyze() {

    GenerateRelativeEffGraph(hspec, relative_eff_graph, sigma_graph, tau_graph , 100);
    GenerateAbsoluteEffGraph();
    //GenerateCrystalEffGraph();
    std::cout << "nevts " << nevts << std::endl;

    std::string out_file_name = "eff_curve_";
    out_file_name +=file_name;

    TFile outfile(out_file_name.c_str(),"recreate");
    relative_eff_graph.Write();
    int ii = 0;
    for(auto& it: cry_relative_eff_graph) {
        std::ofstream outfile;
        outfile.open("test.txt", std::ios_base::app); // append instead of overwrite
        outfile << ii++ << "\n";
        outfile.close();
        it->Write();
    }

    FitEffCurve(relative_eff_graph).Write();

    if (debug_canvas) {
        TCanvas cv;
        relative_eff_graph.Draw();
        FitEffCurve(relative_eff_graph).Draw("same");
        cv.WaitPrimitive();
    }

    absolute_eff_graph.Write();
    sigma_graph.Write();
    tau_graph.Write();

    outfile.Write();
    outfile.Close();
}

void SpectrumAnalyzer::GenerateAbsoluteEffGraph() {
    int nponts=0;
    double Br_threashold = 0.01;
    double Integral_threashold1 = 1000;
    double Integral_threashold2 = 10;
    int Entries_threashold = 1E3;
    TH1D proj_x  = *gg.ProjectionX();
    int cnt{0};
    std::vector<int> meaningful_combinations = {
            1, 5, 24, 31, 35, 46, 47, 99, 158, 167, 168,
            169, 170, 197, 298, 302, 479, 535, 561,
            684, 755, 782, 885, 1240, 1359, 1397
    };
    int fit_x_idx{0};
    int fit_y_idx{0};
    std::vector<double> X, X_err, Y_eff, Y_eff_err;

    bool read_pars_from_file = true;

    std::string pars_file_x_name = "files/fitparams_absolutex_";
    pars_file_x_name += file_name.substr(0, file_name.find(".root"));
    pars_file_x_name += ".txt";
    std::string pars_file_y_name = "files/fitparams_absolutey_";
    pars_file_y_name += file_name.substr(0, file_name.find(".root"));
    pars_file_y_name += ".txt";

    if (read_pars_from_file){
        std::ifstream test_file_x(pars_file_x_name);
        if (!test_file_x.good())
            read_pars_from_file = false;
        else
            test_file_x.close();

        std::ifstream test_file_y(pars_file_y_name);
        if (!test_file_y.good())
            read_pars_from_file = false;
        else
            test_file_y.close();
    }

    if (!read_pars_from_file){
        std::ofstream pars_x_file(pars_file_x_name);
        pars_x_file << "List of parameters X" << std::endl;
        pars_x_file.close();
        std::ofstream pars_y_file(pars_file_y_name);
        pars_y_file << "List of parameters Y" << std::endl;
        pars_y_file.close();
    }

    for (const auto & it: gamma_gamma) {
        cnt++;
        if (it.first->Br.first < Br_threashold || it.second->Br.first < Br_threashold)
            continue;
        if (std::find(meaningful_combinations.begin(), meaningful_combinations.end(), cnt) ==
            meaningful_combinations.end())
            continue;

        std::cout << "\n\n\n------------------------gammma-gamma--->" << *it.first << " -> " << *it.second << "\n\n\n";
        double engamma1 = it.first->Egamma.first / UNITS::keV;
        double engamma2 = it.second->Egamma.first / UNITS::keV;
        double branching2 = it.second->Br.first;
        double alpha = it.second->Alpha.first;

        std::vector<std::pair<double, double>> energies_x = {std::make_pair(engamma1, 1.)};
        std::vector<std::pair<double, double>> energies_y = {std::make_pair(engamma2, 1.)};

        //Fitting X projection
        bool left_tails = !simulation;
        Fitter fitter_x(proj_x, energies_x, left_tails);
        if (debug_canvas)
            fitter_x.EnableCanvas();

        if (read_pars_from_file)
            fitter_x.ReadParsFromFile(pars_file_x_name, cnt);
        else {
            fitter_x.WriteParsOnFile(pars_file_x_name, cnt);
        }

        auto results_x = fitter_x.Fit();
        double  integral_x = results_x[0].integral.first;
        double  integral_x_err = results_x[0].integral.second;

        //Creating Y projection
        double gate_width = 5.;
        double noise_gate_width = 50.;
        double noise_gate_dist = 10.;

        bool subtract_bkg = true;
        TH1D* proj_y = nullptr;

        if (subtract_bkg) {
            proj_y = ProjectAndSubtract(gg,
                                        engamma1 - gate_width / 2.,
                                        engamma1 + gate_width / 2.,
                                        engamma1 - noise_gate_dist - noise_gate_width,
                                        engamma1 - noise_gate_dist,
                                        engamma1 + noise_gate_dist,
                                        engamma1 + noise_gate_dist + noise_gate_width);
        }else{
            proj_y = new TH1D(*(gg.ProjectionY("proj_nog",
                                gg.GetXaxis()->FindBin(engamma1-gate_width/2.),
                                gg.GetXaxis()->FindBin(engamma1+gate_width/2.))));
        }

        //Fitting Y projection
        Fitter fitter_y(*proj_y, energies_y, left_tails);
        if (debug_canvas)
            fitter_y.EnableCanvas();

        if (read_pars_from_file)
            fitter_y.ReadParsFromFile(pars_file_y_name, cnt);
        else {
            fitter_y.WriteParsOnFile(pars_file_y_name, cnt);
        }

        auto results_y = fitter_y.Fit();
        double  integral_y = results_y[0].integral.first;
        double  integral_y_err = results_y[0].integral.second;

        //Computing efficiency
        double eff = integral_y/(integral_x*branching2)*(1+alpha);
        double eff_err = 1./branching2* sqrt(pow(integral_y_err/(integral_x),2.)+
                                                pow(integral_y*integral_x_err/(integral_x*integral_x),2.))*(1+alpha);

        if (eff>1.)
            continue;

        X.push_back(engamma2);
        X_err.push_back(0.);

        Y_eff.push_back(eff);
        Y_eff_err.push_back(eff_err);

        //std::cout << " ener : " << X.back() << " eff : " << Y_eff.back() << std::endl;
        //TCanvas cv;
        //cv.WaitPrimitive();



        delete proj_y;
    }
    absolute_eff_graph = TGraphErrors(X.size(), &X[0], &Y_eff[0], &X_err[0], &Y_eff_err[0]);
    absolute_eff_graph.SetName("absolute_eff_graph");
    absolute_eff_graph.SetTitle("Absolute eff graph");
    absolute_eff_graph.SetMarkerColor(7);
    absolute_eff_graph.SetMarkerStyle(21);
    absolute_eff_graph.SetDrawOption("AP");
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

void SpectrumAnalyzer::GenerateRelativeEffGraph(const TH1D& spectrum, TGraphErrors& relative_eff_graph, TGraphErrors& sigma_graph, TGraphErrors& tau_graph, int nbr) {
    //TH1D projx  = *gg.ProjectionX();
    //TH1D spect = hspec;
    TH1D spect = spectrum;
    spect.Sumw2();

    std::vector<std::vector<int>> seen_transitions;
    double near_peak_threash = 6.5;
    for (unsigned long int ii=0; ii<eu152_intensities.size(); ++ii) {
        if(eu152_intensities[ii].second < 0.001) continue;
        std::vector<int> tmp_vec;
        tmp_vec.push_back(ii);
        while(ii<eu152_intensities.size() &&
                abs(eu152_intensities[ii+1].first-eu152_intensities[ii].first)<near_peak_threash){
            tmp_vec.push_back(++ii);
        }
        seen_transitions.push_back(tmp_vec);
    }

    std::vector<double> X, X_err, Y_eff, Y_eff_err, Y_sigma, Y_sigma_err, Y_tau, Y_tau_err;
//    std::vector<double> blacklisted_energies;
    std::vector<double> blacklisted_energies = {
          295.939, 503.474,586.2648,656.487,719.349,1249.938, 1457.643
    };
    int fit_idx = 0;

    bool read_pars_from_file = true;
    std::string pars_file_name = "files/fitparams_relative_";
    pars_file_name += file_name.substr(0,file_name.find(".root"));
    pars_file_name += ".txt";

    if (read_pars_from_file){
        std::ifstream test_file(pars_file_name);
        if (!test_file.good())
            read_pars_from_file = false;
        else
            test_file.close();
    }

    if (!read_pars_from_file) {
        std::ofstream pars_file(pars_file_name);
        pars_file << "List of parameters" << std::endl;
        pars_file.close();
    }


    for (const auto& it_transition: seen_transitions){
        std::vector<std::pair<double, double>> energies;
        energies.reserve(it_transition.size());
        double max_intensity = 0;
        for (const auto & it: it_transition){
            energies.push_back(eu152_intensities[it]);
            if (energies.back().second > max_intensity)
                max_intensity = energies.back().second;
        }

        bool left_tails = !simulation;
        Fitter fitter(spect, energies, left_tails);
        if (debug_canvas)
            fitter.EnableCanvas();

        if (read_pars_from_file)
            fitter.ReadParsFromFile(pars_file_name,fit_idx++);
        else
            fitter.WriteParsOnFile(pars_file_name, fit_idx++);

        if (max_intensity < 0.01 && use_main_transitions) continue; //skips low intansity transitions
        auto results = fitter.Fit();

        for (unsigned long int ii=0; ii<results.size(); ++ii){
            if (results[ii].integral.first == 0) continue;
            bool skip = false;
            for(const auto& it_blacklist: blacklisted_energies){
                if (abs(energies[ii].first-it_blacklist)<0.1)
                    skip = true;
            }
            double rel_max_err = 0.07;
            double eff = results[ii].integral.first/energies[ii].second;
            double eff_error = results[ii].integral.second/energies[ii].second;

            //if (eff>1.) continue;
            if (abs(eff_error/eff) > rel_max_err) continue;
            if (skip) continue;

            X.push_back(energies[ii].first);
            X_err.push_back(0);

            Y_eff.push_back(eff);
            Y_eff_err.push_back(eff_error);

            Y_sigma.push_back(results[ii].sigma_gauss.first);
            Y_sigma_err.push_back(results[ii].sigma_gauss.second);

            Y_tau.push_back(results[ii].tail.first);
            Y_tau_err.push_back(results[ii].tail.second);

        }

        std::cout << "Writing graphs!!!\n";
        relative_eff_graph = TGraphErrors(X.size(), &X[0], &Y_eff[0], &X_err[0], &Y_eff_err[0]);
        relative_eff_graph.SetName("relative_eff_graph");
        relative_eff_graph.SetTitle("Relative eff graph");
        relative_eff_graph.SetMarkerColor(4);
        relative_eff_graph.SetMarkerStyle(21);

        sigma_graph = TGraphErrors(X.size(), &X[0], &Y_sigma[0], &X_err[0], &Y_sigma_err[0]);
        sigma_graph.SetName("sigma_graph");
        sigma_graph.SetTitle("sigma_graph");
        sigma_graph.SetMarkerColor(5);
        sigma_graph.SetMarkerStyle(21);

        tau_graph = TGraphErrors(X.size(), &X[0], &Y_tau[0], &X_err[0], &Y_tau_err[0]);
        tau_graph.SetName("tau_graph");
        tau_graph.SetTitle("tau_graph");
        tau_graph.SetMarkerColor(6);
        tau_graph.SetMarkerStyle(21);
    }

    if (!simulation){
        double acq_time = start_stop.Y()-start_stop.X();
        //double source_activity = 22296.;
        //double ndays = 2359.;
        double source_activity = 15524.;
        double ndays = -222.;
        double t12 = 13.517 *365.2425; //t12 in years * ndays
        double tau = t12 /log(2); //t12 in years * ndays
        source_activity = source_activity * exp(-ndays/tau);
        nevts = acq_time*source_activity;
    }
    for (int ii=0; ii<relative_eff_graph.GetN(); ++ii){
        relative_eff_graph.SetPoint(ii,
                                   relative_eff_graph.GetPointX(ii),
                                   relative_eff_graph.GetPointY(ii)/nevts);
        relative_eff_graph.SetPointError(ii,
                                        0,
                                        relative_eff_graph.GetErrorY(ii)/nevts);
    }
    std::cout << "number of events : " << nevts << std::endl;
    relative_eff_graph.SaveAs(Form("data_%i.root", nbr));
}

void SpectrumAnalyzer::GenerateCrystalEffGraph() {
    for (int i=0; i<cry_spec.size(); ++i){
        //if (it.GetEntries() == 0 ) continue;

        TGraphErrors rel;
        TGraphErrors sig;
        TGraphErrors tau;
        GenerateRelativeEffGraph(cry_spec[i], rel, sig, tau, i);
        cry_relative_eff_graph.push_back(&rel);
    }
}

TF1 SpectrumAnalyzer::FitEffCurve(TGraphErrors& gr) {
    TF1 effcurve("effcurve",
                 [](const double*x, const double*par){
                     return (par[0]+
                             par[1]*log(x[0])+
                             par[2]*pow(log(x[0]), 2)+
                             par[3]*pow(log(x[0]), 3)+
                             par[4]*pow(log(x[0]), 5)+
                             par[5]*pow(log(x[0]), 7))/x[0];
                 },
                 50,
                 1500,
                 6);
    gr.Fit(&effcurve, "M");
    return effcurve;
}