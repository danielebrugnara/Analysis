#include "Fitter.h"

Fitter::Fitter(const TH1D& spec, std::vector<std::pair<double, double>> energies, const bool& left_tail):
        spec(spec),
        left_tail(left_tail),
        energies(std::move(energies)),
        canvas_enabled(false){
}

std::vector<Fitter::FitRes> Fitter::Fit(){

    if (parameters.sigma.empty())
        FindParameters();
    //return FindFunction();
    // Declare observable x
    RooRealVar x("x", "x", parameters.left_interval, parameters.right_interval);

    // Create a binned dataset that imports contents of TH1 and associates its contents to observable 'x'
    RooDataHist dh("dh", "dh", x, RooFit::Import(spec));

    // Make plot of binned dataset showing Poisson error bars (RooFit default)
    std::string title = "Fit at energies ";
    for(const auto& it: energies){
        title += " ";
        title += std::to_string(it.first);
    }
    RooPlot *frame = x.frame(RooFit::Title(title.c_str()));

    int ngauss=energies.size();
    //RooFitParams params[ngauss];
    std::vector<RooFitParams> params;
    //alternative to decay: RooAbsPdf* beta = bindPdf("beta",ROOT::Math::beta_pdf,x2,a,b) ;
    params.reserve(ngauss);
    for(int i=0; i<ngauss;++i){
        params.emplace_back(
                 new RooRealVar(parameters.ampl[i].name.c_str(),
                                                parameters.ampl[i].name.c_str(),
                                                parameters.ampl[i].initial_value,
                                                parameters.ampl[i].min,
                                                parameters.ampl[i].max),
                 new RooRealVar(parameters.mean[i].name.c_str(),
                                parameters.mean[i].name.c_str(),
                                parameters.mean[i].initial_value,
                                parameters.mean[i].min,
                                parameters.mean[i].max),
                 new RooRealVar(parameters.sigma[i].name.c_str(),
                                parameters.sigma[i].name.c_str(),
                                parameters.sigma[i].initial_value,
                                parameters.sigma[i].min,
                                parameters.sigma[i].max),
                 new RooRealVar(parameters.tau[i].name.c_str(),
                                parameters.tau[i].name.c_str(),
                                parameters.tau[i].initial_value,
                                parameters.tau[i].min,
                                parameters.tau[i].max)
        );
        params.back().fracsignal->setConstant(parameters.ampl[i].fixed);
        params.back().mean->setConstant(parameters.mean[i].fixed);
        //params.back().mean->setConstant(true);
        params.back().sigma->setConstant(parameters.sigma[i].fixed);
        params.back().tau->setConstant(parameters.tau[i].fixed);
    }
    params[ngauss-1].fracsignal->setConstant(false);

    std::vector<RooGaussModel*> gaussians;
    std::vector<RooDecay*> smeared_gaussians;
    RooArgList smeared_gaussians_list;
    RooArgList simple_gaussians_list;
    RooArgList signals_list;
//    std::vector<RooAbsPdf*> partial_model;
    for (int i=0; i<ngauss; ++i){
        gaussians.push_back(new RooGaussModel(Form("gaussian_%i",i),Form("gaussian_%i",i), x, *params[i].mean, *params[i].sigma));
        smeared_gaussians.push_back(new RooDecay(Form("smeared_gaussian_%i", i),Form("smeared_gaussian_%i", i),x,*params[i].tau,*gaussians[i],RooDecay::Flipped));
        smeared_gaussians_list.add(*smeared_gaussians.back());
        simple_gaussians_list.add(*gaussians.back());
        //signals_list.add(*params[i].fracsignal);
    }
    RooRealVar pol1("pol1", "pol1", 0, -0.01, 0.01);
    //RooRealVar pol2("pol2", "pol2", 0, -0.01, 0.01);
    RooPolynomial bkg("bkg", "bkg", x, RooArgList(pol1));
    auto* sign_y = new RooRealVar("sign_y","sign_y", 1E3, 0.,1E7);
    auto* bkg_y = new RooRealVar("bkg_y","bkg_y", 1E3, 0.,1E7);

    for (int i=0; i<ngauss-1; ++i) {
        signals_list.add(*params[i].fracsignal);
    }
    params[ngauss-1].fracsignal = sign_y;

    RooAddPdf* signal = nullptr;
    if (left_tail) {
        signal = new RooAddPdf(Form("signal"),
                            Form("signal"),
                                smeared_gaussians_list,
                                signals_list);
    }else{
        signal = new RooAddPdf(Form("signal"),
                                Form("signal"),
                                simple_gaussians_list,
                                signals_list);
    }

    auto* model = new RooAddPdf(Form("bkg_gaus%i",ngauss-1),
                                Form("bkg_gaus%i", ngauss-1),
                                RooArgList(*signal,bkg),
                                RooArgList(*sign_y,*bkg_y));

    dh.plotOn(frame);
    auto *result = model->fitTo(dh, RooFit::Save());
    model->plotOn(frame, RooFit::VisualizeError(*result));
    //model->plotOn(frame, RooFit::VisualizeError(*result, *sign_y));
    model->plotOn(frame );
    model->paramOn(frame, RooFit::Layout(0.7, 0.9, 0.9));
    dh.plotOn(frame);
    std::vector<TF1* > root_smeared_gaussians;
    for (int i=0; i<ngauss; ++i) {
        //smeared_gaussians[i]->plotOn(frame, RooFit::LineStyle(kDashed), RooFit::LineColor(kRed));
        model->plotOn(frame, RooFit::Components(*smeared_gaussians[i]), RooFit::LineStyle(kDashed),RooFit::LineColor(kGreen));
        root_smeared_gaussians.push_back(smeared_gaussians[i]->asTF(RooArgList(x),
                                                                    *smeared_gaussians[i]->getParameters(RooArgList(x))));
    }
//    model->plotOn(frame, RooFit::Components(bkg), RooFit::LineStyle(kDashed),RooFit::LineColor(kRed));
    parameters.GetParameters(params);


    //result->Print();

    //RooAbsReal
//    RooAbsMoment *aa = smeared_gaussians[0]->sigma(x);
//    std::cout << "sigma : " << aa->sigma(x)->getVal(x);
//    std::cout << " pm. " << aa->mean()->getVal(x) << std::endl;

//    RooAbsMoment * mean= gaussians[0]->moment(RooArgList(x,*params[0].mean, *params[0].sigma),1, false, true);
//    std::cout << "sigma_gauss : " << mean->mean();
//    std::cout << " pm. " << mean->mean()->sigma(x) << std::endl;
    //smeared_gaussians[0]->sigma



    std::vector<Fitter::FitRes> results;
    x.setRange("increased_range", parameters.left_interval-10, parameters.right_interval+10);
    //RooAbsReal* integral_res = model->createIntegral(RooArgSet(x),RooArgSet(x),"x");
    double int_signal = sign_y->getVal()/spec.GetBinWidth(0);
    double int_signal_err = sign_y->getError()/spec.GetBinWidth(0);
    double int_bkg = bkg_y->getVal()/spec.GetBinWidth(0);
    //double integral_err = integral->getPropagatedError(*result,RooArgSet(x));
    //RooAbsReal* model_counts = model->createIntegral(x,RooFit::NormSet(x), RooFit::Range("x"));
    for(int ii=0; ii<ngauss; ++ii) {
        //double integral =int_signal* root_smeared_gaussians[ii]->Integral(parameters.left_interval,
        //                                                                  parameters.right_interval,
        //                                                                  0);
        //double integral_err =int_signal* root_smeared_gaussians[ii]->IntegralError(parameters.left_interval,
        //                                                                           parameters.right_interval,
        //                                                                           0,
        //                                                                           result->covarianceMatrix().GetMatrixArray());

        double ampl = parameters.ampl[ii].fitted_value.first;
        double ampl_err = parameters.ampl[ii].fitted_value.second;
        if (ii == ngauss-1){
            double remaining_ampl = 1.;
            double remaining_ampl_sqerr = 0;
            for(int jj=0; jj<ngauss-1; ++jj) {
                remaining_ampl -= parameters.ampl[jj].fitted_value.first;
                remaining_ampl_sqerr += pow(parameters.ampl[jj].fitted_value.second,2);
            }
            ampl = remaining_ampl;
            ampl_err = sqrt(remaining_ampl_sqerr);
        }
        double integral = int_signal * ampl;
        double integral_err = sqrt(pow(int_signal * ampl_err,2)+pow(int_signal_err*ampl,2)); //TODO: APPTOXIMATED!! https://root-forum.cern.ch/t/getting-number-of-events-or-fraction-of-events-in-a-range/28624 https://root-forum.cern.ch/t/integral-uncertainty-with-roofit/12622/3
        double sigma_gauss = parameters.sigma[ii].fitted_value.first;
        double sigma_gauss_err = parameters.sigma[ii].fitted_value.second;
        double tail = parameters.tau[ii].fitted_value.first;
        double tail_err = parameters.tau[ii].fitted_value.second;
//        double integral_err = sqrt( pow(parameters.ampl[ii].fitted_value.first * parameters.common_ampl.fitted_value.second ,2)+
//                                    pow(parameters.ampl[ii].fitted_value.second * parameters.common_ampl.fitted_value.first ,2));
//
//        double sigma = sqrt(pow(parameters.sigma[ii].fitted_value.first,2)+pow(parameters.tau[ii].fitted_value.first,2));
//        double sigma_err = 1./sigma*
//                           sqrt(    pow(parameters.sigma[ii].fitted_value.first*parameters.sigma[ii].fitted_value.second,2)+
//                                    pow(parameters.tau[ii].fitted_value.first*parameters.tau[ii].fitted_value.second,2)+
//                                    2*parameters.sigma[ii].fitted_value.first*parameters.tau[ii].fitted_value.first*cov[3+npars*ii][2+npars*ii]);
//
//        double sigma_gauss = parameters.sigma[ii].fitted_value.first;
//        double sigma_gauss_err = parameters.sigma[ii].fitted_value.second;
//
//        double tail = parameters.tau[ii].fitted_value.first;
//        double tail_err = parameters.tau[ii].fitted_value.second;
//
        results.emplace_back(
                energies[ii].first,
                integral, integral_err,
                ampl, ampl_err,
                sigma_gauss, sigma_gauss_err,
                tail, tail_err);
    }
    if (canvas_enabled) {
        auto *cv = new TCanvas();
        frame->Draw();
        cv->WaitPrimitive();
        delete cv;
    }
    return results;

}

void Fitter::FindParameters() {
    std::vector<double>fractions;
    double total_int {0};
    for(const auto&it: energies){
        total_int+=it.second;
    }

    for(const auto&it: energies) {
        fractions.push_back(it.second/total_int);
    }

    int number_of_gaussians = energies.size();
    InitializationParameters initialpars(number_of_gaussians);
    initialpars.min_sigma = 0.5;
    initialpars.max_sigma = 0.02*energies.back().first;
    initialpars.min_tau = 0.001;
    initialpars.max_tau = 10.;
    initialpars.left_interval = energies.front().first - 15.;
    initialpars.right_interval = energies.back().first + 15.;

    for (int ii = 0; ii < number_of_gaussians; ++ii) {
        initialpars.ampl[ii].initial_value = fractions[ii];
        //initialpars.ampl[ii].initial_value = 0.5;
        initialpars.ampl[ii].min = 0.;
        initialpars.ampl[ii].max = 1.;
        initialpars.ampl[ii].fixed = true;
        initialpars.ampl[ii].name = "Ampl_"+std::to_string(ii);

        initialpars.mean[ii].initial_value = energies[ii].first;
        initialpars.mean[ii].min = energies[ii].first - 5.;
        initialpars.mean[ii].max = energies[ii].first + 5.;
        initialpars.mean[ii].fixed = false;
        initialpars.mean[ii].name = "Mean_"+std::to_string(ii);

        initialpars.sigma[ii].initial_value = 1.;
        initialpars.sigma[ii].min = initialpars.min_sigma;
        initialpars.sigma[ii].max = initialpars.max_sigma;
        initialpars.sigma[ii].fixed = false;
        initialpars.sigma[ii].name = "Sigma_"+std::to_string(ii);

        initialpars.tau[ii].initial_value = 1.3;
        initialpars.tau[ii].min = initialpars.min_tau;
        initialpars.tau[ii].max = initialpars.max_tau;
        initialpars.tau[ii].fixed = false;
        initialpars.tau[ii].name = "Tau_"+std::to_string(ii);
    }

    initialpars.ampl[number_of_gaussians-1].initial_value = 1E4;
    initialpars.ampl[number_of_gaussians-1].min = 0.;
    initialpars.ampl[number_of_gaussians-1].max = 1E8;
    initialpars.ampl[number_of_gaussians-1].fixed = false;
    initialpars.ampl[number_of_gaussians-1].name = "Ampl_"+std::to_string(number_of_gaussians-1);

    parameters = initialpars;
}



Fitter::~Fitter() = default;

void Fitter::WriteParsOnFile(const std::string& file_name, const int& idx) {
    FindParameters();
    std::ofstream parfile(file_name, std::ios::app);

    if (!parfile.is_open())
        throw std::runtime_error("could not write file\n");

    parfile << "---------------------------"<< std::endl;
    parfile << "idx: " <<  idx << std::endl;
    parfile <<  " name "  << ":\t\t\t"
            <<  " init_val " << "\t\t\t"
            <<  " min " << "\t\t\t"
            <<  " max " << "\t\t\t"
            <<  " fixed "<< std::endl;
    for (unsigned long int ii=0; ii<parameters.ampl.size(); ++ii){
        parfile << std::scientific << parameters.ampl[ii].name << "\t\t\t"
                << std::scientific << parameters.ampl[ii].initial_value << "\t\t"
                << std::scientific << parameters.ampl[ii].min << "\t\t"
                << std::scientific << parameters.ampl[ii].max << "\t\t"
                << std::scientific << parameters.ampl[ii].fixed << std::endl;

        parfile << std::scientific << parameters.mean[ii].name << "\t\t\t"
                << std::scientific << parameters.mean[ii].initial_value << "\t\t"
                << std::scientific << parameters.mean[ii].min << "\t\t"
                << std::scientific << parameters.mean[ii].max << "\t\t"
                << std::scientific << parameters.mean[ii].fixed << std::endl;

        parfile << std::scientific << parameters.sigma[ii].name << "\t\t"
                << std::scientific << parameters.sigma[ii].initial_value << "\t\t"
                << std::scientific << parameters.sigma[ii].min << "\t\t"
                << std::scientific << parameters.sigma[ii].max << "\t\t"
                << std::scientific << parameters.sigma[ii].fixed << std::endl;

        parfile << std::scientific << parameters.tau[ii].name << "\t\t\t"
                << std::scientific << parameters.tau[ii].initial_value << "\t\t"
                << std::scientific << parameters.tau[ii].min << "\t\t"
                << std::scientific << parameters.tau[ii].max << "\t\t"
                << std::scientific << parameters.tau[ii].fixed << std::endl;
    }


    parfile << std::scientific << parameters.common_ampl.name << "\t\t"
            << std::scientific << parameters.common_ampl.initial_value << "\t\t"
            << std::scientific << parameters.common_ampl.min << "\t\t"
            << std::scientific << parameters.common_ampl.max << "\t\t"
            << std::scientific << parameters.common_ampl.fixed << std::endl;

    parfile << std::scientific << parameters.common_offset.name << "\t"
            << std::scientific << parameters.common_offset.initial_value << "\t\t"
            << std::scientific << parameters.common_offset.min << "\t\t"
            << std::scientific << parameters.common_offset.max << "\t\t"
            << std::scientific << parameters.common_offset.fixed << std::endl;

    parfile << std::scientific << "min_max" << "\t\t"
            << std::scientific << parameters.left_interval << "\t\t"
            << std::scientific << parameters.right_interval << std::endl;

    parfile.close();
}

void Fitter::ReadParsFromFile(const std::string & parfile_name, const int &idx) {
    std::ifstream parfile(parfile_name);
    if (!parfile.is_open())
        throw std::runtime_error("Parfile: "+parfile_name+" not opened\n");

    InitializationParameters tmp_pars;
    std::string line;
    std::vector<InitialState> initial_states;
    bool found = false;
    while(std::getline(parfile,line)){
        if (line.find("idx:") == std::string::npos && !found)
            continue;
        if (!found && std::stoi(line.substr(line.find("idx:")+4)) == idx){
            found = true;
            continue;
        }else if (!found)
            continue;
        if (line.find("-----") == 0 )
            break;

        if (line.find("name") != std::string::npos )
            continue;

        std::istringstream linestream(line);
        std::string name;
        double init_val;
        double min;
        double max;
        bool fixed;

        linestream >> name >> init_val >> min >> max >> fixed;
        initial_states.emplace_back(init_val,min,max,fixed);
        initial_states.back().name = name;
    }

    if ((initial_states.size()-3)%4 != 0)
        throw std::runtime_error("Error reading initialization parameters file, size :"+std::to_string(initial_states.size())+"\n");
    InitializationParameters tmp_parameters(((int)initial_states.size()-3)/4);
    for (const auto& it: initial_states){
        if(it.name.find("Ampl_") != std::string::npos){
            int nr = std::stoi(it.name.substr(it.name.find("Ampl_")+5));
            tmp_parameters.ampl[nr] = it;
        }
        if(it.name.find("Mean_") != std::string::npos){
            int nr = std::stoi(it.name.substr(it.name.find("Mean_")+5));
            tmp_parameters.mean[nr] = it;
        }
        if(it.name.find("Sigma_") != std::string::npos){
            int nr = std::stoi(it.name.substr(it.name.find("Sigma_")+6));
            tmp_parameters.sigma[nr] = it;
        }
        if(it.name.find("Tau_") != std::string::npos){
            int nr = std::stoi(it.name.substr(it.name.find("Tau_")+4));
            tmp_parameters.tau[nr] = it;
        }
        if(it.name.find("CommonAmpl") != std::string::npos){
            tmp_parameters.common_ampl = it;
        }
        if(it.name.find("CommonOffset") != std::string::npos){
            tmp_parameters.common_offset = it;
        }
        if(it.name.find("min_max") != std::string::npos){
            tmp_parameters.left_interval = it.initial_value;
            tmp_parameters.right_interval = it.min;
        }

    }
    parameters = tmp_parameters;
    if (tmp_parameters.sigma.empty()){
        std::cerr << "Could not find parameters, using default ones\n";
        FindParameters();
    }
}

void Fitter::EnableCanvas(bool canvas_enabled) {
    this->canvas_enabled = canvas_enabled;
}
