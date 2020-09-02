#include "Fitter.h"

Fitter::Fitter(TH1D spec, std::vector<std::pair<double, double>> energies):
        spec(std::move(spec)),
        energies(std::move(energies)){
}

std::vector<Fitter::FitRes> Fitter::Fit(){
    FindParameters();
    return FindFunction();
}

std::vector<Fitter::FitRes> Fitter::Fit(const std::string& file, const int& idx) {
    parameters = ReadParsFromFile(file, idx);
    return FindFunction();
}

void Fitter::FindParameters() {
    std::vector<double> heights;
    std::vector<double> sigmas;
    for (const auto &it_en: energies) {
        std::cout << "===== Single peak fit :" << it_en.first << "\n";
        double single_fit_interval{3.0};
        double single_min_sigma{0.5};
        double single_max_sigma{4.0};

        TF1 fitfunc(Form("gaussian_%f", it_en.first),
                    [](const double *x, const double *par) {
                        return par[0] * exp(-pow((x[0] - par[1]) / (par[2]), 2));
                    },
                    it_en.first - single_fit_interval,
                    it_en.first + single_fit_interval,
                    3,
                    1);
        fitfunc.SetNpx(400);

        fitfunc.SetParameter(0, std::max(0., spec.GetBinContent(spec.GetXaxis()->FindBin(it_en.first))));
        fitfunc.SetParLimits(0, std::max(0., fitfunc.GetParameter(0) / 4.),
                             std::max(fitfunc.GetParameter(0) * 3., 40.));

        fitfunc.SetParameter(1, it_en.first);
        fitfunc.SetParLimits(1, it_en.first - 10, it_en.first + 10);

        fitfunc.SetParameter(2, 1.8);
        fitfunc.SetParLimits(2, single_min_sigma, single_max_sigma);

        fitfunc.FixParameter(1, it_en.first);

        TFitResultPtr fit_res_ptr(0);
        fit_res_ptr = spec.Fit(&fitfunc,
                               "S",
                               "",
                               fitfunc.GetXmin(),
                               fitfunc.GetXmax());

        heights.push_back(fitfunc.GetParameter(0));
        TFitResult fit_result = *fit_res_ptr;

        if (fitfunc.GetParameter(2) < sqrt(single_max_sigma) &&
            fitfunc.GetParameter(2) > sqrt(single_min_sigma) &&
            fit_result.IsValid())
            sigmas.push_back(fitfunc.GetParameter(2));
        else
            sigmas.push_back((sqrt(single_max_sigma) + sqrt(single_min_sigma)) * 0.5);

    }

    //Multifit
    //https://en.wikipedia.org/wiki/Exponentially_modified_Gaussian_distribution
    // par[0+ii*npars] -> coeff (integral of distribution, since the rest is normalized)
    // par[1+ii*npars] -> gaussian mean value
    // par[2+ii*npars] -> gaussian sigma
    // par[3+ii*npars] -> exponential tau  (exp(-x/tau))

    //integral: par[0+ii*npars] -> actual integral is par[0+ii*npars]*par[number_of_gaussians*npars]
    //mean: par[1+ii*npars]-par[3+ii*npars]
    //sigma: sqrt(par[1+ii*npars]^2 + par[3+ii*npars]^2)

    //double lim_inf = energies.front().first - fit_interval;
    //double lim_sup = energies.back().first + fit_interval;

    int number_of_gaussians = energies.size();
    InitializationParameters initialpars(number_of_gaussians);
    initialpars.min_sigma = 0.5;
    initialpars.max_sigma = 4.0;
    initialpars.min_tau = 0.06;
    initialpars.max_tau = 8.;
    initialpars.left_interval = energies.front().first - 14.;
    initialpars.right_interval = energies.back().first + 10.;

    for (int ii = 0; ii < number_of_gaussians; ++ii) {
        initialpars.ampl[ii].initial_value = energies[ii].second;
        initialpars.ampl[ii].min = energies[ii].second * 0.9;
        initialpars.ampl[ii].max = energies[ii].second * 1.1;
        initialpars.ampl[ii].fixed = true;

        initialpars.mean[ii].initial_value = energies[ii].first;
        initialpars.mean[ii].min = energies[ii].first - 1.3;
        initialpars.mean[ii].max = energies[ii].first + 1.3;
        initialpars.mean[ii].fixed = true;

        initialpars.sigma[ii].initial_value = sigmas[ii];
        initialpars.sigma[ii].min = initialpars.min_sigma;
        initialpars.sigma[ii].max = initialpars.max_sigma;
        initialpars.sigma[ii].fixed = false;

        initialpars.tau[ii].initial_value = 0.1;
        initialpars.tau[ii].min = initialpars.min_tau;
        initialpars.tau[ii].max = initialpars.max_tau;
        initialpars.tau[ii].fixed = false;
    }

    int max_idx = 0;
    for (unsigned long int ii = 0; ii < heights.size(); ++ii) {
        if (heights[ii] > heights[max_idx])
            max_idx = ii;
    }

    initialpars.common_ampl.initial_value = heights[max_idx] / energies[max_idx].second * 1.2;
    initialpars.common_ampl.min = std::max(0., initialpars.common_ampl.initial_value * 0.8);
    initialpars.common_ampl.max = std::max(1., initialpars.common_ampl.initial_value * 3);
    initialpars.common_ampl.fixed = false;

    initialpars.common_offset.initial_value = 0;
    initialpars.common_offset.min = 0.;
    initialpars.common_offset.max = heights[max_idx] * 1E-1;
    initialpars.common_offset.fixed = true;

    parameters = initialpars;
}

std::vector<Fitter::FitRes> Fitter::FindFunction() {
    int number_of_gaussians = energies.size();
    const int npars{4};
    fitfunc = TF1(Form("gaussians_%ito%i", (int)energies.front().first, (int)energies.back().first),
                        [number_of_gaussians](const double *x, const double *par) {
                            long double result{0};
                            for(int ii=0; ii<number_of_gaussians; ++ii){
                                double argument = 1./sqrt(2.)*(par[2+ii*npars]/par[3+ii*npars]+(x[0]-par[1+ii*npars])/par[2+ii*npars]);
                                result +=   par[0+ii*npars] *
                                            exp(-0.5*pow((x[0]-par[1+ii*npars])/par[2+ii*npars],2))*
                                            0.5/par[3+ii*npars]*
                                            exp(pow(argument,2))*erfc(argument);
                            };
                            return static_cast<double>(result*par[number_of_gaussians*npars]+par[number_of_gaussians*npars+1]);
                        },
                        parameters.left_interval,
                        parameters.right_interval,
                    npars*number_of_gaussians+2,
                    1);

    fitfunc.SetNpx(300);
    parameters.SetupParameters(fitfunc);

    TFitResultPtr fit_res_ptr(0);
    std:: cout << "===========> Multi fit 1:\n";
    fit_res_ptr = spec.Fit( &fitfunc,
                            "SR",
                            "",
                            fitfunc.GetXmin(),
                            fitfunc.GetXmax());
    parameters.GetParameters(fitfunc);

    TFitResult fit_result = *fit_res_ptr;
    auto cov = fit_result.GetCovarianceMatrix();
    bool fit_is_valid  = fit_result.IsValid();

    if(!fit_is_valid){
        fit_res_ptr = 0;
        std:: cout << "===========> Multi fit 1:\n";
        fit_res_ptr = spec.Fit( &fitfunc,
                                "SRWW",
                                "",
                                fitfunc.GetXmin(),
                                fitfunc.GetXmax());
        parameters.GetParameters(fitfunc);

        fit_result = *fit_res_ptr;
        cov = fit_result.GetCovarianceMatrix();
        fit_is_valid  = fit_result.IsValid();
    }




    if(!fit_is_valid){//2nd try -> Try to increase tau
        for(int ii=0; ii<number_of_gaussians; ++ii) {
            //Releasing mean to account for shift due to tail
            parameters.mean[ii].fixed = false;

            //Trying to increase tail starting point
            parameters.tau[ii].initial_value = 3.;
            parameters.tau[ii].min = 0.06;
            parameters.tau[ii].max = 8.;

        }
        //parameters.common_ampl.min = std::max(0.,parameters.common_ampl.initial_value*0.5);
        //parameters.common_ampl.max = std::max(1.,parameters.common_ampl.initial_value*5.);

        std::cout << "===========> Multi fit 2, trying fit AGAIN:\n";
        parameters.SetupParameters(fitfunc);
        fit_res_ptr = spec.Fit(&fitfunc,
                               "SR",
                               "",
                               fitfunc.GetXmin(),
                               fitfunc.GetXmax());
        parameters.GetParameters(fitfunc);

        fit_result = *fit_res_ptr;
        cov = fit_result.GetCovarianceMatrix();
        fit_is_valid = fit_result.IsValid();

        for (int ii=23; (ii>5)&&(!fit_is_valid); --ii) {

            parameters.left_interval = energies.front().first - ii;
            parameters.right_interval = energies.back().first + ii;

            std::cout << "===========> Multi fit 2, trying fit AGAIN:\n";
            parameters.SetupParameters(fitfunc);
            fit_res_ptr = spec.Fit(&fitfunc,
                                   "SR",
                                   "",
                                   fitfunc.GetXmin(),
                                   fitfunc.GetXmax());
            parameters.GetParameters(fitfunc);

            fit_result = *fit_res_ptr;
            cov = fit_result.GetCovarianceMatrix();
            fit_is_valid = fit_result.IsValid();

            if (!fit_is_valid){
                fit_res_ptr = spec.Fit(&fitfunc,
                                       "SRWW",
                                       "",
                                       fitfunc.GetXmin(),
                                       fitfunc.GetXmax());
                parameters.GetParameters(fitfunc);

                fit_result = *fit_res_ptr;
                cov = fit_result.GetCovarianceMatrix();
                fit_is_valid = fit_result.IsValid();
            }
        }
    }

    if (fit_is_valid) {//Try to improve precision
        std:: cout << "===========> CONVERGED!, trying to improve precision:\n";
        for(int ii=0; ii<number_of_gaussians; ++ii) {
            fitfunc.SetParLimits(1+ii*npars, fitfunc.GetParameter(1+ii*npars)*0.8,fitfunc.GetParameter(1+ii*npars)*1.2);
            fitfunc.SetParLimits(2+ii*npars, fitfunc.GetParameter(2+ii*npars)*0.8,fitfunc.GetParameter(2+ii*npars)*1.2);
            fitfunc.SetParLimits(3+ii*npars, fitfunc.GetParameter(3+ii*npars)*0.8,fitfunc.GetParameter(3+ii*npars)*1.2);

            //fitfunc.SetParameter(2+ii*npars, abs(fitfunc.GetParameter(2+ii*npars)));
        }
        fitfunc.SetParLimits(number_of_gaussians*npars, fitfunc.GetParameter(number_of_gaussians*npars)*0.8,fitfunc.GetParameter(number_of_gaussians*npars)*1.2);
        fit_res_ptr = spec.Fit( &fitfunc,
                                "SMERI",
                                "",
                                fitfunc.GetXmin(),
                                fitfunc.GetXmax());
        parameters.GetParameters(fitfunc);
        auto new_fit_result = *fit_res_ptr;
        if (new_fit_result.IsValid()) {
            fit_result = new_fit_result;
            cov = fit_result.GetCovarianceMatrix();
            fit_is_valid = fit_result.IsValid();
        }
    }

    //Computing parameters and errors
    if (!fit_is_valid) std::cerr << "Fit is not valid\n";
    std::vector<FitRes> results;

    for(int ii=0; ii<number_of_gaussians; ++ii) {
        double integral = parameters.ampl[ii].fitted_value.first * parameters.common_ampl.fitted_value.first;
        double integral_err = sqrt( pow(parameters.ampl[ii].fitted_value.first * parameters.common_ampl.fitted_value.second ,2)+
                                    pow(parameters.ampl[ii].fitted_value.second * parameters.common_ampl.fitted_value.first ,2));

        double sigma = sqrt(pow(parameters.sigma[ii].fitted_value.first,2)+pow(parameters.tau[ii].fitted_value.first,2));
        double sigma_err = 1./sigma*
                           sqrt(    pow(parameters.sigma[ii].fitted_value.first*parameters.sigma[ii].fitted_value.second,2)+
                                    pow(parameters.tau[ii].fitted_value.first*parameters.tau[ii].fitted_value.second,2)+
                                    2*parameters.sigma[ii].fitted_value.first*parameters.tau[ii].fitted_value.first*cov[3+npars*ii][2+npars*ii]);

        double sigma_gauss = parameters.sigma[ii].fitted_value.first;
        double sigma_gauss_err = parameters.sigma[ii].fitted_value.second;

        double tail = parameters.tau[ii].fitted_value.first;
        double tail_err = parameters.tau[ii].fitted_value.second;

        results.emplace_back(
                energies[ii].first,
                integral, integral_err,
                sigma, sigma_err,
                sigma_gauss, sigma_gauss_err,
                tail, tail_err);
    }
    return results;
}


Fitter::~Fitter() = default;

void Fitter::WriteParsOnFile(const std::string& file_name, const int& idx) {
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

        parfile << std::scientific << parameters.sigma[ii].name << "\t\t\t"
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

    parfile << std::scientific << "min_max" << "\t\t\t"
            << std::scientific << parameters.left_interval << "\t\t"
            << std::scientific << parameters.right_interval << std::endl;

    parfile.close();
}

Fitter::InitializationParameters Fitter::ReadParsFromFile(const std::string & parfile_name, const int &idx) {
    std::ifstream parfile(parfile_name);
    if (!parfile.is_open())
        throw std::runtime_error("Parfile not opened\n");

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
    return tmp_parameters;
}
