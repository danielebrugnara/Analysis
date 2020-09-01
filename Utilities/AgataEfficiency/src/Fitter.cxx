#include "Fitter.h"

Fitter::Fitter(TH1D spec, std::vector<std::pair<double, double>> energies):
        spec(std::move(spec)),
        energies(std::move(energies)){
}

std::vector<Fitter::FitRes> Fitter::Fit(){

    //Fitting individual gaussians to get starting points
    std::vector<double> heights;
    std::vector<double> sigmas;
    for (const auto& it_en: energies) {
        std:: cout << "===== Single peak fit :" <<  it_en.first << "\n";
        double single_fit_interval {5.0};
        double single_min_sigma {0.5};
        double single_max_sigma {4.0};

        TF1 fitfunc(Form("gaussian_%f", it_en.first),
                    [](const double *x, const double *par) {
                        return par[0] * exp(-pow((x[0] - par[1]) / (par[2]), 2));
                    },
                    it_en.first - single_fit_interval,
                    it_en.first + single_fit_interval,
                    3,
                    1);
        fitfunc.SetNpx(400);

        fitfunc.SetParameter(0, std::max(0.,spec.GetBinContent(spec.GetXaxis()->FindBin(it_en.first))));
        fitfunc.SetParLimits(0, std::max(0.,fitfunc.GetParameter(0) / 4.), std::max(fitfunc.GetParameter(0) * 3., 40.));

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

    double fit_interval {14.};
    int number_of_gaussians = energies.size();
    const int npars {4};
    double lim_inf = energies.front().first - fit_interval;
    double lim_sup = energies.back().first + fit_interval;
    double min_sigma {0.5};
    double max_sigma {4.0};

    fitfunc = TF1(Form("gaussians_%ito%i", (int)energies.front().first, (int)energies.back().second),
                        [number_of_gaussians](const double *x, const double *par) {
                            long double result{0};
                            for(int ii=0; ii<number_of_gaussians; ++ii){
                                double argument = 1./sqrt(2.)*(par[2+ii*npars]/par[3+ii*npars]+(x[0]-par[1+ii*npars])/par[2+ii*npars]);
                                result +=   par[0+ii*npars] *
                                            exp(-0.5*pow((x[0]-par[1+ii*npars])/par[2+ii*npars],2))*
                                            0.5/par[3+ii*npars]*
                                            exp(pow(argument,2))*erfc(argument);
                            };
                            return static_cast<double>(result*par[number_of_gaussians*npars]);
                        },
                        lim_inf,
                        lim_sup,
                    npars*number_of_gaussians+1,
                    1);

    fitfunc.SetNpx(200);

    for(int ii=0; ii<number_of_gaussians; ++ii) {
        //fitfunc.SetParameter(1+npars*ii, energies[ii].first);
        fitfunc.FixParameter(0+npars*ii, energies[ii].second);
        fitfunc.SetParName(0+npars*ii, Form("IntAmpl_%i", ii));

        fitfunc.SetParameter(1+npars*ii, energies[ii].first);
        fitfunc.SetParLimits(1+npars*ii, fitfunc.GetParameter(1+npars*ii)-2, fitfunc.GetParameter(1+npars*ii)+2);
        fitfunc.FixParameter(1+npars*ii,fitfunc.GetParameter(1+npars*ii));
        fitfunc.SetParName(1+npars*ii, Form("MeanGaus_%i", ii));

        fitfunc.SetParameter(2+npars*ii, sigmas[ii]);
        fitfunc.SetParLimits(2+npars*ii, min_sigma, max_sigma);
        fitfunc.SetParName(2+npars*ii, Form("sigmaGaus_%i", ii));

        fitfunc.SetParameter(3+npars*ii, 0.08);
        fitfunc.SetParLimits(3+npars*ii, 0.04, 5*max_sigma);
        fitfunc.SetParName(3+npars*ii, Form("ExpoDecay_%i", ii));
    }

    int max_idx = 0;
    for (unsigned long int ii=0; ii<heights.size(); ++ii){
        if (heights[ii]>heights[max_idx])
            max_idx = ii;
    }

    fitfunc.SetParameter(number_of_gaussians*npars, heights[max_idx]/energies[max_idx].second*1.2);
    fitfunc.SetParName(number_of_gaussians*npars, "TotalAmpl");
    fitfunc.SetParLimits(number_of_gaussians*npars,
                         std::max(0.,fitfunc.GetParameter(number_of_gaussians*npars)*0.8),
                         std::max(1.,fitfunc.GetParameter(number_of_gaussians*npars)*7.) );


    TFitResultPtr fit_res_ptr(0);
    std:: cout << "===========> Multi fit 1:\n";
    fit_res_ptr = spec.Fit( &fitfunc,
                            "SR",
                            "",
                            fitfunc.GetXmin(),
                            fitfunc.GetXmax());

    TFitResult fit_result = *fit_res_ptr;
    auto cov = fit_result.GetCovarianceMatrix();
    bool fit_is_valid  = fit_result.IsValid();

    if (fit_is_valid) {//Try to improve precision
        std:: cout << "===========> Multi fit 2, trying to improve precision:\n";
        fit_res_ptr = spec.Fit( &fitfunc,
                                "SMERI",
                                "",
                                fitfunc.GetXmin(),
                                fitfunc.GetXmax());
        fit_result = *fit_res_ptr;
        cov = fit_result.GetCovarianceMatrix();
        fit_is_valid  = fit_result.IsValid();
    }else{//change parameters
        for(int ii=0; ii<number_of_gaussians; ++ii) {
            //Releasing mean to account for shift due to tail
            fitfunc.ReleaseParameter(1+npars*ii);

            //Trying to increase tail starting point
            fitfunc.SetParameter(3+npars*ii, 5.);
            fitfunc.SetParLimits(3+npars*ii, 1, 8);
        }
        fitfunc.SetParLimits(number_of_gaussians*npars,
                             std::max(0.,fitfunc.GetParameter(number_of_gaussians*npars)*0.5),
                             std::max(1.,fitfunc.GetParameter(number_of_gaussians*npars)*5.) );

        std:: cout << "===========> Multi fit 2, trying fit AGAIN:\n";
        fit_res_ptr = spec.Fit( &fitfunc,
                                "S",
                                "",
                                fitfunc.GetXmin(),
                                fitfunc.GetXmax());

        fit_result = *fit_res_ptr;
        cov = fit_result.GetCovarianceMatrix();
        fit_is_valid  = fit_result.IsValid();
    }

    //Computing parameters and errors
    if (!fit_is_valid) std::cerr << "Fit is not valid\n";
    std::vector<FitRes> results;

    for(int ii=0; ii<number_of_gaussians; ++ii) {
        double integral =   fitfunc.GetParameter(0+npars*ii)*fitfunc.GetParameter(number_of_gaussians*npars);
        double integral_err =  fitfunc.GetParameter(0+npars*ii)*fitfunc.GetParError(number_of_gaussians*npars);

        double sigma = sqrt(pow(fitfunc.GetParameter(2+npars*ii),2)+pow(fitfunc.GetParameter(3+npars*ii),2));
        double sigma_err = 1./sigma*
                           sqrt(pow(fitfunc.GetParameter(2+npars*ii)*fitfunc.GetParError(2+npars*ii),2)+
                                pow(fitfunc.GetParameter(3+npars*ii)*fitfunc.GetParError(3+npars*ii),2)+
                                2*fitfunc.GetParameter(2+npars*ii)*fitfunc.GetParameter(3+npars*ii)*cov[3+npars*ii][2+npars*ii]);
        //0);

        double sigma_gauss = fitfunc.GetParameter(2+npars*ii);
        double sigma_gauss_err = fitfunc.GetParError(2+npars*ii);

        double tail = fitfunc.GetParameter(3+npars*ii);
        double tail_err = fitfunc.GetParError(3+npars*ii);

        results.emplace_back(
                energies[ii].first,
                integral, integral_err,
                sigma, sigma_err,
                sigma_gauss, sigma_gauss_err,
                tail, tail_err);
    }

    return results;
}

Fitter::~Fitter() {
}
