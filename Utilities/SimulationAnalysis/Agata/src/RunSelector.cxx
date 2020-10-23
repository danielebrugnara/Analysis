#include <utility>

#include "RunSelector.h"



RunSelector::RunSelector(TTree& tree, std::string  out_file_name):
                nevts(0),
                the_tree(tree),
                out_file_name(std::move(out_file_name)){}

void RunSelector::Run(const std::vector<double>& energies, const std::string& fit_histogram_name) {
    nevts = the_tree.GetEntries();
    auto *selector = new Selector();
    auto nentries = TVirtualTreePlayer::kMaxEntries;
    //int nentries = 10;
    auto firstentry = 0;
    the_tree.Process(selector, "option", nentries,firstentry);
    if (nentries<nevts) nevts = nentries;

    TIter iter (selector->GetOutputList());
    TObject *obj;
    auto* out_file = new TFile(out_file_name.c_str(), "recreate");

    std::string peak_search_histo_name = "addb_spec_DC_scan_";
    std::vector<TGraph2D> centroids;
    std::vector<TGraph2D> sigmas;
    std::vector<TGraph2D> sigmas_tot;
    std::vector<TGraph2D> taus;
    for (double energy : energies){
        centroids.emplace_back();
        centroids.back().SetName(Form("centroids_%i", (int)energy));
        centroids.back().SetTitle(Form("centroids_%i", (int)energy));
        sigmas.emplace_back();
        sigmas.back().SetName(Form("sigmas_%i", (int)energy));
        sigmas.back().SetTitle(Form("sigmas_%i", (int)energy));
        sigmas_tot.emplace_back();
        sigmas_tot.back().SetName(Form("sigmas_tot_%i", (int)energy));
        sigmas_tot.back().SetTitle(Form("sigmas_tot_%i", (int)energy));
        taus.emplace_back();
        taus.back().SetName(Form("taus_%i", (int)energy));
        taus.back().SetTitle(Form("taus_%i", (int)energy));
    }

    while ((obj=iter())){
        obj->Write();
        if (fit_histogram_name==obj->GetName() && !energies.empty()){
            GetIntegral(obj, energies);
        }
        std::string tmp_histo_name = obj->GetName();
        if (tmp_histo_name.find(peak_search_histo_name)==0){
            auto results = SearchCentroid(obj, energies);

            std::string tmp_str = tmp_histo_name.substr(tmp_histo_name.find("pos")+3);
            tmp_str = tmp_str.substr(0, tmp_str.find('_'));
            double position = std::stod(tmp_str);

            double beta = std::stod(tmp_histo_name.substr(tmp_histo_name.find("beta")+4));

            for (unsigned int i=0; i<energies.size(); ++i){
                centroids.at(i).SetPoint(centroids.at(i).GetN(),
                                         position,
                                         beta,
                                         results.at("centroids").at(i));

                sigmas.at(i).SetPoint(sigmas.at(i).GetN(),
                                      position,
                                      beta,
                                      results.at("sigmas").at(i));

                sigmas_tot.at(i).SetPoint(sigmas_tot.at(i).GetN(),
                                          position,
                                          beta,
                                          sqrt(pow(results.at("sigmas").at(i), 2)+pow(results.at("taus").at(i),2)));

                taus.at(i).SetPoint(taus.at(i).GetN(),
                                    position,
                                    beta,
                                    results.at("taus").at(i));
            }
        }
    }
    for(const auto& it:centroids){
        it.Write();
    }
    for(const auto& it:sigmas){
        it.Write();
    }
    for(const auto& it:sigmas_tot){
        it.Write();
    }
    for(const auto& it:taus){
        it.Write();
    }
    TVector2 start_stop(0.,nevts);
    start_stop.Write("start_stop");
    out_file->Close();
    delete selector;
}

std::unordered_map<std::string, std::vector<double>> RunSelector::SearchCentroid(TObject *obj, const std::vector<double>& energies) {
    bool read_pars_from_file = false;
    std::string pars_file_name = "tmp_pars.txt";
    TH1D* histo = dynamic_cast<TH1D *>(obj);
    std::vector<double> centroids;
    std::vector<double> sigmas;
    std::vector<double> taus;

    for(const auto& energy: energies) {
        std::vector<std::pair<double, double>> eners;
        eners.emplace_back(energy, 0);
        Fitter fit(*histo, eners, true);
        fit.SetMeanIntervalAroundCenter(20, 20);
        fit.SetTauInterval(1, 10);
        //fit.EnableCanvas();
        fit.SetMeanFixed(false);
        if (read_pars_from_file)
            fit.ReadParsFromFile(pars_file_name, peak_search_cnt);
        else {
            fit.WriteParsOnFile(pars_file_name, peak_search_cnt);
        }
        auto results = fit.Fit();
        centroids.push_back(results.front().mean_gauss.first);
        sigmas.push_back(results.front().sigma_gauss.first);
        taus.push_back(results.front().tail.first);
        ++peak_search_cnt;
    }
    std::unordered_map<std::string,std::vector<double>> results;
    results.emplace("centroids", centroids);
    results.emplace("sigmas", sigmas);
    results.emplace("taus", taus);
    return results;
}

void RunSelector::GetIntegral(TObject *obj, const std::vector<double> &energies) {
    double eff = 0;
    std::ofstream out_file("Integrals.txt", std::ios_base::app);
    for(const auto& energy: energies) {
        std::vector<std::pair<double, double>> eners;
        eners.emplace_back(energy, 0);
        auto* histo = dynamic_cast<TH1D *>(obj);
        if (histo == nullptr)
            throw std::runtime_error("Unable to find histogram\n");
        fits.emplace_back(*histo, eners, true);
        fits.back().SetMeanIntervalAroundCenter(20, 20);
        fits.back().SetTauInterval(1, 10);
        fits.back().EnableCanvas();
        fits.back().SetMeanFixed(false);
        std::vector<Fitter::FitRes> results = fits.back().Fit();
        eff = results.at(0).integral.first / (double) nevts;
        std::cout   << "\033[1;31m---------->Efficiency for " << energy
                    <<" keV : " << eff << "\033[0m\n";
        if (out_file.is_open())
            out_file    << "---->" << obj->GetName()
                        << " @" << energy
                        << " : " << eff << std::endl;
    }

}
