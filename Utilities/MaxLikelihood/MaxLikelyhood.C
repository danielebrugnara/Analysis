#include <vector>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TFile.h>
#include <TH2D.h>
#include <TH1D.h>
#include <TF1.h>
#include <TGraph.h>
#include <TGraph2D.h>
#include <TObjArray.h>
#include <stdio.h>
#include <TMath.h>
#include <TLine.h>
#include <math.h>
#include <functional>

#include "CreateSimulationFile.C"
#include "AverageDistributions.C"

double binomial(int N, int k, double p);

void drawSigmas(TH2D* histo, double, double, TVirtualPad*, double, int);

void computeChiSquared(TH1D* data, TH1D* simu);

long double factorial(int n, std::vector<long double>& mem);

struct LhResults{
    TH2D* lhHisto{nullptr};
    TH2D* chi2Histo{nullptr};
    TH2D* lhHisto2D{nullptr};
    TH2D* chi2Histo2D{nullptr};
    int chi2ndof;
    TH1D* fit{nullptr};
    TH1D* fit2D{nullptr};
    std::vector<TH1D*> fitComponents;
    std::vector<TH1D*> fitComponents2D;
    std::vector<std::pair<double, TH1D*>> fitSigmas;
    TH1D* discrepancy{nullptr};
    TH1D* discrepancy2D{nullptr};
    TH1D* dataWithErrors{nullptr};
    double maxVal{-1E10};
    double minValchi2{1E10};
    double maxX{0};
    double maxY{0};
    double maxXchi2{0};
    double maxYchi2{0};
    double percentX{0};
    double percentY{0};
    double maxVal2D{-1E10};
    double maxX2D{0};
    double maxY2D{0};
    double percentX2D{0};
    double percentY2D{0};
    int N;
};

struct LhInputs{
    TH1D* data{nullptr};
    TH2D* data2D{nullptr};
    std::vector<TH1D*> simu;
    std::vector<TH2D*> simu2D;
    int nX{500};
    //int nX{500};
    std::vector<double> sigmas;
    std::vector<double> xsections;
    int nY{500};
    double startX{0};
    double endX{1};
    double startY{0};
    double endY{2.5};
    //int nY{500};
    //double startX{0};
    //double endX{1};
    //double startY{0.4};
    //double endY{1.4};
    //std::vector<double (*)(double, double )> tsf;
    std::vector<std::function<double(double, double)>> tsf;
    int N;
    int NSimu;
    bool compute2D{false};
    LhInputs(const int& N, const int& NSimu): N(N), NSimu(NSimu){};
};

LhResults gridSearch(LhInputs& inputs);

LhResults MaximizeLikelyhood(TH1D* data, std::vector<TH1D*> simu, TH2D* data2D, std::vector<TH2D*> simu2D, std::map<std::string, double> xsections) {

    if (simu.size() != 3) throw std::runtime_error("not correct simu size\n");

    //std::vector<int> excludedBins {0, 11, 12, 13, 14};
    //std::vector<int> excludedBins {6,7,8};
    //std::vector<int> excludedBins {0};
    //std::vector<int> excludedBins{0, 9, 10, 11, 12, 13, 14, 15};

    //std::map<int, std::vector<int>>nbinstoexcluded;
    //
    //nbinstoexcluded[70] = {0,16, 17, 18, 19, 20, 21, 22};//70 bins
    //nbinstoexcluded[80] = {0,21, 22, 23, 24};//80 bins
    //nbinstoexcluded[90] = {0,21, 22, 23, 24, 25, 26, 27, 28};//90 bins
    //nbinstoexcluded[100] ={0,26, 27, 28, 29, 30, 31, 79};//100 bins
    //nbinstoexcluded[110] ={0,26, 27, 28, 29, 30, 31, 32, 33, 34};//110 bins
    //nbinstoexcluded[120] ={0,26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37};//120 bins
    //nbinstoexcluded[130] ={0,26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40};//130 bins
    //
    //std::vector<int> excludedBins = nbinstoexcluded[data->GetXaxis()->GetNbins()];
    //if(excludedBins.empty()) throw std::runtime_error("decided bins exluded");
    std::vector<int> excludedBins = {0};
    
    //std::vector<int> excludedBins{0, 16, 17, 18, 19, 20, 21, 22};//70 bins
    //std::vector<int> excludedBins{0,21, 22, 23, 24};//80 bins
    //std::vector<int> excludedBins{0,21, 22, 23, 24, 25, 26, 27, 28};//90 bins
    //std::vector<int> excludedBins{0,26, 27, 28, 29, 30, 31, 79};//100 bins
    //std::vector<int> excludedBins{0,26, 27, 28, 29, 30, 31, 32, 33, 34};//110 bins
    //std::vector<int> excludedBins{0,26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37};//120 bins
    //std::vector<int> excludedBins{0,26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40};//130 bins
    std::vector<double> excludedBinsCenter;

    //remove not wanted bins and normalize

    //Exclude from data
    for (const auto &it: excludedBins) {
        data->SetBinContent(it, 0);
        excludedBinsCenter.push_back(data->GetBinCenter(it));
    }
    //Exclude from data2D
    for (const auto &it: excludedBinsCenter) {
        for (int i{0}; i < data2D->GetYaxis()->GetNbins()+2; ++i) {
            data2D->SetBinContent(data2D->GetXaxis()->FindBin(it), i, 0);
            for(int k{0}; k<simu2D.size(); ++k){
                simu2D[k]->SetBinContent(simu2D[k]->GetXaxis()->FindBin(it), i,  0);//too many counts in empty bin??!!
            }
        }
    }

    //Check
    TH1D* histoCheck = data2D->ProjectionX("histoCheck");
    for (int j = 0; j < data->GetNbinsX(); ++j) {//Start from 1 since 0 is the underflow bin!
        if (histoCheck->GetBinContent(j) != data->GetBinContent(j)) {
            std::cout << "check : " << histoCheck->GetBinContent(j) << std::endl;
            std::cout << "data : " << data->GetBinContent(j) << std::endl;
            auto* cv = new TCanvas();
            histoCheck->SetLineColor(kRed);
            histoCheck->Draw();
            data->Draw("same");
            auto* cv2 = new TCanvas();
            data2D->Draw("colz");
            cv->WaitPrimitive();
            throw std::runtime_error("FATAL ERROR in bin "+std::to_string(j)+"\n");
        }
    }


    //Exclude from simu
    for(int i=0; i<simu.size(); ++i){
        for(const auto&it: excludedBins) {
            simu[i]->SetBinContent(it, 0);//too many counts in empty bin??!!
        }
    }
    std::cout << "Finished removing unwanted bins\n";



    //Computing number of events
    int N{0};
    int N2D{0};
    for (int j = 1; j < data->GetNbinsX(); ++j) {//Start from 1 since 0 is the underflow bin!
        N += data->GetBinContent(j);
    }
    int NSimu = 70000;

    //Fist graph
    LhInputs inputs(N, NSimu);
    inputs.tsf.resize(2, nullptr);
    inputs.tsf[0] = [](double x, double y)->double { return x;};
    inputs.tsf[1] = [](double x, double y)->double { return y;};
    inputs.data = data;
    inputs.data2D = data2D;
    inputs.simu = simu;
    inputs.simu2D = simu2D;
    inputs.sigmas ={0.01,0.02, 003};
    inputs.xsections = {xsections["l0_"], xsections["l2_"], xsections["l3_"]};
    inputs.compute2D = false;
    LhResults results = gridSearch(inputs);
    std::cout << "xmin : " << results.maxX << std::endl;
    std::cout << "ymin : " << results.maxY << std::endl;
    std::cout << "max : " << results.maxVal << std::endl;

    //Fit results
    {//1D
        auto *fitResultsCv = new TCanvas("fitResultsCv", "fitResultsCv");
        fitResultsCv->Divide(2, 1);
        auto *tmpPad = fitResultsCv->cd(1);
        tmpPad->SetLogy();

        results.fit->SetLineColor(kRed);
        results.fit->SetMarkerColor(kRed);
        results.dataWithErrors->SetTitle(Form(";Laboratory Angle [deg]; Counts/%.1f deg", 90./results.dataWithErrors->GetXaxis()->GetNbins()));
        results.discrepancy->SetTitle(Form(";Laboratory Angle [deg]; Discrepancy/%.1f deg", 90./results.discrepancy->GetXaxis()->GetNbins()));
        results.dataWithErrors->SetMarkerColorAlpha(kBlack, 0.8);
        results.dataWithErrors->SetMarkerStyle(8);
        results.discrepancy->SetMarkerStyle(8);
        results.dataWithErrors->SetMarkerSize(0.7);
        results.discrepancy->GetXaxis()->SetRangeUser(116, 172);
        results.dataWithErrors->GetXaxis()->SetRangeUser(116,172);
        results.discrepancy->SetMarkerSize(0.7);
        results.dataWithErrors->SetLineColor(kBlack);
        results.discrepancy->SetLineColor(kBlack);
        results.dataWithErrors->SetLineWidth(2);
        results.discrepancy->SetLineWidth(2);
        results.dataWithErrors->Draw("e1");

        results.fit->SetLineWidth(3);
        results.fit->SetLineColor(kRed+1);
        results.fit->Draw("histo, same");

        results.fitComponents[0]->SetLineColor(kBlue);
        results.fitComponents[1]->SetLineColor(kOrange + 1);
        results.fitComponents[2]->SetLineColor(kGreen+2);

        results.fit->SetMarkerSize(0.04);
        results.fitComponents[0]->SetMarkerSize(0.04);
        results.fitComponents[1]->SetMarkerSize(0.04);
        results.fitComponents[2]->SetMarkerSize(0.04);
        results.fitComponents[0]->SetLineWidth(3);
        results.fitComponents[1]->SetLineWidth(3);
        results.fitComponents[2]->SetLineWidth(3);
        results.fit->SetMarkerSize(0.04);
        results.fitComponents[0]->SetMarkerSize(0.04);
        results.fitComponents[1]->SetMarkerSize(0.04);
        results.fitComponents[2]->SetMarkerSize(0.04);

        //results.fitSigmas[2].second->SetLineStyle(8);

        results.fitComponents[0]->Draw("same, histo");
        results.fitComponents[1]->Draw("same, histo");
        results.fitComponents[2]->Draw("same, histo");
        //results.fitSigmas[2].second->Draw("same, histo");

        TLegend *leg = new TLegend(0.1, 0.1, 0.3, 0.3);
        leg->AddEntry(results.dataWithErrors, "data", "lep");
        leg->AddEntry(results.fit, "fit", "l");
        leg->AddEntry(results.fitComponents[0], "L=0", "l");
        leg->AddEntry(results.fitComponents[1], "L=2", "l");
        leg->AddEntry(results.fitComponents[2], "L=3", "l");
        leg->Draw();

        //Discrepancy plot
        fitResultsCv->cd(2);
        results.discrepancy->Draw("e1");
    }
    if (inputs.compute2D){//2D
        auto *fitResultsCv2D = new TCanvas("fitResultsCv2D", "fitResultsCv2D");
        fitResultsCv2D->Divide(2, 1);
        auto *tmpPad = fitResultsCv2D->cd(1);
        tmpPad->SetLogy();

        results.fit2D->SetLineColor(kRed);
        results.fit2D->SetMarkerColor(kRed);
        results.dataWithErrors->Draw("");
        results.fit2D->Draw("histo, same");

        results.fitComponents2D[0]->SetLineColor(kBlue);
        results.fitComponents2D[1]->SetLineColor(kGreen + 1);
        results.fitComponents2D[2]->SetLineColor(kBlack);

        results.dataWithErrors->SetMarkerSize(3);
        results.fit2D->SetMarkerSize(3);
        results.fitComponents2D[0]->SetMarkerSize(3);
        results.fitComponents2D[1]->SetMarkerSize(3);
        results.fitComponents2D[2]->SetMarkerSize(3);
        results.dataWithErrors->SetMarkerSize(3);
        results.fit2D->SetMarkerSize(3);
        results.fitComponents2D[0]->SetMarkerSize(3);
        results.fitComponents2D[1]->SetMarkerSize(3);
        results.fitComponents2D[2]->SetMarkerSize(3);


        results.fitComponents2D[0]->Draw("same, histo");
        results.fitComponents2D[1]->Draw("same, histo");
        results.fitComponents2D[2]->Draw("same, histo");

        TLegend *leg = new TLegend();
        leg->AddEntry(results.dataWithErrors, "data", "lep");
        leg->AddEntry(results.fit2D, "fit", "l");
        leg->AddEntry(results.fitComponents2D[0], "L=0", "l");
        leg->AddEntry(results.fitComponents2D[1], "L=2", "l");
        leg->AddEntry(results.fitComponents2D[2], "L=3", "l");
        leg->Draw();

        //Discrepancy plot
        fitResultsCv2D->cd(2);
        results.discrepancy2D->Draw();
    }

    //Likelihood plot

    {//1D
        auto *likelihoodCv = new TCanvas("likelihoodCv", "likelihoodCv");
        int nbins = inputs.data->GetXaxis()->GetNbins();
        likelihoodCv->Divide(2, 1);
        TVirtualPad *pad21 = likelihoodCv->cd(1);

        results.lhHisto->SetName("lLH");
        results.chi2Histo->SetName("chi2");
        drawSigmas(results.lhHisto, results.maxX, results.maxY, pad21, 2.,nbins);
        computeChiSquared(results.dataWithErrors, results.fit);

        TVirtualPad *pad22 = likelihoodCv->cd(2);
        drawSigmas(results.chi2Histo, results.maxXchi2, results.maxYchi2, pad22, 1., nbins);

        std::cout << "1D ->Expected gammas @360 from L2: " << N * (1-results.percentX-results.percentY) * 0.027 << std::endl;
        std::cout << "1D ->Expected gammas @360 from L3: " << N * results.percentY * 0.013 << std::endl;
        std::cout << "1D ->Expected gammas @1660 from L3: " << N * results.percentY * 0.0099 << std::endl;
    }
    if (inputs.compute2D){//2D
        auto *likelihoodCv2D = new TCanvas("likelihoodCv2D", "likelihoodCv2D");
        likelihoodCv2D->Divide(2, 1);
        TVirtualPad *pad21 = likelihoodCv2D->cd(1);

        results.lhHisto2D->SetName("lLH2D");
        results.chi2Histo2D->SetName("chi22D");
        drawSigmas(results.lhHisto2D, results.maxX2D, results.maxY2D, pad21, 2., 0);
        computeChiSquared(results.dataWithErrors, results.fit2D);

        TVirtualPad *pad22 = likelihoodCv2D->cd(2);
        drawSigmas(results.chi2Histo2D, results.maxX2D, results.maxY2D, pad22, 1.,0);
    }

    results.N = inputs.N;
    return results;
}

LhResults gridSearch(LhInputs& inputs) {
    std::vector<long double> factorialMem;
    LhResults results;
    TH1D *pVals = (TH1D *) inputs.data->Clone("pVals");
    TH2D *pVals2D = (TH2D *) inputs.data2D->Clone("pVals2D");
    results.lhHisto = new TH2D("logLhPlot", "- Log L;C^{2}S[L=2]/C^{2}S[L=0]; C^{2}S[L=3]/C^{2}S[L=0]", inputs.nX, inputs.startX,
                               inputs.endX, inputs.nY, inputs.startY, inputs.endY);
    results.chi2Histo = new TH2D("chi2Plot", "#Chi^{2};C^{2}S[L=2]/C^{2}S[L=0]; C^{2}S[L=3]/C^{2}S[L=0]", inputs.nX, inputs.startX,
                                 inputs.endX, inputs.nY, inputs.startY, inputs.endY);
    results.lhHisto2D = new TH2D("logLhPlot2D", "- Log L;C^{2}S[L=2]/C^{2}S[L=0]; C^{2}S[L=3]/C^{2}S[L=0]", inputs.nX,
                                 inputs.startX, inputs.endX, inputs.nY, inputs.startY, inputs.endY);
    results.chi2Histo2D = new TH2D("chi2Plot2D", "#Chi^{2};C^{2}S[L=2]/C^{2}S[L=0]; C^{2}S[L=3]/C^{2}S[L=0]", inputs.nX, inputs.startX,
                                   inputs.endX, inputs.nY, inputs.startY, inputs.endY);

    results.chi2ndof = 0;
    for (int i = 1; i < inputs.data->GetNbinsX(); ++i) {
        if (inputs.data->GetBinContent(i)) results.chi2ndof++;
    }

    new TCanvas();
    TH1D *projEx = inputs.data2D->ProjectionY();
    projEx->Draw();
    TF1 *exPeaks = new TF1("exPeaks", "gaus(0)+gaus(3)", -5, 5);
    exPeaks->FixParameter(1, 0);
    exPeaks->FixParameter(4, 2.02);
    exPeaks->SetParameter(0, 50);
    exPeaks->SetParameter(3, 50);
    exPeaks->SetParameter(2, 1.2);
    exPeaks->SetParameter(5, 1.2);

    projEx->Fit(exPeaks, "RI+");
    exPeaks->Draw("same");
    std::vector<std::pair<double, double>> energies{{0,    abs(exPeaks->GetParameter(2))},
                                                    {0.36, abs(exPeaks->GetParameter(2))},
                                                    {2.02, abs(exPeaks->GetParameter(5))}};
    std::vector<TF1 *> exFunctions;
    std::vector<TH1D *> exHistograms;
    for (const auto &it: energies) {
        exFunctions.push_back(new TF1(("ex" + std::to_string(it.first)).c_str(), "gaus(0)", -100, 100));
        exFunctions.back()->FixParameter(0, 0.387);
        exFunctions.back()->FixParameter(1, it.first);
        exFunctions.back()->FixParameter(2, it.second);
        exFunctions.back()->SetNpx(1000);
        exHistograms.push_back(new TH1D(("exhisto" + std::to_string(it.first)).c_str(),
                                        ("exhisto" + std::to_string(it.first)).c_str(),
                                        inputs.data2D->GetYaxis()->GetNbins(),
                                        inputs.data2D->GetYaxis()->GetXmin(),
                                        inputs.data2D->GetYaxis()->GetXmax()));
        double samples{1E6};
        for (int i{0}; i < samples; ++i) {
            exHistograms.back()->Fill(exFunctions.back()->GetRandom(inputs.data2D->GetYaxis()->GetXmin(),
                                                                    inputs.data2D->GetYaxis()->GetXmax(),
                                                                    new TRandom1()));
        }
        exHistograms.back()->Scale(1. / samples);
        //new TCanvas(); //Plot for debug
        //exHistograms.back()->Draw();
        //exFunctions.back()->Draw("same");
    }


    std::vector<double> coeffs;
    coeffs.resize(3, 0);
    results.fitSigmas.resize(inputs.sigmas.size(), {0, nullptr});

    results.fitComponents.resize(3);
    results.fitComponents[0] = nullptr;
    results.fitComponents[1] = nullptr;
    results.fitComponents[2] = nullptr;
    results.fitComponents2D.resize(3);
    results.fitComponents2D[0] = nullptr;
    results.fitComponents2D[1] = nullptr;
    results.fitComponents2D[2] = nullptr;

    std::vector<TH1D> simuNormalized;
    std::vector<TH2D> simuNormalized2D;
    //Normalize simulations
    for (unsigned long k = 0; k < coeffs.size(); ++k) {//Loop on L transfers
        std::string name = inputs.simu[k]->GetName();
        simuNormalized.push_back(*(static_cast<TH1D *>(inputs.simu[k]->Clone(name.c_str()))));
        simuNormalized2D.push_back(*(static_cast<TH2D *>(inputs.simu2D[k]->Clone((name+"2D").c_str()))));
        simuNormalized.back().SetBinContent(0, 0); //Remove underflow, just in case
        double normalization{0};
        double normalization2D{0};
        for (int i = 0; i < simuNormalized.back().GetNbinsX(); ++i) {
            normalization += simuNormalized.back().GetBinContent(i);
        }
        for (int ix = 0; ix < simuNormalized2D.back().GetXaxis()->GetNbins(); ++ix) {
            for (int iy = 0; iy < simuNormalized2D.back().GetYaxis()->GetNbins(); ++iy) {
                normalization2D += simuNormalized2D.back().GetBinContent(ix, iy);
            }
        }
        //Normalize
        for (int i = 0; i < simuNormalized.back().GetNbinsX(); ++i) {
            simuNormalized.back().SetBinContent(i, simuNormalized.back().GetBinContent(i) / normalization);
        }
        for (int ix = 0; ix < simuNormalized2D.back().GetXaxis()->GetNbins(); ++ix) {
            for (int iy = 0; iy < simuNormalized2D.back().GetYaxis()->GetNbins(); ++iy) {
                simuNormalized2D.back().SetBinContent(ix, iy, simuNormalized2D.back().GetBinContent(ix, iy)/normalization);
            }
        }

    }

    inputs.xsections[0] *= inputs.simu[0]->GetEntries()/70000.; 
    inputs.xsections[1] *= inputs.simu[1]->GetEntries()/70000.; 
    inputs.xsections[2] *= inputs.simu[2]->GetEntries()/70000.; 

    std::cout << "--------------L=0 :" <<  inputs.simu[0]->GetEntries()<< std::endl; 
    std::cout << "--------------L=2 :" <<  inputs.simu[1]->GetEntries()<< std::endl; 
    std::cout << "--------------L=3 :" <<  inputs.simu[2]->GetEntries()<< std::endl; 
    std::cout << "\n";

    //Start of grid search
    std::cout << "-----------------------------------------------------------------------\n";
    std::cout << "--------------starting grid search-------------------------------------\n";
    int cnt{0};
    for (int x = 0; x < inputs.nX; ++x) {
        double xval = inputs.startX + (inputs.endX - inputs.startX) / (inputs.nX - 1) * x;
        for (int y = 0; y < inputs.nY; ++y) {
            std::cout << "\riteration " << cnt++ << "/" << inputs.nX * inputs.nY;
            double yval = inputs.startY + (inputs.endY - inputs.startY) / (inputs.nY - 1) * y;
            //std::cout << "new pt\n";
            coeffs[0] = 1;
            coeffs[1] = inputs.tsf[0](xval, yval);
            coeffs[2] = inputs.tsf[1](xval, yval);

            //std::cout   << "(p0 : " << coeffs[0] << " p2 " << coeffs[1] << ")->("
            //            << "x: " << xval << " y: " << yval << ")" << std::endl;

            if (coeffs[0] < 0 || coeffs[1] < 0 || coeffs[2] < 0) {
                results.lhHisto->SetBinContent(results.lhHisto->FindBin(xval, yval), 1);
                results.chi2Histo->SetBinContent(results.chi2Histo->FindBin(xval, yval), 0);
                results.lhHisto2D->SetBinContent(results.lhHisto->FindBin(xval, yval), 1);
                results.chi2Histo2D->SetBinContent(results.chi2Histo->FindBin(xval, yval), 0);
                continue;
            }
            long double logLh{0};
            long double chi2{0};
            long double logLh2D{0};
            long double chi22D{0};



            //Initialize pvals at zero
            for (int j = 0; j < pVals->GetNbinsX(); ++j) {
                pVals->SetBinContent(j, 0);
            }
            pVals2D->Reset("ICE");

            ////////////////////////////////////////////////////////////////////////////////////////////////////////////
            //1 dimensional calculation
            //Sets p values
            for (int j = 0; j < pVals->GetNbinsX(); ++j) {
                double normalization{0};
                double probVal{0};
                for (unsigned long k = 0; k < coeffs.size(); ++k) {//Loop on L transfers
                    //Sigma ration changes normalization
                    probVal += coeffs[k] * inputs.xsections[k] * simuNormalized[k].GetBinContent(j);
                    normalization += coeffs[k] * inputs.xsections[k];
                }
                pVals->SetBinContent(j, probVal / normalization);
            }
            double check{0};
            for (int j = 1; j < pVals->GetNbinsX(); ++j) {
                check += pVals->GetBinContent(j);
            }
            if (abs(check - 1) > 1E-4)
                throw std::runtime_error("Not normalized, value is : " + std::to_string(check) + "\n");

            //Compute Likelihood based on data
            logLh += logl((long double) factorial(inputs.N, factorialMem));
            for (int j = 1; j < pVals->GetNbinsX(); ++j) {//Starts from 1 to skip overflow (bin 0)
                //logLh += logl((long double) binomial(inputs.N, inputs.data->GetBinContent(j), pVals->GetBinContent(j)));
                logLh += logl((long double) powl(pVals->GetBinContent(j), inputs.data->GetBinContent(j)) /
                              factorial((int) inputs.data->GetBinContent(j), factorialMem));
                if (inputs.data->GetBinContent(j) > 0) {
                    double error = sqrt(inputs.N * pVals->GetBinContent(j) * (1 - pVals->GetBinContent(j)));
                    chi2 += pow((inputs.data->GetBinContent(j) - inputs.N * pVals->GetBinContent(j)) / error, 2);
                }
            }
            results.lhHisto->SetBinContent(results.lhHisto->FindBin(xval, yval), -logLh);
            results.chi2Histo->SetBinContent(results.chi2Histo->FindBin(xval, yval), chi2);

            //Save max vals
            if (logLh > results.maxVal) {
                results.maxVal = logLh;
                results.maxX = xval;
                results.maxY = yval;
                double norm = (1*inputs.xsections[0]+xval*inputs.xsections[1]+yval*inputs.xsections[2]);
                results.percentX = xval*inputs.xsections[1]/norm;
                results.percentY = yval*inputs.xsections[2]/norm;
                if (results.fit != nullptr) {
                    results.fit->Delete();
                }
                results.fit = (TH1D *) pVals->Clone("fitResult");
                results.fit->Scale(inputs.N);

                if (results.dataWithErrors != nullptr) {
                    results.dataWithErrors->Delete("");
                }
                results.dataWithErrors = (TH1D *) inputs.data->Clone("DataWithErrors");
                for (int j = 0; j < pVals->GetNbinsX(); ++j) {//Sets multivariate variance as error (approximation)
                    results.dataWithErrors->SetBinError(j, sqrt(inputs.N * pVals->GetBinContent(j) *
                                                                (1 - pVals->GetBinContent(j))));
                }

                for (unsigned int s = 0; s < inputs.sigmas.size(); ++s) {
                    if (abs(yval - inputs.sigmas[s]) < abs(yval - results.fitSigmas[s].first)) {
                        results.fitSigmas[s].first = yval;
                        results.fitSigmas[s].second = (TH1D *) pVals->Clone(Form("%i_Sigma", s));
                        results.fitSigmas[s].second->Scale(inputs.N);
                    }
                }
            }//End of save max vals
            if (chi2 < results.minValchi2){
                results.minValchi2 = chi2;
                results.maxXchi2 = xval;
                results.maxYchi2 = yval;
            }

            ////////////////////////////////////////////////////////////////////////////////////////////////////////////
            //2 dimensional calculation
            //Sets p values
            if (!inputs.compute2D) continue;
            for (int dx = 1; dx < pVals2D->GetNbinsX(); ++dx) {
                for (int dy = 1; dy < pVals2D->GetNbinsY(); ++dy) {
                    double th = pVals2D->GetXaxis()->GetBinCenter(dx);
                    double ex = pVals2D->GetYaxis()->GetBinCenter(dy);
                    double normalization{0};
                    double probVal{0};

                    for (unsigned long k = 0; k < coeffs.size(); ++k) {//Loop on L transfers
                        //Sigma ration changes normalization
                        probVal += coeffs[k] * inputs.xsections[k] *
                                   simuNormalized2D[k].GetBinContent(   simuNormalized2D[k].GetXaxis()->FindBin(th), 
                                                                        simuNormalized2D[k].GetYaxis()->FindBin(ex));
                        normalization += coeffs[k] * inputs.xsections[k];
                    }
                    pVals2D->SetBinContent(dx, dy, probVal / normalization);
                }
            }
            double check2D{0};
            for (int dx = 1; dx < pVals2D->GetNbinsX(); ++dx) {
                for (int dy = 1; dy < pVals2D->GetNbinsY(); ++dy) {
                    check2D += pVals2D->GetBinContent(dx, dy);
                }
            }
            //if (abs(check2D - 1) > 1E-2) ///WARNING!!!!!!!!!!!!!!!!
            //    throw std::runtime_error("Not normalized, value is : " + std::to_string(check2D) + "\n");

            //Compute Likelihood based on data
            logLh2D += logl((long double) factorial(inputs.N, factorialMem));
            //std::cout << "start loglh2d : " << logLh2D << std::endl;
            try {
                for (int dx = 1; dx < pVals2D->GetNbinsX(); ++dx) {
                    for (int dy = 1; dy < pVals2D->GetNbinsY(); ++dy) {
                        double th = pVals2D->GetXaxis()->GetBinCenter(dx);
                        double ex = pVals2D->GetYaxis()->GetBinCenter(dy);
                        long double partialLogLh = logl(
                                (long double) powl(pVals2D->GetBinContent(dx, dy),
                                                   inputs.data2D->GetBinContent(dx, dy)) /
                                factorial((int) inputs.data2D->GetBinContent(dx, dy), factorialMem));
                        logLh2D += partialLogLh;
                        //std::cout << "new loglh2d : " << logLh2D << std::endl;
                        if (std::isnan(logLh2D) || std::isinf(logLh2D)) {
                            if (!std::isinf(partialLogLh)) std::cout << "partial is not inf!!!!!!!!!!!!!!!!!\n";
                            std::cout << "giving parameters pVals2d " << pVals2D->GetBinContent(dx, dy) << std::endl;
                            std::cout << "giving parameters data2d " << inputs.data2D->GetBinContent(dx, dy)
                                      << std::endl;
                            std::cout << "ex value " << ex << std::endl;
                            std::cout << "theta value " << th << std::endl;
                            std::cout << "coeff 1 " << coeffs[1] << std::endl;
                            std::cout << "coeff 2 " << coeffs[2] << std::endl;
                            std::cout << "factorial "
                                      << factorial((int) inputs.data2D->GetBinContent(dx, dy), factorialMem)
                                      << std::endl;
                            //throw std::runtime_error("stop\n"); //WARNING!!!!!!!!!!!!
                        }
                        if (inputs.data2D->GetBinContent(dx, dy) > 0) {
                            double error2D = sqrt(
                                    inputs.N * pVals2D->GetBinContent(dx, dy) * (1 - pVals->GetBinContent(dx, dy)));
                            chi22D += pow(
                                    (inputs.data2D->GetBinContent(dx, dy) - inputs.N * pVals2D->GetBinContent(dx, dy)) /
                                    error2D, 2);
                        }
                    }
                }
            } catch (const std::runtime_error &e) {
                auto *tmpcv = new TCanvas("pvals");
                pVals2D->Draw("colz");
                auto *tmpcv22 = new TCanvas("data");
                inputs.data2D->Draw("colz");
                tmpcv->WaitPrimitive();
                logLh2D = 0;
            }
            results.lhHisto2D->SetBinContent(results.lhHisto->FindBin(xval, yval), -logLh2D);
            results.chi2Histo2D->SetBinContent(results.chi2Histo->FindBin(xval, yval), chi22D);

            //Save max vals
            if (logLh2D > results.maxVal2D) {
                results.maxVal2D = logLh2D;
                results.maxX2D = xval;
                results.maxY2D = yval;
                double norm = (1*inputs.xsections[0]+xval*inputs.xsections[1]+yval*inputs.xsections[2]);
                results.percentX2D = xval*inputs.xsections[1]/norm;
                results.percentY2D = yval*inputs.xsections[2]/norm;
                if (results.fit2D != nullptr) {
                    results.fit2D->Delete();
                }
                results.fit2D = (TH1D *) pVals2D->ProjectionX("fitResult2D");
                results.fit2D->Scale(inputs.N);
            }//End of save max vals
        }//End of y scan
    }//End of x scan

    //Update errors based on multivariate distribution
    results.discrepancy = (TH1D *) results.dataWithErrors->Clone("discrepancies");
    results.discrepancy2D = (TH1D *) results.dataWithErrors->Clone("discrepancies2D");
    for (int i = 0; i < results.discrepancy->GetNbinsX(); ++i) {
        results.discrepancy->SetBinContent(i, results.dataWithErrors->GetBinContent(i) - results.fit->GetBinContent(i));
        results.discrepancy->SetBinError(i, results.dataWithErrors->GetBinError(i));
        if(inputs.compute2D){
            results.discrepancy2D->SetBinContent(i, results.dataWithErrors->GetBinContent(i) - results.fit2D->GetBinContent(i));
            results.discrepancy2D->SetBinError(i, results.dataWithErrors->GetBinError(i));
        }
    }

    //double scale = results.fit->Integral(results.fit->FindBin(90), results.fit->FindBin(180));
    double scale = 2000;

    //Save fot components
    {//1D
        results.fitComponents[0] = static_cast<TH1D *>( simuNormalized[0].Clone(Form("partialL%i", 0)));
        results.fitComponents[1] = static_cast<TH1D *>( simuNormalized[1].Clone(Form("partialL%i", 2)));
        results.fitComponents[2] = static_cast<TH1D *>( simuNormalized[2].Clone(Form("partialL%i", 3)));
        double normalization =
                1. * inputs.xsections[0] + results.maxX * inputs.xsections[1] + results.maxY * inputs.xsections[2];
        //Rescale fit components
        results.fitComponents[0]->Scale(inputs.N * 1. * inputs.xsections[0] / normalization);
        results.fitComponents[1]->Scale(inputs.N * results.maxX * inputs.xsections[1] / normalization);
        results.fitComponents[2]->Scale(inputs.N * results.maxY * inputs.xsections[2] / normalization);
    }
    if (inputs.compute2D){//2D
        results.fitComponents2D[0] = static_cast<TH1D *>( simuNormalized[0].Clone(Form("partialL%i", 0)));
        results.fitComponents2D[1] = static_cast<TH1D *>( simuNormalized[1].Clone(Form("partialL%i", 2)));
        results.fitComponents2D[2] = static_cast<TH1D *>( simuNormalized[2].Clone(Form("partialL%i", 3)));
        double normalization =
                1. * inputs.xsections[0] + results.maxX2D * inputs.xsections[1] + results.maxY2D * inputs.xsections[2];
        //Rescale fit components
        results.fitComponents2D[0]->Scale(inputs.N * 1. * inputs.xsections[0] / normalization);
        results.fitComponents2D[1]->Scale(inputs.N * results.maxX2D * inputs.xsections[1] / normalization);
        results.fitComponents2D[2]->Scale(inputs.N * results.maxY2D * inputs.xsections[2] / normalization);
    }
    return results;
}

void drawSigmas(TH2D* histo, double minx, double miny, TVirtualPad* cv, double sigmaCoeff, int nbins){
    bool divide{true};
    if(divide){
        cv->Divide(1,3);
        cv->cd(1);
    }else{
        cv->cd();
    }
    double val = histo->GetBinContent(histo->FindBin(minx, miny));
    std::string name = histo->GetName();
    name += "_clone";
    auto* histoclone = (TH2D*)histo->Clone(name.c_str());

    std::vector<double> contours;
    for (int i=0; i<100; ++i){
        contours.push_back(val+(double)i*i/(sigmaCoeff));
    }
    histo->SetContour(contours.size(), &contours[0]);

    histo->Draw("CONT Z LIST");
    cv->Update();

    TObjArray *conts = (TObjArray*)gROOT->GetListOfSpecials()->FindObject("contours");
    TList* contLevel = NULL;
    TGraph* curv     = NULL;
    std::vector<TGraph*> gc;
    Int_t nGraphs    = 0;
    Int_t TotalConts = 0;

    if (conts == NULL){
        printf("*** No Contours Were Extracted!\n");
        TotalConts = 0;
        return 0;
    } else {
        TotalConts = conts->GetSize();
    }
    for(int i = 0; i < TotalConts; i++){
        contLevel = (TList*)conts->At(i);
        //printf("Contour %d has %d Graphs\n", i, contLevel->GetSize());
        nGraphs += contLevel->GetSize();
    }
    nGraphs = 0;
    Double_t xval0, yval0, zval0;

    for(int i = 0; i < TotalConts; i++){
        contLevel = (TList*)conts->At(i);
        if (i<3) zval0 = contours[2-i];
        else     zval0 = contours[i];
        //printf("Z-Level Passed in as:  Z = %f\n", zval0);

        // Get first graph from list on curves on this level
        curv = (TGraph*)contLevel->First();
        for(int j = 0; j < contLevel->GetSize(); j++){
            curv->GetPoint(0, xval0, yval0);
            if (zval0<0) curv->SetLineColor(kRed);
            if (zval0>0) curv->SetLineColor(kBlue);
            nGraphs ++;
            //printf("\tGraph: %d  -- %d Elements\n", nGraphs,curv->GetN());
            gc.push_back((TGraph*)curv->Clone());
            curv = (TGraph*)contLevel->After(curv); // Get Next graph
        }
    }

    //histoclone->SetMaximum(0.3);
    //histoclone->SetMinimum(0.01);
    histoclone->Draw("COLZ");

    for (const auto&it: gc){
        //it->Draw("C PMC PLC");
        it->SetLineColor(kRed);
        it->Draw("C");
    }
    auto* minPlot = new TGraph();
    minPlot->SetPoint(0, minx, miny);
    minPlot->SetMarkerSize(1);
    minPlot->SetMarkerColor(kRed);
    minPlot->Draw("same, P*");
    if(sigmaCoeff==2){
        gc[1]->SetName(Form("cont_%i", nbins));
        minPlot->SetName(Form("min_%i", nbins));
        gc[1]->SaveAs(Form("cont_%i_2.root", nbins));
        minPlot->SaveAs(Form("min_%i.root", nbins));
    }

    int x, y, z;
    double nsigma = 3;
    if(divide){
    histoclone->GetBinXYZ(histoclone->FindBin(minx, miny), x, y, z);
    {
        auto* tmpPad = cv->cd(2);
        //cv->Update();
        std::string name = histoclone->GetName();
        name += "_pX";
        auto* proj = histoclone->ProjectionX(name.c_str(), y, y);
        proj->SetMaximum(val+nsigma*nsigma/2.);
        proj->SetMinimum(val-0.5);
        proj->Draw();
        for (int i=0; i<nsigma; ++i){
            double dist = val+i*i/2.;
            auto* line = new TLine(proj->GetXaxis()->GetXmin(), dist, proj->GetXaxis()->GetXmax(),dist );
            line->SetLineColor(kRed);
            line->SetLineWidth(i+1);
            line->Draw();

            //double xLow, xHigh,y,z;
            //int binNr;
            //proj->GetBinWithContent(dist,
            //                              binNr,
            //                              proj->FindBin(proj->GetXaxis()->GetXmin())+0.001,
            //                              proj->FindBin(minx));
            //xLow = proj->GetBinCenter(binNr);

            //proj->GetBinWithContent(dist,
            //                              binNr,
            //                              proj->FindBin(minx),
            //                              proj->FindBin(proj->GetXaxis()->GetXmax()));
            //xHigh = proj->GetBinCenter(binNr);

            //auto* lineXLow = new TLine(xLow, proj->GetYaxis()->GetXmin(), xLow, proj->GetYaxis()->GetXmax());
            //auto* lineXHigh = new TLine(xHigh, proj->GetYaxis()->GetXmin(), xHigh, proj->GetYaxis()->GetXmax());
            //lineXLow->Dump();
            //lineXHigh->Dump();

            //lineXLow->SetLineColor(kRed);
            //lineXLow->SetLineWidth(2);
            //lineXHigh->SetLineColor(kRed);
            //lineXHigh->SetLineWidth(2);
            //lineXLow->Draw();
            //lineXHigh->Draw();
        }
    }
    {
        cv->cd(3);
        //cv->Update();
        std::string name = histoclone->GetName();
        name += "_pY";
        auto* proj = histoclone->ProjectionY(name.c_str(), x, x);
        proj->SetMaximum(val+nsigma*nsigma/2.);
        proj->SetMinimum(val-0.5);
        proj->Draw();
        for (int i=0; i<nsigma; ++i){
            double dist = val+i*i/2.;
            auto* line = new TLine(proj->GetXaxis()->GetXmin(), dist, proj->GetXaxis()->GetXmax(),dist );
            line->SetLineColor(kRed);
            line->SetLineWidth(i+1);
            line->Draw();
        }
    }
    }
    cv->Update();

}

double binomial(int N, int k, double p){
    return TMath::Binomial(N, k) * pow(p, k)* pow(1-p, N-k);
}

std::map<std::string, std::string> SumThicknesses(std::map<std::string, double>& xsections){
    std::map<int, int> thicknessToNDeuterons{
            {15, 62},
            {20, 0},
            {25, 207},
            {30, 151},
            {35, 151},
            {40, 60},
            {45, 41},
            {50, 0},
            {55, 0}
    };
    int maxVal = std::max_element(thicknessToNDeuterons.begin(), thicknessToNDeuterons.end(),
                                  [](const std::pair<int, int>& p1, const std::pair<int, int>& p2){
                                      return p1.second < p2.second;
                                  }
    )->second;

    //Root files sum over different thicknesses
    std::map<std::string, std::string> prefixes{
            {"s12", "selector_46Ar3Hed47K_0keV_s12_"},
            {"d32", "selector_46Ar3Hed47K_360keV_d32_"},
            {"f72", "selector_46Ar3Hed47K_2020keV_f72_"},
            {"flat0", "selector_46Ar3Hed47K_0keV_flat_"},
            {"flat2", "selector_46Ar3Hed47K_360keV_flat_"},
            {"flat3", "selector_46Ar3Hed47K_2020keV_flat_"}
    };

    std::map<int, double> percentOfThickness;
    for (const auto& it: thicknessToNDeuterons){
        percentOfThickness.emplace(it.first, static_cast<double>(it.second)/maxVal);
    }

    std::map<std::string, std::string> outputFiles;
    for (const auto& itReaction:prefixes){
        std::map<std::string, double> files;
        for (const auto& itThickness: percentOfThickness){
            std::string filePath{"./../../DataAnalyzed/simu/"};
            filePath += itReaction.second;
            //if(itReaction.second.find("flat")!=itReaction.second.npos){
            //    filePath += "0_0_";
            //}
            filePath += std::to_string(itThickness.first);
            filePath += "um_analyzed.root";
            files[filePath] = itThickness.second;
        }
        outputFiles.emplace(itReaction.first, itReaction.second+"sumthickness.root");
        std::ifstream testFile(outputFiles[itReaction.first]);
        if (!testFile.is_open())
            CreateSimulationFile(files, outputFiles[itReaction.first]);
    }
    //Angular distributions sum over different thicknesses
    std::string distrFolder{"./AngularDistributions/"};
    std::vector<std::string> distrPrefixes{"l0_", "l2_", "l3_"};
    std::string distrSuffix{"um.in"};

    for(const auto& itDistr: distrPrefixes){
        std::map<std::string, double> filesToAverage;
        for(const auto& itThickness: thicknessToNDeuterons){
            std::string fileName{distrFolder+itDistr+std::to_string(itThickness.first)+distrSuffix};
            filesToAverage[fileName] = itThickness.second;
        }
        xsections[itDistr] = AverageDistributions(filesToAverage, itDistr+"averaged.txt");
    }


    return outputFiles;
}

void SumSimulationsWithPercentages(const LhResults& results, std::map<std::string, std::string> files){
    {
        std::map<std::string, double> filesToSum;
        filesToSum.emplace(files["s12"], 1-results.percentX-results.percentY);
        filesToSum.emplace(files["d32"], results.percentX);
        filesToSum.emplace(files["f72"], results.percentY);
        //double p0 = 0.40*2.48;
        //double p2 = 0.24*2.65;
        //double p3 = 0.7*3.2;
        //double ptot = p0+p2+p3;
        //filesToSum.emplace(files["s12"], p0/ptot);
        //filesToSum.emplace(files["d32"], p2/ptot);
        //filesToSum.emplace(files["f72"], p3/ptot);
        CreateSimulationFile(filesToSum, "simulationsum.root");
    }
    {
        std::map<std::string, double> filesToSum;
        filesToSum.emplace(files["flat0"], 1-results.percentX-results.percentY);
        filesToSum.emplace(files["flat2"], results.percentX);
        filesToSum.emplace(files["flat3"], results.percentY);
        CreateSimulationFile(filesToSum, "simulationsumflat.root");
    }
}

void NormalizeHistograms(const std::string& fileName, const LhResults& results){
    TFile* outFile = new TFile(fileName.c_str(),"update");
    if (!outFile->IsOpen()) throw std::runtime_error("histofile not present!!");
    std::cout << "Normalizing histograms\n";

    double (*integral)(TH1D*) = [](TH1D* h){
        return h->Integral(h->FindBin(1), h->FindBin(89));
    };

    void (*normalize)(TH1D*, TH1D*, double  (*)(TH1D*)) = [](TH1D* h, TH1D* data, double (*integral)(TH1D*)){
        h->Scale(integral(data)/integral(h));
    };

    TH1D* (*cloneHisto)(const std::string&, TFile*) = [](const std::string& name, TFile* outFile){
        return  static_cast<TH1D*>(static_cast<TH1D*>(outFile->Get(name.c_str()))->Clone((name+"Normalized").c_str()));
    };

    auto* dataCM = cloneHisto("dataCM", outFile);
    double N = integral(dataCM);

    std::vector<TH1D* > histos;
    auto* normalizeds12CM = cloneHisto("s12CM", outFile);
    auto* normalizedd32CM = cloneHisto("d32CM", outFile);
    auto* normalizedf72CM = cloneHisto("f72CM", outFile);
    auto* normalizedflat0CM = cloneHisto("flat0CM", outFile);
    auto* normalizedflat2CM = cloneHisto("flat2CM", outFile);
    auto* normalizedflat3CM = cloneHisto("flat3CM", outFile);

    normalize(normalizeds12CM  , dataCM, integral);
    normalize(normalizedd32CM  , dataCM, integral);
    normalize(normalizedf72CM  , dataCM, integral);
    normalize(normalizedflat0CM, dataCM, integral);
    normalize(normalizedflat2CM, dataCM, integral);
    normalize(normalizedflat3CM, dataCM, integral);

    TH1D normalizationTotal = (1-results.percentX-results.percentY)* *normalizedflat0CM + results.percentX* *normalizedflat2CM + results.percentY* *normalizedflat3CM;
    normalizationTotal.SetName("normalizationTotal");

    *dataCM = *dataCM / normalizationTotal;
    TH1D* simuCM = new TH1D(((1-results.percentX-results.percentY)* *normalizeds12CM + results.percentX* *normalizedd32CM + results.percentY* *normalizedf72CM) / normalizationTotal);
    simuCM->SetName("simuCMNormalized");


    normalizeds12CM  ->Write();
    normalizedd32CM  ->Write();
    normalizedf72CM  ->Write();
    normalizedflat0CM->Write();
    normalizedflat2CM->Write();
    normalizedflat3CM->Write();
    dataCM->Write();
    simuCM->Write();
    normalizationTotal.Write();
    outFile->Write();
    outFile->Close();
}


TGraph* gsum{nullptr};
double func(double* x, double* p ){
    return gsum->Eval(x[0]);
}
double funcSin(double* x, double* p ){
    return sin(x[0]) * gsum->Eval(x[0]);
}

void SaveTheoryDistributions(const std::string& fileName, const LhResults& results) {
    TFile *outFile = new TFile(fileName.c_str(), "update");
    if (!outFile->IsOpen()) throw std::runtime_error("histofile not present!!");
    std::cout << "Normalizing histograms\n";

    TGraph gs12("l0_averaged.txt");
    gs12.SetName("gs12");
    TGraph gd32("l2_averaged.txt");
    gd32.SetName("gd32");
    TGraph gf72("l3_averaged.txt");
    gf72.SetName("gf72");

    gs12.Write();
    gd32.Write();
    gf72.Write();

    TGraph gsum;
    gsum.SetName("gsum");

    for (int i{0}; i<gs12.GetN(); ++i){
        gsum.SetPoint(gsum.GetN(),
                      gs12.GetPointX(i),
                      gs12.GetPointY(i)* (1-results.percentX-results.percentY)
                      + gd32.GetPointY(i)* results.percentX
                      + gf72.GetPointY(i) *results.percentY);
    }
    gsum.Write();
    //Correction factor
    //std::function<double(double*, double*)> funcSin = [&gsum](double* x, double* p){return sin(x[0]) * gsum.Eval(x[0]);};
    //std::function<double(double*, double*)> func = [&gsum](double* x, double* p){return gsum.Eval(x[0]);};
    //TF1 fsumsin("fsum", *funcSin.target<double(*)(double*, double*)>() , 0, 180);
    //TF1 fsum("fsum", *func.target<double(*)(double*, double*)>() , 0, 180);
    ::gsum = &gsum;

    TF1 fsumsin("fsumsim", funcSin , 0, 180);
    TF1 fsum("fsum", func , 0, 180);

    TH1D* (*cloneHisto)(const std::string&, TFile*) = [](const std::string& name, TFile* outFile){
        return  static_cast<TH1D*>(static_cast<TH1D*>(outFile->Get(name.c_str()))->Clone((name+"Corrected").c_str()));
    };

    TH1D* data = cloneHisto("dataCMNormalized", outFile);
    TH1D* simu = cloneHisto("simuCMNormalized", outFile);

    for (auto& it: std::vector<TH1D*>{data, simu}) {
        for (int i{0}; i < it->GetNbinsX(); ++i) {
            if (it->GetBinContent(i) == 0) continue;
            double theta = it->GetBinCenter(i);
            double thetaMin = it->GetBinLowEdge(i);
            double thetaMax = thetaMin + it->GetBinWidth(i);
            double factor = fsum.Eval(theta) * (cos(thetaMin) - cos(thetaMax)) / fsumsin.Integral(thetaMin, thetaMax);
            it->SetBinContent(i, it->GetBinContent(i) * factor);//NEED TO UPDATE ERRORS!!!
        }
        it->Write();
    }

    outFile->Write();
    outFile->Close();
}

void MaxLikelyhood(){
    gStyle->SetOptStat(0);

    std::map<std::string, double> xsections;
    std::map<std::string, std::string> files = SumThicknesses(xsections);
    files["data"]   = "./../../DataAnalyzed/sum.root";
    //files["data"]   = "./../../build/Out/sum.root";

    std::map<std::string, std::string> conditions;
    //conditions["data"] = "MugastData.M ==2 && MugastData.Z == 1 && VamosData.id_Z == 19 && VamosData.id_M == 47";
    //conditions["data"] = "VamosData.id_Z == 19 && VamosData.id_M == 47 && MugastData.T_proton<-2 && MugastData.T_proton>-13";
    conditions["data"] = "";
    conditions["s12"] = "";
    conditions["d32"] = "";
    conditions["f72"] = "";

    //Ex condition on data
    //conditions["data"] += "MugastData.Ex>-3 && MugastData.Ex<6.";

    //LEAVE THIS ALWAYS to avoid overflow
    //for (auto& it: conditions){
    //    std::string cond = "MugastData.Ex > -4 && MugastData.Ex <1 && MugastData.EmissionDirection.Theta()*180./TMath::Pi()>92 && MugastData.EmissionDirection.Theta()*180./TMath::Pi()<178";
    //    if (!it.second.empty() ) it.second += (" && "+cond);
    //    else it.second += cond;
    //}
    //for (auto& it: conditions){
    //    std::string cond = "MugastData.Ex > -3 && MugastData.Ex <3 && MugastData.EmissionDirection.Theta()*180./TMath::Pi()>92 && MugastData.EmissionDirection.Theta()*180./TMath::Pi()<178";
    //    if (!it.second.empty() ) it.second += (" && "+cond);
    //    else it.second += cond;
    //}
    //Add for 2D
    double Exmin{-5};
    double Exmax{4};
    double Thetamin{118.};
    double Thetamax{174};
    for (auto& it: conditions){
        std::string cond =  "MugastData.Ex > "+std::to_string(Exmin) +
                            " && MugastData.Ex < "+std::to_string(Exmax)+
                            " && MugastData.EmissionDirection.Theta()*180./TMath::Pi()> "+std::to_string(Thetamin)+ 
                            " && MugastData.EmissionDirection.Theta()*180./TMath::Pi()< "+std::to_string(Thetamax);
        if (!it.second.empty() ) it.second += (" && "+cond);
        else it.second += cond;
    }

    TFile* file = new TFile("./histofile.root", "read");
    if (!file->IsOpen()){
        TFile* outFile = new TFile("./histofile.root","recreate");
        for (const auto&it: files){
            TFile* file = new TFile(it.second.c_str(), "read");
            TTree* tree{nullptr};
            if(it.first == "data")
                tree = (TTree*) file->Get("AnalyzedTreeDeuterons");
            else
                tree = (TTree*) file->Get("AnalyzedTree");

            TH1D* histo{nullptr};
            TH2D* histo2D{nullptr};

            if(true){//constant width bins
                //histo parameters
                //int nbinsTheta{129};//GOOD LH
                //int nbinsTheta{113};
                //int nbinsTheta{100};
                int nbinsTheta{90};
                //int nbinsTheta{82};
                //int nbinsTheta{69};
                //int nbinsTheta{60};
                double startTheta{90};
                double stopTheta{180};
                int nbinsEx{30};
                double startEx{-10};
                double stopEx{10};

                //1 dimension
                histo = new TH1D(it.first.c_str(), it.first.c_str(),
                                 nbinsTheta, startTheta, stopTheta);

                //2 dimensions
                histo2D = new TH2D((it.first + "2D").c_str(), (it.first + "2D").c_str(),
                                   nbinsTheta, startTheta, stopTheta,
                                   nbinsEx, startEx, stopEx);
            }else{//Different width bins
                std::vector<double> binsTheta{106.5, 107.25, 108, 108.75, 109.5, 110.25, 111, 111.75, 112.5, 113.25, 114,
                                              114.75, 115.5, 116.25, 117, 117.75, 118.5, 119.25, 120, 120.75, 121.5,
                                              122.25, 123, 123.75, 124.5, 126, 127.5, 129,
                                              130.5, 132, 132.75, 133.5, 134.25, 135, 135.75, 136.5,
                                              137.25, 138, 138.75, 139.5, 140.25, 141, 141.75, 142.5, 143.25, 144,
                                              144.75, 146.25, 147.75, 149.25, 150, 150.75, 151.5,
                                              152.25, 153, 153.75, 154.5, 155.25, 156, 156.75, 157.5, 158.25, 159,
                                              159.75, 160.5, 161.25, 162, 162.75, 164.25, 165.75,
                                              167.25, 168.75, 170.25, 171.75, 173.25, 174.75};

                int nbinsEx{100};
                double startEx{-10};
                double stopEx{10};
                std::vector<double> binsEx;
                for (int i{0}; i<nbinsEx; ++i){
                    binsEx.push_back(startEx+(stopEx-startEx)/nbinsEx*i);
                }

                //1 dimension
                histo = new TH1D(it.first.c_str(), it.first.c_str(),
                                 binsTheta.size()-1, &binsTheta[0]);

                //2 dimensions
                histo2D = new TH2D((it.first + "2D").c_str(), (it.first + "2D").c_str(),
                                   binsTheta.size()-1, &binsTheta[0],
                                    binsEx.size()-1, &binsEx[0]);
            }

            tree->Draw(Form("MugastData.EmissionDirection.Theta()*180./TMath::Pi()>>%s", it.first.c_str()),
                       conditions[it.first].c_str());
            tree->Draw(Form("MugastData.Ex:MugastData.EmissionDirection.Theta()*180./TMath::Pi()>>%s",
                            (it.first + "2D").c_str()), conditions[it.first].c_str());


            TH1D* histoCM = new TH1D((it.first+"CM").c_str(), (it.first+"CM").c_str(), 90, 0, 90);
            if(it.first != "data")
                tree->Draw(Form("MugastData.Theta_CM*180./TMath::Pi()>>%s", (it.first+"CM").c_str()), conditions[it.first].c_str());
            else
                tree->Draw(Form("180 - MugastData.Theta_CM*180./TMath::Pi()>>%s", (it.first+"CM").c_str()), conditions[it.first].c_str());

            outFile->cd();
            histo->Write();
            histo2D->Write();
            histoCM->Write();
        }
        outFile->Write();
        outFile->Close();
        file = new TFile("./histofile.root", "read");
    }

    TH1D* data = (TH1D*) file->Get("data");
    std::vector<TH1D*> simu{(TH1D*) file->Get("s12"),
                            (TH1D*) file->Get("d32"),
                            (TH1D*) file->Get("f72")};
    TH2D* data2D = (TH2D*) file->Get("data2D");
    std::vector<TH2D*> simu2D{(TH2D*) file->Get("s122D"),
                            (TH2D*) file->Get("d322D"),
                            (TH2D*) file->Get("f722D")};

    LhResults results =  MaximizeLikelyhood(data, simu, data2D, simu2D, xsections);

    SumSimulationsWithPercentages(results, files);
    NormalizeHistograms("histofile.root", results);
    SaveTheoryDistributions("histofile.root", results);


}

void computeChiSquared(TH1D* data, TH1D* simu) {
    double xsq{0};
    int nonEmptyBins{0};
    for (int i = 0; i < data->GetNbinsX(); ++i) {
        double dt = data->GetBinContent(i);
        double dterr = data->GetBinError(i);
        double sm = simu->GetBinContent(i);
        if (dt>0 && sm>0){
            xsq += pow((dt-sm)/dterr,2);
            nonEmptyBins ++;
        }
    }
    double ndof = nonEmptyBins - 2;
    std::cout << "x squared : " << xsq << std::endl;
    std::cout << "DOF : " << ndof << std::endl;
    std::cout << "reduced x squared : " << xsq/ndof << std::endl;
    std::cout << "prob : " << TMath::Prob(xsq, ndof) << std::endl;

}

//using double for large numbers
long double factorial(int n, std::vector<long double>& mem){
    if (n == 1 || n == 0)
        return 1;

    if (mem.size()<2)
        mem.resize(2, 1);

    if (mem.size() == n)
        mem.push_back(mem[n-1]*n);

    if (mem.size() > n)
        return mem[n];

    return n*factorial(n-1, mem);
}
