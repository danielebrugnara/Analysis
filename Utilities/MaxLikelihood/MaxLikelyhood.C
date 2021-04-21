#include <vector>
#include <TH1D.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TH2D.h>
#include <TH1D.h>
#include <TGraph.h>
#include <TGraph2D.h>
#include <TObjArray.h>
#include <stdio.h>
#include <TMath.h>
#include <TLine.h>
#include <math.h>

#include "CreateSimulationFile.C"

double binomial(int N, int k, double p);

void drawSigmas(TH2D* histo, double, double);

void computeChiSquared(TH1D* data, TH1D* simu);

long double factorial(int n, std::vector<long double>& mem);

struct LhResults{
    TH2D* lhHisto{nullptr};
    TH2D* chi2Histo{nullptr};
    TH1D* fit{nullptr};
    TH1D* discrepancy{nullptr};
    TH1D* dataWithErrors{nullptr};
    double maxVal{-1E10};
    double maxX{0};
    double maxY{0};
};

struct LhInputs{
    TH1D* data{nullptr};
    std::vector<TH1D*> simu;
    std::vector<double> sigmasRatio;
    int nX{500};
    int nY{500};
    double startX{0};
    double endX{1};
    double startY{0};
    double endY{1};
    std::vector<double (*)(double, double )> tsf;
    int N;
    LhInputs(const int& N): N(N){};
};

LhResults gridSearch(LhInputs& inputs);

std::vector<double> MaximizeLikelyhood(TH1D* data, std::vector<TH1D*> simu){

    if (simu.size() != 3) return std::vector<double>();
    int N{0};

    for (int j = 0; j < data->GetNbinsX(); ++j) {
        N += data->GetBinContent(j);
    }
    //std::vector<int> excludedBins {0, 12, 13, 14};
    //std::vector<int> excludedBins {6,7,8};
    std::vector<int> excludedBins {0, 7, 8, 9};

    //remove not wanted bins and normalize
    for(const auto&it: excludedBins){
        data->SetBinContent(it, 0);
    }
    for(int i=0; i<simu.size(); ++i){
        for(const auto&it: excludedBins) {
            simu[i]->SetBinContent(it, 0);//too many counts in empty bin??!!
        }
    }

    std::vector<double> sigmasRatio{1./2.48, 1./2.65, 1./3.2};

    LhInputs inputs(N);
    inputs.tsf.resize(2, nullptr);
    inputs.tsf[0] = [](double x, double y)->double { return x;};
    inputs.tsf[1] = [](double x, double y)->double { return y;};
    inputs.sigmasRatio = sigmasRatio;
    inputs.data = data;
    inputs.simu = simu;
    LhResults results = gridSearch(inputs);
    std::cout << "xmin : " << results.maxX << std::endl;
    std::cout << "ymin : " << results.maxY << std::endl;
    std::cout << "max : " << results.maxVal << std::endl;

    new TCanvas();
    results.fit->SetLineColor(kRed);
    results.fit->SetMarkerColor(kRed);
    results.dataWithErrors->Draw("");
    results.fit->Draw("histo, same");
    drawSigmas(results.lhHisto, results.maxX, results.maxY);
    computeChiSquared(results.dataWithErrors, results.fit);

    auto* tmpCv = new TCanvas();
    results.discrepancy->Draw();

    std::cout << "Expected gammas @360 from L2: " << N*results.maxY*0.027  << std::endl;
    std::cout << "Expected gammas @360 from L3: " << N*(1-results.maxX-results.maxY)*0.013  << std::endl;
    std::cout << "Expected gammas @1660 from L3: " << N*(1-results.maxX-results.maxY)*0.0099  << std::endl;

    LhInputs inputsPartial(N);
    inputsPartial.tsf.resize(2, nullptr);
    double centerX = results.maxY/(results.maxX+results.maxY);
    double centerY = (results.maxX+results.maxY);
    inputsPartial.startX = centerX*0.01;
    inputsPartial.endX = centerX*4;
    inputsPartial.startY = centerY*0.7;
    inputsPartial.endY = centerY*1.3;
    inputsPartial.sigmasRatio = sigmasRatio;
    inputsPartial.tsf[0] = [](double x, double y)->double { return -y*(x-1);};//WHY is htis negative??!?!!?!?!??!!?!?!?
    inputsPartial.tsf[1] = [](double x, double y)->double { return x*y;};
    inputsPartial.data = data;
    inputsPartial.simu = simu;

    LhResults resultsPartial = gridSearch(inputsPartial);
    std::cout << "xmin : " << resultsPartial.maxX << std::endl;
    std::cout << "ymin : " << resultsPartial.maxY << std::endl;
    std::cout << "max : " << resultsPartial.maxVal << std::endl;
    resultsPartial.lhHisto->SetTitle("-Log LH;Percent of L=2 vs L=0+2; Percent of L=1+2");
    drawSigmas(resultsPartial.lhHisto, resultsPartial.maxX, resultsPartial.maxY);
    std::vector<double> sf {0.24, 0.4};

    TLine* theoricValue = new TLine(sf[1]/(sf[1]+sf[0]), inputsPartial.startY,sf[1]/(sf[1]+sf[0]), inputsPartial.endY);
    theoricValue->Draw();

    return std::vector<double>{results.maxX, results.maxY, 1-results.maxX-results.maxY};
}

LhResults gridSearch(LhInputs& inputs) {
    std::vector<long double> factorialMem;
    LhResults results;
    TH1D* pVals = (TH1D*)inputs.data->Clone("pVals");
    results.lhHisto = new TH2D("logLhPlot","- Log L;percent of L=0; percent of L=2", inputs.nX, inputs.startX, inputs.endX, inputs.nY, inputs.startY, inputs.endY);
    std::vector<double> coeffs;
    coeffs.resize(3, 0);
    for (int x = 0; x < inputs.nX; ++x) {
        double xval = inputs.startX + (inputs.endX - inputs.startX) / (inputs.nX - 1) * x;
        for (int y = 0; y < inputs.nY; ++y) {
            double yval = inputs.startY + (inputs.endY - inputs.startY) / (inputs.nY - 1) * y;
            coeffs[0] = inputs.tsf[0](xval, yval);
            coeffs[1] = inputs.tsf[1](xval, yval);
            coeffs[2] = 1-coeffs[0]-coeffs[1];
            if (coeffs[0]<0 || coeffs[1] < 0 || coeffs[2] < 0) {
                results.lhHisto->SetBinContent(results.lhHisto->FindBin(xval, yval), 1);
                continue;
            }
            long double logLh{0};
            long double chi2{0};

            for (int j = 0; j < pVals->GetNbinsX(); ++j) {
                pVals->SetBinContent(j, 0);
            }

            double normalization{0};
            for (int j = 0; j < pVals->GetNbinsX(); ++j) {
                for (int k = 0; k < coeffs.size(); ++k) {
                    if (isinf(logLh)) {
                        results.lhHisto->SetBinContent(results.lhHisto->FindBin(xval, yval), 1);
                        break;
                    }

                    //Sigma ration changes normalization
                    pVals->SetBinContent(j, pVals->GetBinContent(j) + inputs.sigmasRatio[k] * coeffs[k] * inputs.simu[k]->GetBinContent(j));
                }
                normalization +=pVals->GetBinContent(j);
            }
            for (int j = 0; j < pVals->GetNbinsX(); ++j) {
                pVals->SetBinContent(j, 1./normalization*pVals->GetBinContent(j));
            }

            logLh += logl((long double) factorial(inputs.N, factorialMem));
            //std::cout <<  "loglh : " <<logLh << std::endl;
            for (int j = 0; j < pVals->GetNbinsX(); ++j) {
                //logLh += logl((long double) binomial(inputs.N, inputs.data->GetBinContent(j), pVals->GetBinContent(j)));
                logLh += logl((long double) powl(pVals->GetBinContent(j),inputs.data->GetBinContent(j))/factorial((int)inputs.data->GetBinContent(j), factorialMem));
            }
            //std::cout <<  "loglh : " <<logLh << std::endl;
            results.lhHisto->SetBinContent(results.lhHisto->FindBin(xval, yval), -logLh);
            if (logLh > results.maxVal) {
                results.maxVal = logLh;
                results.maxX = xval;
                results.maxY = yval;
                if (results.fit != nullptr){
                    results.fit->Delete();
                }
                results.fit = (TH1D *) pVals->Clone("fitResult");
                results.fit->Scale(inputs.N);

                if (results.dataWithErrors != nullptr){
                    results.dataWithErrors->Delete("");
                }
                results.dataWithErrors = (TH1D*) inputs.data->Clone("DataWithErrors");
                for (int j = 0; j < pVals->GetNbinsX(); ++j) {
                    results.dataWithErrors->SetBinError(j, sqrt(inputs.N * pVals->GetBinContent(j) * (1 - pVals->GetBinContent(j))));
                }
            }
        }
    }
    results.discrepancy = (TH1D*) results.dataWithErrors->Clone("discrepancies");
    for (int i = 0; i < results.discrepancy->GetNbinsX(); ++i) {
        results.discrepancy->SetBinContent(i, results.dataWithErrors->GetBinContent(i)-results.fit->GetBinContent(i));
        results.discrepancy->SetBinError(i, results.dataWithErrors->GetBinError(i));
    }

        return results;
}

void drawSigmas(TH2D* histo, double minx, double miny){

    double val = histo->GetBinContent(histo->FindBin(minx, miny));

    std::vector<double> contours;
    for (int i=0; i<100; ++i){
        contours.push_back(val+(double)i*i/2.);
    }
    histo->SetContour(contours.size(), &contours[0]);

    auto* cv = new TCanvas();
    //histo->Draw("colz");
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


    histo->Draw("COLZ");
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

    cv->Update();
}

double binomial(int N, int k, double p){
    return TMath::Binomial(N, k) * pow(p, k)* pow(1-p, N-k);
}

void MaxLikelyhood(){

    TFile* file = new TFile("./histofile.root", "read");
    std::map<std::string, std::string> files;
    std::map<std::string, std::string> conditions;
    files["data"]   = "./../../DataAnalyzed/sum.root";
    conditions["data"] = "MugastData.M ==2 && MugastData.Z == 1 && VamosData.id_Z == 19 && VamosData.id_M == 47";
    files["s12"]    = "./../../DataAnalyzed/simu/46Ar3Hed47K_0keV_s12_0_0_analyzed.root";
    conditions["s12"] = "";
    files["d32"]    = "./../../DataAnalyzed/simu/46Ar3Hed47K_360keV_d32_0_0_analyzed.root";
    conditions["d32"] = "";
    files["f72"]    = "./../../DataAnalyzed/simu/46Ar3Hed47K_2020keV_f72_0_0_analyzed.root";
    conditions["f72"] = "";
    if (!file->IsOpen()){
        TFile* outFile = new TFile("./histofile.root","recreate");
        for (const auto&it: files){
            TFile* file = new TFile(it.second.c_str(), "read");
            TTree* tree = (TTree*) file->Get("AnalyzedTree");
            TH1D* histo = new TH1D(it.first.c_str(), it.first.c_str(), 45, 90, 180);
            tree->Draw(Form("MugastData.Pos.Theta()*180./TMath::Pi()>>%s", it.first.c_str()), conditions[it.first].c_str());

            outFile->cd();
            histo->Write();
        }
        outFile->Write();
        outFile->Close();
        file = new TFile("./histofile.root", "read");
    }

    TH1D* data = (TH1D*) file->Get("data");
    std::vector<TH1D*> simu{(TH1D*) file->Get("s12"),
                            (TH1D*) file->Get("d32"),
                            (TH1D*) file->Get("f72")};

    std::vector<double>results =  MaximizeLikelyhood(data, simu);

    bool createRootFile{true};
    if (createRootFile){
        std::vector<std::pair<std::string, double>> filesToSum;
        filesToSum.emplace_back(files["s12"], results[0]);
        filesToSum.emplace_back(files["d32"], results[1]);
        filesToSum.emplace_back(files["f72"], results[2]);
        CreateSimulationFile(filesToSum);
    }
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
