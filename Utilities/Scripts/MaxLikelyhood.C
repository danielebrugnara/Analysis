#include <vector>
#include <TH1D.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TH2D.h>
#include <TGraph.h>
#include <TGraph2D.h>
#include <TObjArray.h>
#include <stdio.h>
#include <TMath.h>
#include <math.h>

double binomial(int N, int k, double p);

void drawSigmas(TH2D* histo, double, double);

void computeChiSquared(TH1D* data, TH1D* simu);

struct LhResults{
    TH2D* histo{nullptr};
    TH1D* fit{nullptr};
    TH1D* dataWithErrors{nullptr};
    double maxVal{-1E10};
    double maxX{0};
    double maxY{0};
};

struct LhInputs{
    TH1D* data{nullptr};
    std::vector<TH1D*> simu;
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

void MaximizeLikelyhood(TH1D* data, std::vector<TH1D*> simu){

    if (simu.size() != 3) return;
    int N{0};

    for (int j = 0; j < data->GetNbinsX(); ++j) {
        N += data->GetBinContent(j);
    }
    std::vector<int> excludedBins {0, 12, 13, 14};

    //remove not wanted bins and normalize
    for(const auto&it: excludedBins){
        data->SetBinContent(it, 0);
    }
    for(int i=0; i<simu.size(); ++i){
        for(const auto&it: excludedBins) {
            simu[i]->SetBinContent(it, 0);//too many counts in empty bin??!!
        }
        double integral{0};
        for (int j = 0; j < simu[i]->GetNbinsX(); ++j) {
            integral += simu[i]->GetBinContent(j);
        }
        //std::cout << "before integral " << integral << std::endl;
        for (int j = 0; j < simu[i]->GetNbinsX(); ++j) {
            simu[i]->SetBinContent(j, simu[i]->GetBinContent(j)/integral);
        }
        integral = 0;
        for (int j = 0; j < simu[i]->GetNbinsX(); ++j) {
            integral += simu[i]->GetBinContent(j);
        }
    }

    LhInputs inputs(N);
    inputs.tsf.resize(2, nullptr);
    inputs.tsf[0] = [](double x, double y)->double { return x;};
    inputs.tsf[1] = [](double x, double y)->double { return y;};
    //inputs.tsf[2] = [](double x, double y)->double { return 1-x-y;};
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
    drawSigmas(results.histo, results.maxX, results.maxY);
    computeChiSquared(results.dataWithErrors, results.fit);


    LhInputs inputsPartial(N);
    inputsPartial.tsf.resize(2, nullptr);
    double centerX = results.maxY/(results.maxX+results.maxY);
    double centerY = (results.maxX+results.maxY);
    inputsPartial.startX = centerX*0.3;
    inputsPartial.endX = centerX*1.7;
    inputsPartial.startY = centerY*0.7;
    inputsPartial.endY = centerY*1.3;
    //inputsPartial.tsf[0] = [](double x, double y)->double { return x*(y-1);};
    //inputsPartial.tsf[1] = [](double x, double y)->double { return (x-1)*(y-1);};
    //inputsPartial.tsf[2] = [](double x, double y)->double { return y;};
    inputsPartial.tsf[0] = [](double x, double y)->double { return -y*(x-1);};//WHY is htis negative??!?!!?!?!??!!?!?!?
    inputsPartial.tsf[1] = [](double x, double y)->double { return x*y;};
    //inputsPartial.tsf[2] = [](double x, double y)->double { return 1-2*x-2*y;};
    inputsPartial.data = data;
    inputsPartial.simu = simu;

    LhResults resultsPartial = gridSearch(inputsPartial);
    std::cout << "xmin : " << resultsPartial.maxX << std::endl;
    std::cout << "ymin : " << resultsPartial.maxY << std::endl;
    std::cout << "max : " << resultsPartial.maxVal << std::endl;
    resultsPartial.histo->SetTitle("-Log LH; L=2 vs L=0+2; L=1+3");
    drawSigmas(resultsPartial.histo, resultsPartial.maxX, resultsPartial.maxY);
}

LhResults gridSearch(LhInputs& inputs) {
    LhResults results;
    TH1D* pVals = (TH1D*)inputs.data->Clone("pVals");
    results.histo = new TH2D("logLhPlot","- Log L;percent of L=0; percent of L=2", inputs.nX, inputs.startX, inputs.endX, inputs.nY, inputs.startY, inputs.endY);
    std::vector<double> coeffs;
    coeffs.resize(3, 0);
    for (int x = 0; x < inputs.nX; ++x) {
        double xval = inputs.startX + (inputs.endX - inputs.startX) / (inputs.nX - 1) * x;
        for (int y = 0; y < inputs.nY; ++y) {
            double yval = inputs.startY + (inputs.endY - inputs.startY) / (inputs.nY - 1) * y;
            coeffs[0] = inputs.tsf[0](xval, yval);
            coeffs[1] = inputs.tsf[1](xval, yval);
            //std::cout << "------------\n";
            //std::cout  << "coeff 0 " << coeffs[0] << std::endl;
            //std::cout  << "xval " << xval << std::endl;
            //std::cout  << "coeff 1 " << coeffs[1] << std::endl;
            //std::cout  << "yval " << yval << std::endl;
            //coeffs[2] = inputs.tsf[2](xval, yval);
            coeffs[2] = 1-coeffs[0]-coeffs[1];
            if (coeffs[0]<0 || coeffs[1] < 0 || coeffs[2] < 0) {
                results.histo->SetBinContent(results.histo->FindBin(xval, yval), 1);
                continue;
            }
            long double logLh{0};

            for (int j = 0; j < pVals->GetNbinsX(); ++j) {
                pVals->SetBinContent(j, 0);
            }

            for (int j = 0; j < pVals->GetNbinsX(); ++j) {
                if (isinf(logLh)) {
                    results.histo->SetBinContent(results.histo->FindBin(xval, yval), 1);
                    break;
                }
                for (int k = 0; k < coeffs.size(); ++k) {
                    if (isinf(logLh)) {
                        results.histo->SetBinContent(results.histo->FindBin(xval, yval), 1);
                        break;
                    }

                    pVals->SetBinContent(j, pVals->GetBinContent(j) + coeffs[k] * inputs.simu[k]->GetBinContent(j));
                    if (pVals->GetBinContent(j) > 1) {
                        std::cerr << "error in probability\n";
                    }
                }
            }

            for (int j = 0; j < pVals->GetNbinsX(); ++j) {
                logLh += logl((long double) binomial(inputs.N, inputs.data->GetBinContent(j), pVals->GetBinContent(j)));
            }
            results.histo->SetBinContent(results.histo->FindBin(xval, yval), -logLh);
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
    return results;
}

void drawSigmas(TH2D* histo, double minx, double miny){

    double val = histo->GetBinContent(histo->FindBin(minx, miny));

    std::vector<double> contours;
    for (int i=0; i<100; ++i){
        contours.push_back(val+(double)i);
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
    TH1D* data = (TH1D*) file->Get("his_data");
    std::vector<TH1D*> simu{(TH1D*) file->Get("his_602"),
                            (TH1D*) file->Get("his_618"),
                            (TH1D*) file->Get("his_633")};

    MaximizeLikelyhood(data, simu);
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