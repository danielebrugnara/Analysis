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
    double percentX{0};
    double percentY{0};
    int N;
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

LhResults MaximizeLikelyhood(TH1D* data, std::vector<TH1D*> simu){

    if (simu.size() != 3) return LhResults();
    int N{0};

    for (int j = 0; j < data->GetNbinsX(); ++j) {
        N += data->GetBinContent(j);
    }
    //std::vector<int> excludedBins {0, 12, 13, 14};
    //std::vector<int> excludedBins {6,7,8};
    //std::vector<int> excludedBins {0, 7, 8, 9};
    std::vector<int> excludedBins {0};

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

    results.N = inputs.N;

    return results;
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
    double sigmasRatioTotal = inputs.sigmasRatio[0]+inputs.sigmasRatio[1]+inputs.sigmasRatio[2];
    results.percentX = inputs.sigmasRatio[0]/sigmasRatioTotal * results.maxX;
    results.percentY = inputs.sigmasRatio[1]/sigmasRatioTotal * results.maxY;
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

std::map<std::string, std::string> SumThicknesses(){
    std::map<int, int> thicknessToNDeuterons{
            {20, 10},
            {25, 36},
            {30, 30},
            {35, 141},
            {40, 108},
            {45, 80},
            {50, 120},
            {55, 20}
    };
    int maxVal = std::max_element(thicknessToNDeuterons.begin(), thicknessToNDeuterons.end(),
                          [](const std::pair<int, int>& p1, const std::pair<int, int>& p2){
                                return p1.second < p2.second;
                            }
                        )->second;

    std::map<std::string, std::string> prefixes{
            {"s12", "selector_46Ar3Hed47K_0keV_s12"},
            {"d32", "selector_46Ar3Hed47K_360keV_d32"},
            {"f72", "selector_46Ar3Hed47K_2020keV_f72"},
            {"flat0", "selector_46Ar3Hed47K_0keV_flat"},
            {"flat2", "selector_46Ar3Hed47K_360keV_flat"},
            {"flat3", "selector_46Ar3Hed47K_2020keV_flat"}
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
            filePath += "_0_0_";
            filePath += std::to_string(itThickness.first);
            filePath += "um_analyzed.root";
            files[filePath] = itThickness.second;
        }
        outputFiles.emplace(itReaction.first, itReaction.second+"sumthickness.root");
        std::ifstream testFile(outputFiles[itReaction.first]);
        if (!testFile.is_open())
            CreateSimulationFile(files, outputFiles[itReaction.first]);
    }
    return outputFiles;
}

void SumSimulationsWithPercentages(const LhResults& results, std::map<std::string, std::string> files){
    {
        std::map<std::string, double> filesToSum;
        filesToSum.emplace(files["s12"], results.percentX);
        filesToSum.emplace(files["d32"], results.percentY);
        filesToSum.emplace(files["f72"], 1 - results.percentX - results.percentY);
        CreateSimulationFile(filesToSum, "simulationsum.root");
    }
    {
        std::map<std::string, double> filesToSum;
        filesToSum.emplace(files["flat0"], results.percentX);
        filesToSum.emplace(files["flat2"], results.percentY);
        filesToSum.emplace(files["flat3"], 1 - results.percentX - results.percentY);
        CreateSimulationFile(filesToSum, "simulationsum.root");
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

    TH1D normalizationTotal = results.percentX* *normalizedflat0CM + results.percentY* *normalizedflat2CM + (1-results.percentX - results.percentY)* *normalizedflat3CM;
    normalizationTotal.SetName("normalizationTotal");

    *dataCM = *dataCM / normalizationTotal;
    TH1D* simuCM = new TH1D((results.percentX* *normalizeds12CM + results.percentY* *normalizedd32CM + (1 - results.percentX - results.percentY)* *normalizedf72CM) / normalizationTotal);
    simuCM->SetName("simuCM");


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

void SaveTheoryDistributions(const std::string& fileName, const LhResults& results) {
    TFile *outFile = new TFile(fileName.c_str(), "update");
    if (!outFile->IsOpen()) throw std::runtime_error("histofile not present!!");
    std::cout << "Normalizing histograms\n";

    TGraph gs12("xsec_s12.txt");
    gs12.SetName("gs12");
    TGraph gd32("xsec_d32.txt");
    gd32.SetName("gd32");
    TGraph gf72("xsec_f72.txt");
    gf72.SetName("gf72");

    gs12.Write();
    gd32.Write();
    gf72.Write();

    TGraph gsum;
    gs12.SetName("gsum");

    for (int i{0}; i<gs12.GetN(); ++i){
        gsum.SetPoint(gsum.GetN(),
                      gs12.GetPointX(i),
                      gs12.GetPointY(i)* results.percentX
                        + gd32.GetPointY(i)* results.percentY
                        + gf72.GetPointY(i) *(1- results.percentX-results.percentY));
    }
    gsum.Write();
    outFile->Write();
    outFile->Close();
}

void MaxLikelyhood(){


    std::map<std::string, std::string> files = SumThicknesses();
    files["data"]   = "./../../DataAnalyzed/sum.root";

    std::map<std::string, std::string> conditions;
    conditions["data"] = "MugastData.M ==2 && MugastData.Z == 1 && VamosData.id_Z == 19 && VamosData.id_M == 47";
    conditions["s12"] = "";
    conditions["d32"] = "";
    conditions["f72"] = "";

    TFile* file = new TFile("./histofile.root", "read");
    if (!file->IsOpen()){
        TFile* outFile = new TFile("./histofile.root","recreate");
        for (const auto&it: files){
            TFile* file = new TFile(it.second.c_str(), "read");
            TTree* tree = (TTree*) file->Get("AnalyzedTree");
            TH1D* histo = new TH1D(it.first.c_str(), it.first.c_str(), 45, 90, 180);
            tree->Draw(Form("MugastData.Pos.Theta()*180./TMath::Pi()>>%s", it.first.c_str()), conditions[it.first].c_str());

            TH1D* histoCM = new TH1D((it.first+"CM").c_str(), (it.first+"CM").c_str(), 45, 0, 90);
            if(it.first != "data")
                tree->Draw(Form("MugastData.Theta_CM*180./TMath::Pi()>>%s", (it.first+"CM").c_str()), conditions[it.first].c_str());
            else
                tree->Draw(Form("180 - MugastData.Theta_CM*180./TMath::Pi()>>%s", (it.first+"CM").c_str()), conditions[it.first].c_str());

            outFile->cd();
            histo->Write();
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

    LhResults results =  MaximizeLikelyhood(data, simu);

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
