#include <vector>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>
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
#include <functional>

#include "CreateSimulationFile.C"
#include "AverageDistributions.C"

double binomial(int N, int k, double p);

void drawSigmas(TH2D* histo, double, double, TVirtualPad*);

void computeChiSquared(TH1D* data, TH1D* simu);

long double factorial(int n, std::vector<long double>& mem);

struct LhResults{
    TH2D* lhHisto{nullptr};
    TH2D* chi2Histo{nullptr};
    TH1D* fit{nullptr};
    std::vector<TH1D*> fitComponents;
    std::vector<std::pair<double, TH1D*>> fitSigmas;
    TH1D* discrepancy{nullptr};
    TH1D* dataWithErrors{nullptr};
    double maxVal{-1E10};
    double maxX{0};
    double maxY{0};
    int N;
};

struct LhInputs{
    TH1D* data{nullptr};
    std::vector<TH1D*> simu;
    int nX{500};
    std::vector<double> sigmas;
    int nY{500};
    double startX{0};
    double endX{1};
    double startY{0};
    double endY{1};
    //std::vector<double (*)(double, double )> tsf;
    std::vector<std::function<double(double, double)>> tsf;
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
    //std::vector<int> excludedBins {0, 11, 12, 13, 14};
    //std::vector<int> excludedBins {6,7,8};
    //std::vector<int> excludedBins {0};
    //std::vector<int> excludedBins {0, 24, 25, 26};
    std::vector<int> excludedBins {0, 40, 41, 42, 43, 44, 45, 46};
    //std::vector<int> excludedBins {0, 7, 8, 9};
    //std::vector<int> excludedBins {0, 13, 14, 15};
    //std::vector<int> excludedBins {0, 7, 8};

    //remove not wanted bins and normalize
    for(const auto&it: excludedBins){
        data->SetBinContent(it, 0);
    }
    for(int i=0; i<simu.size(); ++i){
        for(const auto&it: excludedBins) {
            simu[i]->SetBinContent(it, 0);//too many counts in empty bin??!!
        }
    }

    //Fist graph
    LhInputs inputs(N);
    inputs.tsf.resize(2, nullptr);
    inputs.tsf[0] = [](double x, double y)->double { return x;};
    inputs.tsf[1] = [](double x, double y)->double { return y;};
    inputs.data = data;
    inputs.simu = simu;
    inputs.sigmas ={0.016, 0.058, 0.105};
    LhResults results = gridSearch(inputs);
    std::cout << "xmin : " << results.maxX << std::endl;
    std::cout << "ymin : " << results.maxY << std::endl;
    std::cout << "max : " << results.maxVal << std::endl;

    auto* cv1 = new TCanvas();
    cv1->Divide(2,1);
    auto* tmpPad = cv1->cd(1);
    tmpPad->SetLogy();

    results.fit->SetLineColor(kRed);
    results.fit->SetMarkerColor(kRed);
    results.dataWithErrors->Draw("");
    results.fit->Draw("histo, same");

    results.fitComponents[0]->SetLineColor(kBlue);
    results.fitComponents[1]->SetLineColor(kGreen+1);
    results.fitComponents[2]->SetLineColor(kBlack);

    results.dataWithErrors->SetMarkerSize(3);
    results.fit->SetMarkerSize(3);
    results.fitComponents[0]->SetMarkerSize(3);
    results.fitComponents[1]->SetMarkerSize(3);
    results.fitComponents[2]->SetMarkerSize(3);
    results.dataWithErrors->SetMarkerSize(3);
    results.fit->SetMarkerSize(3);
    results.fitComponents[0]->SetMarkerSize(3);
    results.fitComponents[1]->SetMarkerSize(3);
    results.fitComponents[2]->SetMarkerSize(3);

    results.fitSigmas[2].second->SetLineStyle(8);

    results.fitComponents[0]->Draw("same, histo");
    results.fitComponents[1]->Draw("same, histo");
    results.fitComponents[2]->Draw("same, histo");
    //results.fitSigmas[2].second->Draw("same, histo");

    TLegend* leg = new TLegend();
    leg->AddEntry(results.dataWithErrors, "data", "lep");
    leg->AddEntry(results.fit, "fit", "l");
    leg->AddEntry(results.fitComponents[0], "L=0", "l");
    leg->AddEntry(results.fitComponents[1], "L=2", "l");
    leg->AddEntry(results.fitComponents[2], "L=3", "l");
    leg->Draw();


    cv1->cd(2);
    results.discrepancy->Draw();

    auto* cv2 = new TCanvas();
    cv2->Divide(2,1);
    TVirtualPad* pad21 = cv2->cd(1);

    results.lhHisto->SetName("lLH");
    drawSigmas(results.lhHisto, results.maxX, results.maxY, pad21);
    computeChiSquared(results.dataWithErrors, results.fit);

    std::cout << "Expected gammas @360 from L2: " << N*results.maxY*0.027  << std::endl;
    std::cout << "Expected gammas @360 from L3: " << N*(1-results.maxX-results.maxY)*0.013  << std::endl;
    std::cout << "Expected gammas @1660 from L3: " << N*(1-results.maxX-results.maxY)*0.0099  << std::endl;

    //2nd graph
    LhInputs inputsPartial(N);
    inputsPartial.tsf.resize(2, nullptr);
    double centerX = results.maxY/(results.maxX+results.maxY);
    double centerY = (results.maxX+results.maxY);
    inputsPartial.startX = 0.0;
    inputsPartial.endX = 0.30;
    inputsPartial.startY = centerY*0.7;
    inputsPartial.endY = centerY*1.3;

    //l2/l0+l2
    //inputsPartial.tsf[0] = [](double x, double y)->double { return -y*(x-1);};//WHY is htis negative??!?!!?!?!??!!?!?!?
    //inputsPartial.tsf[1] = [](double x, double y)->double { return x*y;};

    ////l2/l0
    inputsPartial.tsf[0] = [](double x, double y)->double { return y/(x+1);};    
    inputsPartial.tsf[1] = [](double x, double y)->double { return x*y/(x+1);};

    inputsPartial.data = data;
    inputsPartial.simu = simu;

    LhResults resultsPartial = gridSearch(inputsPartial);
    std::cout << "xmin : " << resultsPartial.maxX << std::endl;
    std::cout << "ymin : " << resultsPartial.maxY << std::endl;
    std::cout << "max : " << resultsPartial.maxVal << std::endl;
    resultsPartial.lhHisto->SetTitle("-Log L; L=2 / L=0; Percent of L=0+2");
    resultsPartial.lhHisto->SetName("partialLH");

    //cv2->cd(2);
    TVirtualPad* pad22 = cv2->cd(2);
    drawSigmas(resultsPartial.lhHisto, resultsPartial.maxX, resultsPartial.maxY, pad22);

    results.N = inputs.N;


    {
        //3rd graph
        LhInputs inputsPartial(N);
        inputsPartial.tsf.resize(2, nullptr);
        double centerX = results.maxY/(results.maxX+results.maxY);
        double centerY = (results.maxX+results.maxY);
        inputsPartial.startX = 0.0;
        inputsPartial.endX = 0.70;
        inputsPartial.startY = centerY*0.7;
        inputsPartial.endY = centerY*1.3;


        ////s2/s0
        std::vector<double> xsec{2.388, 1.534, 2.886 };
        std::vector<double> degeneracy{2, 4, 8};
        inputsPartial.tsf[0] = [xsec, degeneracy](double x, double y)->double { return y*(degeneracy[0]*xsec[0])/(x*degeneracy[1]*xsec[1]+degeneracy[0]*xsec[0]);};    
        inputsPartial.tsf[1] = [xsec, degeneracy](double x, double y)->double { return x*y*(degeneracy[1]*xsec[1])/(x*degeneracy[1]*xsec[1]+degeneracy[0]*xsec[0]);};

        inputsPartial.data = data;
        inputsPartial.simu = simu;

        LhResults resultsPartial = gridSearch(inputsPartial);
        TCanvas* cv = new TCanvas();
        drawSigmas(resultsPartial.lhHisto, resultsPartial.maxX, resultsPartial.maxY, cv);


        std::vector<double> sf {0.24, 0.4};
        TLine* line = new TLine(sf[1]/(sf[0]), inputsPartial.startY,sf[1]/(sf[0]), inputsPartial.endY);
        line->SetLineColor(kRed);
        line->SetLineWidth(3);
        cv->cd(1);
        line->Draw();
        cv->Update();

    }


    return results;
}

LhResults gridSearch(LhInputs& inputs) {
    std::vector<long double> factorialMem;
    LhResults results;
    TH1D* pVals = (TH1D*)inputs.data->Clone("pVals");
    results.lhHisto = new TH2D("logLhPlot","- Log L;percent of L=0; percent of L=2", inputs.nX, inputs.startX, inputs.endX, inputs.nY, inputs.startY, inputs.endY);
    std::vector<double> coeffs;
    coeffs.resize(3, 0);
    results.fitSigmas.resize(inputs.sigmas.size(), {0, nullptr});

    results.fitComponents.resize(3);
    results.fitComponents[0]= nullptr;
    results.fitComponents[1]= nullptr;
    results.fitComponents[2]= nullptr;

    for (int x = 0; x < inputs.nX; ++x) {
        double xval = inputs.startX + (inputs.endX - inputs.startX) / (inputs.nX - 1) * x;
        for (int y = 0; y < inputs.nY; ++y) {
            double yval = inputs.startY + (inputs.endY - inputs.startY) / (inputs.nY - 1) * y;
            //std::cout << "new pt\n";
            coeffs[0] = inputs.tsf[0](xval, yval);
            coeffs[1] = inputs.tsf[1](xval, yval);
            coeffs[2] = 1-coeffs[0]-coeffs[1];

            //std::cout   << "(p0 : " << coeffs[0] << " p2 " << coeffs[1] << ")->("
            //            << "x: " << xval << " y: " << yval << ")" << std::endl;

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
                for (unsigned long k = 0; k < coeffs.size(); ++k) {//Loop on L transfers
                    if (isinf(logLh)) {
                        results.lhHisto->SetBinContent(results.lhHisto->FindBin(xval, yval), 1);
                        break;
                    }

                    //Sigma ration changes normalization
                    pVals->SetBinContent(j, pVals->GetBinContent(j) + coeffs[k] * inputs.simu[k]->GetBinContent(j));
                }
                normalization +=pVals->GetBinContent(j);
            }
            if (normalization == 0 ) throw std::runtime_error("Null normalization!!!!\n");
            for (int j = 0; j < pVals->GetNbinsX(); ++j) {
                pVals->SetBinContent(j, 1./normalization*pVals->GetBinContent(j));
            }

            logLh += logl((long double) factorial(inputs.N, factorialMem));
            for (int j = 0; j < pVals->GetNbinsX(); ++j) {
                //logLh += logl((long double) binomial(inputs.N, inputs.data->GetBinContent(j), pVals->GetBinContent(j)));
                logLh += logl((long double) powl(pVals->GetBinContent(j),inputs.data->GetBinContent(j))/factorial((int)inputs.data->GetBinContent(j), factorialMem));
            }
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

                for(unsigned int s=0; s<inputs.sigmas.size(); ++s){
                    if (abs(yval - inputs.sigmas[s])<abs(yval - results.fitSigmas[s].first)){
                        results.fitSigmas[s].first = yval;
                        results.fitSigmas[s].second = (TH1D*) pVals->Clone(Form("%i_Sigma", s));
                        results.fitSigmas[s].second->Scale(inputs.N);
                    }
                }

            }
        }
    }
    results.discrepancy = (TH1D*) results.dataWithErrors->Clone("discrepancies");
    for (int i = 0; i < results.discrepancy->GetNbinsX(); ++i) {
        results.discrepancy->SetBinContent(i, results.dataWithErrors->GetBinContent(i)-results.fit->GetBinContent(i));
        results.discrepancy->SetBinError(i, results.dataWithErrors->GetBinError(i));
    }

    //double scale = results.fit->Integral(results.fit->FindBin(90), results.fit->FindBin(180));
    double scale = 2000;

    results.fitComponents[0] = static_cast<TH1D*>( inputs.simu[0]->Clone(Form("partialL%i", 0)));
    results.fitComponents[1] = static_cast<TH1D*>( inputs.simu[1]->Clone(Form("partialL%i", 2)));
    results.fitComponents[2] = static_cast<TH1D*>( inputs.simu[2]->Clone(Form("partialL%i", 3)));

    std::vector<int> normalization{0,0,0};
    for (int j = 0; j < pVals->GetNbinsX(); ++j) {
        normalization[0] += results.fitComponents[0]->GetBinContent(j);
        normalization[1] += results.fitComponents[1]->GetBinContent(j);
        normalization[2] += results.fitComponents[2]->GetBinContent(j);
    }

    results.fitComponents[0]->Scale(inputs.N*results.maxX/normalization[0]);
    results.fitComponents[1]->Scale(inputs.N*results.maxY/normalization[1]);
    results.fitComponents[2]->Scale(inputs.N*(1-results.maxX-results.maxY)/normalization[2]);
    return results;
}

void drawSigmas(TH2D* histo, double minx, double miny, TVirtualPad* cv){

    cv->Divide(1,3);
    cv->cd(1);

    double val = histo->GetBinContent(histo->FindBin(minx, miny));
    std::string name = histo->GetName();
    name += "_clone";
    auto* histoclone = (TH2D*)histo->Clone(name.c_str());

    std::vector<double> contours;
    for (int i=0; i<100; ++i){
        contours.push_back(val+(double)i*i/2.);
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

    int x, y, z;
    double nsigma = 3;
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
    cv->Update();

}

double binomial(int N, int k, double p){
    return TMath::Binomial(N, k) * pow(p, k)* pow(1-p, N-k);
}

std::map<std::string, std::string> SumThicknesses(){
    std::map<int, int> thicknessToNDeuterons{
        {15, 46},
            {20, 0},
            {25, 134},
            {30, 93},
            {35, 101},
            {40, 89},
            {45, 108},
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
        AverageDistributions(filesToAverage, itDistr+"averaged.txt");
    }


    return outputFiles;
}

void SumSimulationsWithPercentages(const LhResults& results, std::map<std::string, std::string> files){
    {
        std::map<std::string, double> filesToSum;
        filesToSum.emplace(files["s12"], results.maxX);
        filesToSum.emplace(files["d32"], results.maxY);
        filesToSum.emplace(files["f72"], 1 - results.maxX - results.maxY);
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
        filesToSum.emplace(files["flat0"], results.maxX);
        filesToSum.emplace(files["flat2"], results.maxY);
        filesToSum.emplace(files["flat3"], 1 - results.maxX - results.maxY);
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

    TH1D normalizationTotal = results.maxX* *normalizedflat0CM + results.maxY* *normalizedflat2CM + (1-results.maxX - results.maxY)* *normalizedflat3CM;
    normalizationTotal.SetName("normalizationTotal");

    *dataCM = *dataCM / normalizationTotal;
    TH1D* simuCM = new TH1D((results.maxX* *normalizeds12CM + results.maxY* *normalizedd32CM + (1 - results.maxX - results.maxY)* *normalizedf72CM) / normalizationTotal);
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
                gs12.GetPointY(i)* results.maxX
                + gd32.GetPointY(i)* results.maxY
                + gf72.GetPointY(i) *(1- results.maxX-results.maxY));
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


    std::map<std::string, std::string> files = SumThicknesses();
    files["data"]   = "./../../DataAnalyzed/sum.root";

    std::map<std::string, std::string> conditions;
    conditions["data"] = "MugastData.M ==2 && MugastData.Z == 1 && VamosData.id_Z == 19 && VamosData.id_M == 47 ";
    //conditions["data"] = "MugastData.M ==2 && MugastData.Z == 1 && VamosData.id_Z == 19 && VamosData.id_M == 47";
    conditions["s12"] = "";
    conditions["d32"] = "";
    conditions["f72"] = "";

    for (auto& it: conditions){
        if (!it.second.empty() ) it.second += " && MugastData.Ex < 2.5";
        else it.second += "MugastData.Ex < 2.5";
    }

    TFile* file = new TFile("./histofile.root", "read");
    if (!file->IsOpen()){
        TFile* outFile = new TFile("./histofile.root","recreate");
        for (const auto&it: files){
            TFile* file = new TFile(it.second.c_str(), "read");
            TTree* tree = (TTree*) file->Get("AnalyzedTree");
            //TH1D* histo = new TH1D(it.first.c_str(), it.first.c_str(), 90, 90, 180);
            TH1D* histo = new TH1D(it.first.c_str(), it.first.c_str(), 150, 90, 180);
            //TH1D* histo = new TH1D(it.first.c_str(), it.first.c_str(), 25, 90, 180);
            tree->Draw(Form("MugastData.EmissionDirection.Theta()*180./TMath::Pi()>>%s", it.first.c_str()), conditions[it.first].c_str());

            TH1D* histoCM = new TH1D((it.first+"CM").c_str(), (it.first+"CM").c_str(), 80, 0, 90);
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
