#pragma once
#include <map>
#include <string>
#include <fstream>
#include <functional>

#include <TGraph.h>
#include <TF1.h>
#include <Math/ParamFunctor.h>

TGraph* distr{nullptr};

double AverageDistributions(std::map<std::string, double> files,
                            const std::string& outputFileName){

    //updates weights
    double sum{0};
    for (const auto& it: files){
        sum+=it.second;
    }
    for (auto& it: files){
        it.second/=sum;
    }

    std::map<double, double> values;
    std::map<TGraph*, double> graphs;
    for (const auto& it: files){
        graphs.emplace(new TGraph(it.first.c_str()), it.second);
    }

    auto graphReference = graphs.begin()->first;
    for(int i{0}; i<graphReference->GetN(); ++i){
        values[graphReference->GetPointX(i)] = 0;
    }

    for (auto& itValue: values){
        for (auto& itGraph: graphs){
            itValue.second += itGraph.first->Eval(itValue.first)*itGraph.second;
        }
    }

    std::ofstream outputFile(outputFileName);
    for (const auto& it: values){
        outputFile << it.first << "\t" << it.second << "\n";
    }
    outputFile.close();

    distr = new TGraph(outputFileName.c_str());

    TF1 integrand("integrand", [](double* x, double*p){return 2*TMath::Pi()*distr->Eval(x[0]*180./TMath::Pi())* sin(x[0]);}, 0, TMath::Pi());
    double integral = integrand.Integral(0, TMath::Pi());
    std::cout << "Integral of " << outputFileName << " has value: " << integral<< std::endl;
    return integral;
}