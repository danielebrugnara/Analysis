#ifndef EFFANALYSIS_PLOTTER_H
#define EFFANALYSIS_PLOTTER_H

#include "iostream"

#include "Globals.h"

#include "TROOT.h"
#include "TGClient.h"
#include "TCanvas.h"
#include "TPaveLabel.h"
#include "TVirtualPad.h"
#include "TGraph.h"
#include "TF1.h"
#include "TH1.h"

class Plotter {
public:
    Plotter()=default;


    template<typename T>
    void PlotOnCanvas(T& h, const std::string& opt) {
        gSystem->ProcessEvents();

        int ww = 1800;
        int wh = 1000;
        //Close canvas is not plotted on the same
        if (opt.find("same") == std::string::npos) {
            if (canvas != nullptr) {
                canvas->Close();
                gSystem->ProcessEvents();
                delete canvas;
                canvas = nullptr;
            }
            canvas = new TCanvas("plotter_canvas", "", ww, wh);
        } else {
            if (canvas == nullptr) {
                canvas = new TCanvas("plotter_canvas", "", ww, wh);
            }
        }
        UnZoom();

        gSystem->ProcessEvents();

        //Change line color if necessary
        TIter it(canvas->GetListOfPrimitives());
        TObject *obj = nullptr;
        std::vector<Color_t> colors;
        while ((obj = it.Next())) {
            auto *ptr = dynamic_cast<TAttLine *>(obj);
            if (ptr == nullptr)
                continue;
            colors.push_back(ptr->GetLineColor());
        }
        while (std::find(colors.begin(), colors.end(), h.GetLineColor()) != colors.end()) {
            h.SetLineColor(h.GetLineColor() + 10);
        }

        //Plotting histo on canvas
        canvas->cd();
        h.Draw(opt.c_str());
        //canvas->Modified();
        canvas->Update();
        gSystem->ProcessEvents();
        canvas->WaitPrimitive();

        gSystem->ProcessEvents();
    }
    static void WriteOnCanvas(const std::string& str, const Color_t& fillcolor=16){
        if (canvas == nullptr)
            return;

        UnZoom();

        canvas->cd();
        double x1 = canvas->GetPad(0)->GetX1();
        double x2 = canvas->GetPad(0)->GetX2();
        double y1 = canvas->GetPad(0)->GetY1();
        double y2 = canvas->GetPad(0)->GetY2();

        x1 = x1 + (3./4.)*(x2-x1);
        y1 = y1 + (3./4.)*(y2-y1);


        auto *title = new TPaveLabel(x1, y1, x2, y2,str.c_str());
        title->SetFillColor(fillcolor);
        title->SetTextFont(52);
        title->Draw();
        canvas->Update();

        gSystem->ProcessEvents();
        canvas->WaitPrimitive();

    }

    static void UnZoom(){

        if (canvas == nullptr)
            return;

        TObject* obj= nullptr;
        TIter it (canvas->GetListOfPrimitives());
        while ((obj = it.Next())) {
            if (dynamic_cast<TH1 *>(obj) != nullptr) {
                auto *ptr = dynamic_cast<TH1 *>(obj);
                ptr->GetXaxis()->UnZoom();
                ptr->GetYaxis()->UnZoom();
                continue;
            }
            if (dynamic_cast<TGraph *>(obj) != nullptr) {
                auto *ptr = dynamic_cast<TGraph *>(obj);
                ptr->GetXaxis()->UnZoom();
                ptr->GetYaxis()->UnZoom();
                continue;
            }
            if (dynamic_cast<TF1 *>(obj) != nullptr) {
                auto *ptr = dynamic_cast<TF1 *>(obj);
                ptr->GetXaxis()->UnZoom();
                ptr->GetYaxis()->UnZoom();
                continue;
            }

        }

    }

    static void SetRange(const double& valmin, const double& valmax){

        TObject* obj = nullptr;
        TIter it (canvas->GetListOfPrimitives());
        while ((obj = it.Next())) {
            if (dynamic_cast<TH1 *>(obj) != nullptr) {
                auto *ptr = dynamic_cast<TH1 *>(obj);
                ptr->GetXaxis()->SetRangeUser(valmin, valmax);
                continue;
            }
            if (dynamic_cast<TGraph *>(obj) != nullptr) {
                auto *ptr = dynamic_cast<TGraph *>(obj);
                ptr->GetXaxis()->SetRangeUser(valmin, valmax);
                continue;
            }
            if (dynamic_cast<TF1 *>(obj) != nullptr) {
                auto *ptr = dynamic_cast<TF1 *>(obj);
                ptr->GetXaxis()->SetRangeUser(valmin, valmax);
                continue;
            }
        }
        gSystem->ProcessEvents();
        canvas->WaitPrimitive();
    }
};


#endif //EFFANALYSIS_PLOTTER_H
