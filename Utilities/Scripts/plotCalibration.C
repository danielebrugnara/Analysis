#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TH2D.h>

//The script must be run in the agata cloud due to root file compatibility issues
void plotCalibration(){
    gStyle->SetOptStat(kFALSE);
    gStyle->SetPalette(kLightTemperature);
    gStyle->SetCanvasPreferGL(kTRUE);


    TFile* f1 = new TFile("run_0501.root", "read");
    TFile* f2 = new TFile("testdaniele.root", "read");

    auto* t1 = (TTree*)f1->Get("AutoTree");
    auto* t2 = (TTree*)f2->Get("PhysicsTree");

    auto* h1 = new TH2D("h1", ";Strip Number; Channel", 128, 0, 128, 500, 8200, 9000);
    auto* h2 = new TH2D("h2", ";Strip Number; Energy [MeV]", 128, 0, 128, 500, 0, 6);

    h1->GetXaxis()->SetLabelSize(0.05);
    h1->GetXaxis()->SetTitleSize(0.05);
    h1->GetYaxis()->SetLabelSize(0.05);
    h1->GetYaxis()->SetNdivisions(503);
    h1->GetYaxis()->SetTitleSize(0.05);
    h1->GetYaxis()->SetTitleOffset(0.9);
    h2->GetXaxis()->SetLabelSize(0.05);
    h2->GetXaxis()->SetTitleSize(0.05);
    h2->GetYaxis()->SetLabelSize(0.05);
    h2->GetYaxis()->SetTitleSize(0.05);
    h2->GetYaxis()->SetTitleOffset(0.9);

    auto* cv1 = new TCanvas();
    cv1->cd();
    t1->Draw("Mugast.fMG_DSSDXE_Energy:Mugast.fMG_DSSDXE_StripNbr>>h1", "Mugast.fMG_DSSDXE_DetectorNbr==1", "col");
    cv1->SaveAs("pre_energy_calibration.C");
    cv1->SaveAs("pre_energy_calibration.pdf");
    auto* cv2 = new TCanvas();
    cv2->cd();
    t2->Draw("Mugast.DSSD_E:Mugast.DSSD_X>>h2", "Mugast.TelescopeNumber==1", "col");
    cv2->SaveAs("post_energy_calibration.C");
    cv2->SaveAs("post_energy_calibration.pdf");

}
