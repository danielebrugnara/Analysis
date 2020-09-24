#define Selector_cxx
// The class definition in Selector.h has been generated automatically
// by the ROOT utility TTree::MakeSelector(). This class is derived
// from the ROOT class TSelector. For more information on the TSelector
// framework see $ROOTSYS/README/README.SELECTOR or the ROOT User Manual.


// The following methods are defined in this file:
//    Begin():        called every time a loop on the tree starts,
//                    a convenient place to create your histograms.
//    SlaveBegin():   called after Begin(), when on PROOF called only on the
//                    slave servers.
//    Process():      called for each event, in this function you decide what
//                    to read and fill your histograms.
//    SlaveTerminate: called at the end of the loop on the tree, when on PROOF
//                    called only on the slave servers.
//    Terminate():    called at the end of the loop on the tree,
//                    a convenient place to draw/fit your histograms.
//
// To use this file, try the following session on your Tree T:
//
// root> T->Process("Selector.C")
// root> T->Process("Selector.C","some options")
// root> T->Process("Selector.C+")
//


#include "Selector.h"
#include <TH2.h>
#include <TStyle.h>

void Selector::Begin(TTree * /*tree*/)
{
    // The Begin() function is called at the start of the query.
    // When running with PROOF Begin() is only called on the client.
    // The tree argument is deprecated (on PROOF 0 is passed).

    TString option = GetOption();
}

void Selector::SlaveBegin(TTree * /*tree*/)
{
    // The SlaveBegin() function is called after the Begin() function.
    // When running with PROOF SlaveBegin() is called on each slave server.
    // The tree argument is deprecated (on PROOF 0 is passed).

    TString option = GetOption();

    data_addb_spec = new TH1D("data_addb_spec", "data_addb_spec",
                              4000, 0, 4000);
    fOutput->Add(data_addb_spec);
    data_addb_gg = new TH2D("data_addb_gg", "data_addb_gg",
                            4000, 0, 4000,
                            4000, 0, 4000);
    fOutput->Add(data_addb_gg);

    data_core_spec = new TH1D("data_core_spec", "data_core_spec",
                              4000, 0, 4000);
    fOutput->Add(data_core_spec);
    data_core_gg = new TH2D("data_core_gg", "data_core_gg",
                            4000, 0, 4000,
                            4000, 0, 4000);
    fOutput->Add(data_core_gg);

    mspec_core = new TH2D("mspec_core", "mspec_core",
                          50, 0, 50,
                          4000, 0, 4000);
    fOutput->Add(mspec_core);

    geom   = new TH3D("geom", "geom",
                      100, -350, 350,
                      100, -350, 350,
                      100, -350, 350);
    fOutput->Add(geom);


    //Individual crystal core spectra
    for(int i=0; i<=45; ++i){
        crystal_spectra.push_back(new TH1D(Form("crtstal_spectra_%i", i),
                                           Form("crtstal_spectra_%i", i),
                                           4000, 0, 4000));
        fOutput->Add(crystal_spectra.back());;
    }
}

Bool_t Selector::Process(Long64_t entry)
{
    // The Process() function is called for each entry in the tree (or possibly
    // keyed object in the case of PROOF) to be processed. The entry argument
    // specifies which entry in the currently loaded tree is to be processed.
    // When processing keyed objects with PROOF, the object is already loaded
    // and is available via the fObject pointer.
    //
    // This function should contain the \"body\" of the analysis. It can contain
    // simple or elaborate selection criteria, run algorithms on the data
    // of the event and typically fill histograms.
    //
    // The processing can be stopped by calling Abort().
    //
    // Use fStatus to set the return value of TTree::Process().
    //
    // The return value is currently not used.

    fReader.SetLocalEntry(entry);

    std::unordered_map<int, double> coreEn;
    std::unordered_map<int, double> coreX;
    std::unordered_map<int, double> coreY;
    std::unordered_map<int, double> coreZ;
    for (unsigned int i=0; i<hitE.GetSize(); ++i){
        int id = hitId.At(i);
        if (coreEn.find(id) == coreEn.end()){
            coreEn.emplace(id, hitE.At(i));
            coreX.emplace(id, hitGX.At(i));
            coreY.emplace(id, hitGY.At(i));
            coreZ.emplace(id, hitGZ.At(i));
        }else{
            if (coreEn.at(id)<hitE.At(i)){
                coreX.at(id) = hitGX.At(i);
                coreY.at(id) = hitGY.At(i);
                coreZ.at(id) = hitGZ.At(i);
            }
            coreEn.at(id) += hitE.At(i);
        }
    }

    for (const auto& it: coreEn){
        geom->Fill(coreX.at(it.first),coreY.at(it.first), coreZ.at(it.first));
    }

    for (const auto & en1: AddE){
        data_addb_spec->Fill(en1);
        for (const auto & en2: AddE){
            if (en1 == en2) continue;
            data_addb_gg->Fill(en1, en2);
            data_addb_gg->Fill(en2, en1);
        }
    }
    for (const auto & en1: coreE0){
        data_core_spec->Fill(en1);
        for (const auto & en2: AddE){
            if (en1 == en2) continue;
            data_core_gg->Fill(en1, en2);
            data_core_gg->Fill(en2, en1);
        }
    }

    for (unsigned long ii=0; ii<coreE0.GetSize(); ++ii) {
        mspec_core->Fill(coreId[ii], coreE0[ii]);
        crystal_spectra[coreId[ii]]->Fill(coreE0[ii]);
    }
    return kTRUE;
}

void Selector::SlaveTerminate()
{
    // The SlaveTerminate() function is called after all entries or objects
    // have been processed. When running with PROOF SlaveTerminate() is called
    // on each slave server.

}

void Selector::Terminate()
{
    // The Terminate() function is the last function to be called during
    // a query. It always runs on the client, it can be used to present
    // the results graphically or save the results to file.

    std::vector<TF1> threasholds;
    bool show_canvas = false;
    for (unsigned int i=0; i<crystal_spectra.size(); ++i){
        RooRealVar x("x", "x", 0, 150);
        x.setRange("fullrange", 0, 150);
        x.setRange("partial",1, 149);

        RooDataHist dh("dh", "dh", x, RooFit::Import(*(crystal_spectra[i])));

        //RooRealVar scale("scale", "scale", 1000/2., 1, 1000/2.*10);
        RooRealVar threashold("threashold", "threashold", 80., 10, 180);
        RooRealVar lambda("lambda", "lambda", 20., 1, 60.);
        RooGenericPdf model(Form("model_%i", i),
                                Form("model_%i", i),
                                "(1+erf((x-threashold)/lambda))",
                                RooArgSet(x, threashold, lambda));

        //Plotting
        RooPlot* frame = x.frame(RooFit::Title(Form("Threashold_%i", i)));
        dh.plotOn(frame);
        auto *result = model.fitTo(dh, RooFit::Save());
        model.plotOn(frame);

        RooArgList pars(*model.getParameters(RooArgSet(x)));
        RooArgSet prodset(model);
        RooProduct normPdf(Form("normmodel_%i", i),
                           Form("normmodel_%i", i),
                           prodset);

        threasholds.push_back(*normPdf.asTF(RooArgList(x), pars));
        threasholds.back().SaveAs(Form("tmp_threa_%i.root", i));

        if (show_canvas) {
            auto* cv = new TCanvas();
            frame->Draw();
            cv->WaitPrimitive();
            delete cv;
        }
    }
    //Note that write does not behave correctly, so saveas is needed
    system("hadd threasholds.root tmp_threa_*.root");
    system("rm tmp_threa_*.root");
}
