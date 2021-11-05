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

    outfile = new TFile(outputFileName.c_str(), "recreate");
    outfile->cd();
    tree = new TTree("AnalyzedTree", "AnalyzedTree");
    tree->Branch("MugastData", &mugastData);

    targetPos = TVector3(0, 0, 25.);

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

    if(     (*VamosAcceptedth == 0)
        ||  (*VamosAcceptedph == 0)
        ||  (*VamosAccepteddelta == 0)) return kTRUE;

    for (int i=0; i<Mugast->EventMultiplicity; ++i) {
        {
            Id id_tmp = {Mugast->DSSD_X[i], Mugast->DSSD_Y[i], Mugast->TelescopeNumber[i]};

            Pos pos_tmp = {Mugast->PosX[i], Mugast->PosY[i], Mugast->PosZ[i]};

            auto it = transform.find(id_tmp);

            if (it == transform.end()) {
                transform[id_tmp][pos_tmp] = 1;
            } else {
                auto it2 = it->second.find(pos_tmp);
                if (it2 == it->second.end()) {
                    it->second[pos_tmp] = 1;
                } else {
                    it2->second++;
                }
            }
        }
    }

    mugastData.~MugastData();
    new (&mugastData) MugastData(*Nev, 0);

    for (int i=0; i<*Nev; ++i){
        if (threasholdStripE.empty()) continue;
        if (threasholdStripT.empty()) continue;

        Id idE_tmp = {DSSD_X[i], -1, TelescopeNr[i]};
        Id idT_tmp = {-1, DSSD_Y[i], TelescopeNr[i]};
        auto itE = threasholdStripE.find(idE_tmp);
        auto itT = threasholdStripT.find(idT_tmp);
        if(itE != threasholdStripE.end() && itT != threasholdStripT.end() && EDep[i] > itE->second){
            mugastData.MG[i]      = TelescopeNr[i];
            mugastData.Pos[i]      = TVector3(X[i], Y[i], Z[i]);
            mugastData.EmissionDirection[i]  = mugastData.Pos[i] - targetPos;
            mugastData.SI_X[i]     = DSSD_X[i];
            mugastData.SI_Y[i]     = DSSD_Y[i];
            mugastData.SI_E[i]     = EDep[i];
            mugastData.E[i]        = ELab[i];
            mugastData.Ex[i]       = Ex[i];
            mugastData.Theta_CM[i] = ThetaCM[i];
        }else{
            //Initialize with -1
            mugastData.Pos[i]      = TVector3(0, 0, 0);
            mugastData.SI_X[i]     =-1000;
            mugastData.SI_Y[i]     =-1000;
            mugastData.SI_E[i]     =-1000;
            mugastData.E[i]        =-1000;
            mugastData.Ex[i]       =-1000;
            mugastData.Theta_CM[i] =-1000;
        }
    }

    //}

    tree->Fill();
    return kTRUE;
}

void Selector::SlaveTerminate()
{
    // The SlaveTerminate() function is called after all entries or objects
    // have been processed. When running with PROOF SlaveTerminate() is called
    // on each slave server.

    outfile->Write();
    outfile->Close();

}

void Selector::Terminate()
{
    // The Terminate() function is the last function to be called during
    // a query. It always runs on the client, it can be used to present
    // the results graphically or save the results to file.

}

void Selector::Init(TTree *tree)
{
    // The Init() function is called when the selector needs to initialize
    // a new tree or chain. Typically here the reader is initialized.
    // It is normally not necessary to make changes to the generated
    // code, but the routine can be extended by the user if needed.
    // Init() will be called many times when running on PROOF
    // (once per file to be processed).

    fReader.SetTree(tree);
}

Bool_t Selector::Notify()
{
    // The Notify() function is called when a new file is opened. This
    // can be either for a new TTree in a TChain or when when a new TTree
    // is started when using PROOF. It is normally not necessary to make changes
    // to the generated code, but the routine can be extended by the
    // user if needed. The return value is currently not used.

    return kTRUE;
}
