#define SelectorData_cxx
// The class definition in SelectorData.h has been generated automatically
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
// root> T->Process("SelectorData.C")
// root> T->Process("SelectorData.C","some options")
// root> T->Process("SelectorData.C+")
//


#include "SelectorData.h"
#include <TH2.h>
#include <TStyle.h>

void SelectorData::Begin(TTree * /*tree*/)
{
   // The Begin() function is called at the start of the query.
   // When running with PROOF Begin() is only called on the client.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();
}

void SelectorData::SlaveBegin(TTree * /*tree*/)
{
   // The SlaveBegin() function is called after the Begin() function.
   // When running with PROOF SlaveBegin() is called on each slave server.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();

}

Bool_t SelectorData::Process(Long64_t entry)
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

    for (int i=0; i<Mugast->multiplicity; ++i) {
        Id id_tmp = {Mugast->SI_X[i], Mugast->SI_Y[i], Mugast->MG[i]};

        Pos pos_tmp = {Mugast->Pos[i].x(), Mugast->Pos[i].y(), Mugast->Pos[i].z()};

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
        //Threasholds X
        {
            Id id_tmp = {Mugast->SI_X[i], -1, Mugast->MG[i]};
            auto itE = threasholdStripE.find(id_tmp);
            if (itE == threasholdStripE.end()) {
                threasholdStripE[id_tmp] = Mugast->SI_E[i];
            }else if(Mugast->SI_E[i] < itE->second){
                itE->second = Mugast->SI_E[i];
            }
            if (Mugast->M[i] != 0){
                threasholdStripT[id_tmp] = 1;
            }
        }
        //Threasholds Y
        {
            Id id_tmp = {-1, Mugast->SI_Y[i], Mugast->MG[i]};
            auto it = threasholdStripE.find(id_tmp);
            if (it == threasholdStripE.end()) {
                threasholdStripE[id_tmp] = Mugast->SI_E[i];
            }else if(Mugast->SI_E[i] < it->second){
                it->second = Mugast->SI_E[i];
            }
            if (Mugast->M[i] != 0){
                threasholdStripT[id_tmp] = 1;
            }
        }
    }
   return kTRUE;
}

void SelectorData::SlaveTerminate()
{
   // The SlaveTerminate() function is called after all entries or objects
   // have been processed. When running with PROOF SlaveTerminate() is called
   // on each slave server.

}

void SelectorData::Terminate()
{
   // The Terminate() function is the last function to be called during
   // a query. It always runs on the client, it can be used to present
   // the results graphically or save the results to file.

}