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
   
    std::unordered_map<int,Hit> hits_analyzed;
    double cr;
   for (const auto & h: Hits){
       if (hits_analyzed.find(cr = h.GetCrystal())== hits_analyzed.end()){
            hits_analyzed[cr] = h; 
        }else{
//            hits_analyzed[cr] = Hit(cr, 
//                                    h.GetEnergy()+hits_analyzed[cr].GetEnergy(),
//                                    (h.GetX()*h.GetEnergy()+hits_analyzed[cr].GetX()*hits_analyzed[cr].GetEnergy())/(h.GetEnergy()+hits_analyzed[cr].GetEnergy()),
//                                    (h.GetY()*h.GetEnergy()+hits_analyzed[cr].GetY()*hits_analyzed[cr].GetEnergy())/(h.GetEnergy()+hits_analyzed[cr].GetEnergy()),
//                                    (h.GetZ()*h.GetEnergy()+hits_analyzed[cr].GetZ()*hits_analyzed[cr].GetEnergy())/(h.GetEnergy()+hits_analyzed[cr].GetEnergy()),
//                                    h.GetDetector(),
//                                    h.GetEnergy()>hits_analyzed[cr].GetEnergy() ? h.GetSegment() : hits_analyzed[cr].GetSegment());
        }
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

}
