//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Apr 13 18:59:38 2021 by ROOT version 6.20/04
// from TTree AnalyzedTree/AnalyzedTree
// found on file: ../../../DataAnalyzed/sum.root
//////////////////////////////////////////////////////////

#ifndef SelectorData_h
#define SelectorData_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>

#include "MugastData.h"
#include "Selector.h"

// Headers needed by this particular selector


class SelectorData : public TSelector {
public :
    TTreeReader fReader;  //!the tree reader
    TTree *fChain = 0;   //!pointer to the analyzed TTree or TChain

    // Readers to access the data (delete the ones you do not need).
    //TTreeReaderValue<VamosData> VamosData = {fReader, "VamosData"};
    TTreeReaderValue<MugastData> Mugast = {fReader, "MugastData"};
    //TTreeReaderValue<Must2Data> Must2Data = {fReader, "Must2Data"};
    //TTreeReaderValue<CatsData> CatsData = {fReader, "CatsData"};
    //TTreeReaderValue<AgataData> AgataData = {fReader, "AgataData"};


    SelectorData(TTree * /*tree*/ = 0) {}
   virtual ~SelectorData() { }
   virtual Int_t   Version() const { return 2; }
   virtual void    Begin(TTree *tree);
   virtual void    SlaveBegin(TTree *tree);
   virtual void    Init(TTree *tree);
   virtual Bool_t  Notify();
   virtual Bool_t  Process(Long64_t entry);
   virtual Int_t   GetEntry(Long64_t entry, Int_t getall = 0) { return fChain ? fChain->GetTree()->GetEntry(entry, getall) : 0; }
   virtual void    SetOption(const char *option) { fOption = option; }
   virtual void    SetObject(TObject *obj) { fObject = obj; }
   virtual void    SetInputList(TList *input) { fInput = input; }
   virtual TList  *GetOutputList() const { return fOutput; }
   virtual void    SlaveTerminate();
   virtual void    Terminate();

   typedef Selector::Pos Pos;
   typedef Selector::Id Id;
   typedef Selector::Transform Transform;

   Transform transform;
    std::map<Id, double> threasholdStripE;
    std::map<Id, double> threasholdStripT;

   Transform GetTsf(){return transform;};
   std::map<Id, double> GetThsE(){return threasholdStripE;};
   std::map<Id, double> GetThsT(){return threasholdStripT;};


   ClassDef(SelectorData,0);

};

#endif

#ifdef SelectorData_cxx
void SelectorData::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the reader is initialized.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   fReader.SetTree(tree);
}

Bool_t SelectorData::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}


#endif // #ifdef SelectorData_cxx
