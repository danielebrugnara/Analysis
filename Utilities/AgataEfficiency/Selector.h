//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Jul 10 16:26:18 2020 by ROOT version 6.20/04
// from TTree TreeMaster/TreeMaster
// found on file: ../../Data/AgataStandalone/Tree_0000.root
//////////////////////////////////////////////////////////

#ifndef Selector_h
#define Selector_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>

#include <TH1D.h>
#include <TH2D.h>

// Headers needed by this particular selector


class Selector : public TSelector {
public :
   TTreeReader     fReader;  //!the tree reader
   TTree          *fChain = 0;   //!pointer to the analyzed TTree or TChain

   // Readers to access the data (delete the ones you do not need).
   TTreeReaderValue<Int_t> nbHits = {fReader, "nbHits"};
   TTreeReaderArray<Float_t> hitE = {fReader, "hitE"};
   TTreeReaderArray<Float_t> hitX = {fReader, "hitX"};
   TTreeReaderArray<Float_t> hitY = {fReader, "hitY"};
   TTreeReaderArray<Float_t> hitZ = {fReader, "hitZ"};
   TTreeReaderArray<Float_t> hitGX = {fReader, "hitGX"};
   TTreeReaderArray<Float_t> hitGY = {fReader, "hitGY"};
   TTreeReaderArray<Float_t> hitGZ = {fReader, "hitGZ"};
   TTreeReaderArray<Int_t> hitId = {fReader, "hitId"};
   TTreeReaderArray<Int_t> hitSg = {fReader, "hitSg"};
   TTreeReaderValue<Int_t> nbCores = {fReader, "nbCores"};
   TTreeReaderArray<Int_t> nbHitsperCry = {fReader, "nbHitsperCry"};
   TTreeReaderArray<Int_t> coreId = {fReader, "coreId"};
   TTreeReaderArray<Float_t> coreE0 = {fReader, "coreE0"};
   TTreeReaderArray<Float_t> coreE1 = {fReader, "coreE1"};
   TTreeReaderArray<Float_t> coreT0 = {fReader, "coreT0"};
   TTreeReaderArray<ULong64_t> coreTS = {fReader, "coreTS"};
   TTreeReaderValue<Int_t> nbAdd = {fReader, "nbAdd"};
   TTreeReaderArray<Int_t> AddId = {fReader, "AddId"};
   TTreeReaderArray<Float_t> AddE = {fReader, "AddE"};
   TTreeReaderArray<Float_t> AddX = {fReader, "AddX"};
   TTreeReaderArray<Float_t> AddY = {fReader, "AddY"};
   TTreeReaderArray<Float_t> AddZ = {fReader, "AddZ"};
   TTreeReaderArray<ULong64_t> AddTS = {fReader, "AddTS"};
   TTreeReaderValue<ULong64_t> TSHit = {fReader, "TSHit"};


   Selector(TTree * /*tree*/ =0) { }
   virtual ~Selector() { }
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

   TH1D* hspec;
   TH2D* mgamma_gamma;

   ClassDef(Selector,0);

};

#endif

#ifdef Selector_cxx
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


#endif // #ifdef Selector_cxx
