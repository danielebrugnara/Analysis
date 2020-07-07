//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Apr 29 12:46:37 2020 by ROOT version 6.20/04
// from TTree PhysicsTree/Data created / analysed with the NPTool package
// found on file: simu_ana.root
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

#include <TVector3.h>
#include <TH1D.h>
#include <TH2D.h>

// Headers needed by this particular selector
#include "TCATSPhysics.h"

#include "TMust2Physics.h"

#include "TMugastPhysics.h"



class Selector : public TSelector {
public :
   TTreeReader     fReader;  //!the tree reader
   TTree          *fChain = 0;   //!pointer to the analyzed TTree or TChain

   // Readers to access the data (delete the ones you do not need).
   TTreeReaderValue<TCATSPhysics> CATS = {fReader, "CATS"};
   TTreeReaderValue<TMust2Physics> MUST2 = {fReader, "MUST2"};
   TTreeReaderValue<TMugastPhysics> Mugast = {fReader, "Mugast"};
   TTreeReaderValue<Int_t> Nev = {fReader, "Nev"};
   TTreeReaderArray<Double_t> Ex = {fReader, "Ex"};
   TTreeReaderArray<Double_t> ELab = {fReader, "ELab"};
   TTreeReaderArray<Double_t> EDep = {fReader, "EDep"};
   TTreeReaderArray<Double_t> ThetaLab = {fReader, "ThetaLab"};
   TTreeReaderArray<Double_t> ThetaCM = {fReader, "ThetaCM"};
   TTreeReaderValue<Int_t> Run = {fReader, "Run"};
   TTreeReaderArray<Double_t> X = {fReader, "X"};
   TTreeReaderArray<Double_t> Y = {fReader, "Y"};
   TTreeReaderArray<Double_t> Z = {fReader, "Z"};
   TTreeReaderArray<Double_t> dE = {fReader, "dE"};
   TTreeReaderArray<Double_t> dTheta = {fReader, "dTheta"};

   //Needed for this particular selector
   TH1D* angdistr;
   TH1D* ex;
   TH2D* ex_theta;
   TH1D* theta_cm;
   TH2D* kinematicline;
   TVector3 position;
   std::vector<int> MG_nr = {1, 3, 4, 5, 7, 11};
   std::map<int, TH2D *> strip_E;


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
