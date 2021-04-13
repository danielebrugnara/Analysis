//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Apr 13 16:23:12 2021 by ROOT version 6.20/00
// from TTree PhysicsTree/Data created / analysed with the NPTool package
// found on file: Data/anaout/tmp_all.root
//////////////////////////////////////////////////////////

#ifndef Selector_h
#define Selector_h

#include <unordered_map>
#include <map>

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>

// Headers needed by this particular selector
#include "TMust2Physics.h"

#include "TMugastPhysics.h"



class Selector : public TSelector {
public :
   TTreeReader     fReader;  //!the tree reader
   TTree          *fChain = 0;   //!pointer to the analyzed TTree or TChain

   // Readers to access the data (delete the ones you do not need).
   TTreeReaderValue<TMust2Physics> MUST2 = {fReader, "MUST2"};
   TTreeReaderValue<TMugastPhysics> Mugast = {fReader, "Mugast"};
   TTreeReaderValue<Int_t> Nev = {fReader, "Nev"};
   TTreeReaderArray<Double_t> Ex = {fReader, "Ex"};
   TTreeReaderArray<Double_t> ELab = {fReader, "ELab"};
   TTreeReaderArray<Double_t> EDep = {fReader, "EDep"};
   TTreeReaderArray<Double_t> ThetaLab = {fReader, "ThetaLab"};
   TTreeReaderArray<Double_t> PhiLab = {fReader, "PhiLab"};
   TTreeReaderArray<Double_t> ThetaCM = {fReader, "ThetaCM"};
   TTreeReaderValue<Int_t> Run = {fReader, "Run"};
   TTreeReaderArray<Double_t> X = {fReader, "X"};
   TTreeReaderArray<Double_t> Y = {fReader, "Y"};
   TTreeReaderArray<Double_t> Z = {fReader, "Z"};
   TTreeReaderArray<Double_t> dE = {fReader, "dE"};
   TTreeReaderArray<Double_t> dTheta = {fReader, "dTheta"};
   TTreeReaderArray<Double_t> DSSD_X = {fReader, "DSSD_X"};
   TTreeReaderArray<Double_t> DSSD_Y = {fReader, "DSSD_Y"};
   TTreeReaderArray<Double_t> TelescopeNr = {fReader, "TelescopeNr"};

   struct Id{
    int x;
    int y;
    int mg;
   };

   struct PosN{
    struct Pos{
      double x;
      double y;
      double z;
    };
    Pos pos;
    int n;
   };

   typedef std::map<Id, std::map<PosN::Pos, int>> Transform;
 
   Transform transform;

   Transform GetTsf(){return transform; };


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


#endif // #ifdef Selector_cxx
