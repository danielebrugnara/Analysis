//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Jul  6 12:53:07 2020 by ROOT version 6.20/04
// from TTree SimulatedAgata/SimulatedAgata
// found on file: Gamma1.root
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
//#include <TVirtualFFT.h>
#include <TRandom.h>



#include <TVector3.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TF1.h>

#include <RooRealVar.h>
#include <RooPlot.h>
#include <RooProduct.h>
#include <RooExtendPdf.h>
#include <RooGenericPdf.h>
#include <RooDataHist.h>
#include <RooDataSet.h>
#include <RooChebychev.h>
#include <RooAddPdf.h>
#include <RooProdPdf.h>
#include <RooPolynomial.h>

// Headers needed by this particular selector
#include <utility>
#include <vector>
#include <unordered_map>
#include <memory>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <unordered_map>

#include "Particle.h"
#include "Hit.h"
#include <TH2.h>
#include <TStyle.h>

#include <DiaGraph.h>

class Selector : public TSelector {
public :
   TTreeReader     fReader;  //!the tree reader
   TTree          *fChain = 0;   //!pointer to the analyzed TTree or TChain

   // Readers to access the data (delete the ones you do not need).
   TTreeReaderArray<Particle> Part = {fReader, "Particle"};
   TTreeReaderArray<Hit> Hits = {fReader, "Hits"};


   Selector(TTree * /*tree*/ =0):   beta(0.112953),
                                    velocity(33.8625),
                                    t12_1(1.1),
                                    t12_2(6.3+t12_1),
                                    em_position_1(0., 0., velocity*t12_1),
                                    em_position_2(0., 0., velocity*t12_2),
                                    agata_shift(0., 0., 0){ }
                                    //agata_shift(0., 0., 26){ }

    virtual ~Selector() {delete addback_graph; }//Do not delete histograms
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

    double ComputeDoppler(const TVector3 &, const double &) const;
    double ComputeDoppler(const TVector3 &, const TVector3 &, const double &) const;
    double ComputeDoppler(const TVector3 &, const TVector3 &, const double &, const double&) const;
    double ComputeDoppler(const TVector3 &, const double &, const double &) const;

    std::map<int, double> crystalEfficiency;
    double averageCrystalEfficiency{0};
    std::map<int, double> crystalEfficiencyDeclared;
    double averageDeclaredCrystalEfficiency{0};
    void ReadEfficiencyMap();
   //ClassDef(Selector,0);

    TH3D * geom;
    TH3D * geom_abovethr;
    TH1D * cal_spec;
    TH1D * core_spec;
    TH1D * core_spec_DC;
    TH1D * core_spec_DC_pos_1;
    TH1D * core_spec_DC_pos_2;
    TH2D * core_gg;
    TH2D * core_gg_DC;
    TH2D * core_gg_DC_pos_1;
    TH2D * core_gg_DC_pos_2;
    TH1D * addb_spec;
    TH1D * addb_spec_DC;
    TH1D * addb_spec_DC_pos_1;
    TH1D * addb_spec_DC_pos_2;
    std::unordered_map<double, std::unordered_map<double,TH1D *>> addb_spec_DC_scan;
    TH2D * addb_gg;
    TH2D * addb_gg_DC;
    TH2D * addb_gg_DC_pos_1;
    TH2D * addb_gg_DC_pos_2;
    TH1D * target_spec;
    std::vector<TH1D*> crystal_spectra;

    const double beta;
    const double velocity; //mm/ns
    const double t12_1; //ns
    const double t12_2; //ns
    const TVector3 em_position_1; //mm/ns
    const TVector3 em_position_2; //mm/ns
    const TVector3 agata_shift; //mm

    struct Crystal{
        int ID;
        std::string name;
        Crystal(int ID, std::string name): ID(ID), name(std::move(name)){};
        Crystal(): ID(-1){};
        friend std::ostream& operator<<(std::ostream& os, const Crystal& cr)
        {
            os << "{" << cr.ID << ", " << cr.name << "}";
            return os;
        }
    };
    DiaGraph<bool,Crystal>* addback_graph;
    std::unordered_map<int, bool> active_crystals;

    std::unordered_map<int, int> simu_to_data;
    std::unordered_map<int, int> data_to_simu;
    bool use_threasholds{false};
    bool use_intrinsic_efficiencies{true};
    bool compute_threasholds{false};
    std::unordered_map<int,TH1*> threasholds;
    TRandom rand;


    double energy_check_tolerance {0.01};
    double hit_pattern_thr {300};
    //double hit_pattern_thr {1000};
private:
    void ComputeThreasholds();
};


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
#endif


#endif // #ifdef Selector_cxx
