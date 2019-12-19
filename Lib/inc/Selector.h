//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Jul 23 16:53:27 2019 by ROOT version 6.16/00
// from TChain PhysicsTree/
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

#include <TMath.h>
// Headers needed by this particular selector
#include "TCATSPhysics.h"
#include "TMust2Physics.h"
#include "TMugastPhysics.h"

#include "Interpolation.h"

#include <vector>
#include <string>
#include <map>
#include <TMath.h>

#include <fstream>

#include <TH2.h>
#include <TH3D.h>
#include <TH2D.h>
#include <TH1D.h>
#include <TCutG.h>
#include <TStyle.h>

#include <TLorentzVector.h>

#include <TSystem.h>
#include <TGClient.h>

class Selector : public TSelector {
public :
   TTreeReader     fReader;  //!the tree reader
   TTree          *fChain = 0;   //!pointer to the analyzed TTree or TChain

   // Readers to access the data (delete the ones you do not need).
   TTreeReaderValue<TCATSPhysics> CATS = {fReader, "CATS"};
   TTreeReaderValue<TMust2Physics> MUST2 = {fReader, "MUST2"};
   TTreeReaderValue<Double_t> CONFDEC_AGATA = {fReader, "CONFDEC_AGATA"};
   TTreeReaderValue<Double_t> CONFDEC_VAMOS = {fReader, "CONFDEC_VAMOS"};
   TTreeReaderValue<Double_t> DATATRIG_CATS1 = {fReader, "DATATRIG_CATS1"};
   TTreeReaderValue<Double_t> DATATRIG_CATS2 = {fReader, "DATATRIG_CATS2"};
   TTreeReaderValue<Double_t> ECATS_AGATA = {fReader, "ECATS_AGATA"};
   TTreeReaderValue<Double_t> EVAMOS_AGATA = {fReader, "EVAMOS_AGATA"};
   TTreeReaderValue<Double_t> E_AGATA = {fReader, "E_AGATA"};
   TTreeReaderValue<Double_t> GATCONF_MASTER = {fReader, "GATCONF_MASTER"};
   TTreeReaderValue<Double_t> GATCONF_SLAVE1 = {fReader, "GATCONF_SLAVE1"};
   TTreeReaderValue<Double_t> GATCONF_SLAVE2 = {fReader, "GATCONF_SLAVE2"};
   TTreeReaderValue<Double_t> PLQ_L = {fReader, "PLQ_L"};
   TTreeReaderValue<Double_t> PLQ_R = {fReader, "PLQ_R"};
   TTreeReaderValue<Double_t> TCATS_MUGAST_HF = {fReader, "TCATS_MUGAST_HF"};
   TTreeReaderValue<Double_t> TMG_MUGAST_HF = {fReader, "TMG_MUGAST_HF"};
   TTreeReaderValue<Double_t> TVAMOS_MUGAST_HF = {fReader, "TVAMOS_MUGAST_HF"};
   TTreeReaderValue<Double_t> T_Agata_HF = {fReader, "T_Agata_HF"};
   TTreeReaderValue<Double_t> T_Cats1_2 = {fReader, "T_Cats1_2"};
   TTreeReaderValue<Double_t> T_Cats1_HF = {fReader, "T_Cats1_HF"};
   TTreeReaderValue<Double_t> T_Fag_HF = {fReader, "T_Fag_HF"};
   TTreeReaderValue<Double_t> T_Si_Agata = {fReader, "T_Si_Agata"};
   TTreeReaderValue<Double_t> T_Si_Cats1 = {fReader, "T_Si_Cats1"};
   TTreeReaderValue<Double_t> T_Si_Vamos = {fReader, "T_Si_Vamos"};
   TTreeReaderValue<Double_t> T_Vamos_Agata = {fReader, "T_Vamos_Agata"};
   TTreeReaderValue<Double_t> T_Vamos_Cats1 = {fReader, "T_Vamos_Cats1"};
   TTreeReaderValue<TMugastPhysics> Mugast = {fReader, "Mugast"};
   TTreeReaderArray<double> Ex = {fReader, "Ex"};
   TTreeReaderArray<int> DetID = {fReader, "DetID"};
   TTreeReaderArray<double> EDC = {fReader, "EDC"};
   TTreeReaderValue<Double_t> AddBack_EDC = {fReader, "AddBack_EDC"};
   TTreeReaderValue<Double_t> EAgata = {fReader, "EAgata"};
   TTreeReaderArray<double> ELab = {fReader, "ELab"};
   TTreeReaderArray<double> Ecm = {fReader, "Ecm"};
   TTreeReaderArray<double> ThetaLab = {fReader, "ThetaLab"};
   TTreeReaderArray<double> PhiLab = {fReader, "PhiLab"};
   TTreeReaderArray<double> ThetaCM = {fReader, "ThetaCM"};
   TTreeReaderValue<Int_t> Run = {fReader, "Run"};
   TTreeReaderArray<double> X = {fReader, "X"};
   TTreeReaderArray<double> Y = {fReader, "Y"};
   TTreeReaderArray<double> Z = {fReader, "Z"};
   TTreeReaderValue<Double_t> dE = {fReader, "dE"};
   TTreeReaderValue<ULong64_t> LTS = {fReader, "LTS"};
   TTreeReaderValue<UShort_t> T_FPMW_CATS1 = {fReader, "T_FPMW_CATS1"};
   TTreeReaderValue<Float_t> T_FPMW_CATS2_C = {fReader, "T_FPMW_CATS2_C"};
   //TTreeReaderArray<UShort_t> T_FPMW__HF = {fReader, "T_FPMW_HF"};
   //TTreeReaderArray<Float_t> T_FPMW_HF_C = {fReader, "T_FPMW_HF_C"};
   TTreeReaderArray<Float_t> IC = {fReader, "IC"};
   TTreeReaderArray<UShort_t> ICRawPUMult01 = {fReader, "ICRawPUMult01"};
   TTreeReaderValue<UShort_t> T_CATS2_HF = {fReader, "T_CATS2_HF"};
   TTreeReaderValue<UShort_t> T_MUGAST_FPMW = {fReader, "T_MUGAST_FPMW"};
   TTreeReaderValue<UShort_t> T_FPMW_CATS2 = {fReader, "T_FPMW_CATS2"};
   TTreeReaderValue<Float_t> DC0_X = {fReader, "DC0_X"};
   TTreeReaderValue<Float_t> DC0_Y = {fReader, "DC0_Y"};
   TTreeReaderValue<Float_t> DC1_X = {fReader, "DC1_X"};
   TTreeReaderValue<Float_t> DC1_Y = {fReader, "DC1_Y"};
   TTreeReaderValue<Float_t> DC2_X = {fReader, "DC2_X"};
   TTreeReaderValue<Float_t> DC2_Y = {fReader, "DC2_Y"};
   TTreeReaderValue<Float_t> DC3_X = {fReader, "DC3_X"};
   TTreeReaderValue<Float_t> DC3_Y = {fReader, "DC3_Y"};
   TTreeReaderValue<Float_t> Xf = {fReader, "Xf"};
   TTreeReaderValue<Float_t> Tf = {fReader, "Tf"};
   TTreeReaderValue<Float_t> PhiL = {fReader, "PhiL"};
   TTreeReaderValue<Float_t> ThetaL = {fReader, "ThetaL"};
   TTreeReaderValue<Float_t> Brho = {fReader, "Brho"};
   TTreeReaderValue<Float_t> TW = {fReader, "TW"};
   TTreeReaderValue<Float_t> Theta = {fReader, "Theta"};
   TTreeReaderValue<Float_t> Phi = {fReader, "Phi"};
   TTreeReaderValue<Float_t> Path = {fReader, "Path"};
   TTreeReaderValue<UShort_t> EWIRE_1_1 = {fReader, "EWIRE_1_1"};
   TTreeReaderValue<UShort_t> EWIRE_1_2 = {fReader, "EWIRE_1_2"};
   TTreeReaderValue<UShort_t> EWIRE_2_1 = {fReader, "EWIRE_2_1"};
   TTreeReaderValue<UShort_t> EWIRE_2_2 = {fReader, "EWIRE_2_2"};
   TTreeReaderValue<Int_t> MW_Nr = {fReader, "MW_Nr"};
   TTreeReaderArray<Float_t> MW_T = {fReader, "MW_T"};
   TTreeReaderArray<UShort_t> MW_N = {fReader, "MW_N"};
   TTreeReaderValue<ULong64_t> AGAVA_VAMOSTS = {fReader, "AGAVA_VAMOSTS"};
   TTreeReaderValue<Double_t> mT = {fReader, "mT"};
   TTreeReaderValue<Double_t> mV = {fReader, "mV"};
   TTreeReaderValue<Double_t> mD = {fReader, "mD"};
   TTreeReaderValue<Double_t> mBeta = {fReader, "mBeta"};
   TTreeReaderValue<Double_t> mGamma = {fReader, "mGamma"};
   TTreeReaderValue<Double_t> mM_Q = {fReader, "mM_Q"};
   TTreeReaderValue<Double_t> mM = {fReader, "mM"};
   TTreeReaderValue<Double_t> mE = {fReader, "mE"};
   TTreeReaderValue<Double_t> mdE = {fReader, "mdE"};
   TTreeReaderValue<Int_t> nbAdd = {fReader, "nbAdd"};
   TTreeReaderValue<ULong64_t> TSHit = {fReader, "TSHit"};
   TTreeReaderValue<ULong64_t> AddTS = {fReader, "AddTS"};
   TTreeReaderArray<Float_t> AddE = {fReader, "AddE"};
   TTreeReaderArray<Float_t> AddX = {fReader, "AddX"};
   TTreeReaderArray<Float_t> AddY = {fReader, "AddY"};
   TTreeReaderArray<Float_t> AddZ = {fReader, "AddZ"};


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

   //My methods
   inline void IdentifyFragment();
   inline void FillMugastConfHistograms();
   void CheckCutsPresence();

   //Constants
   const Long64_t TS_TO_S = 1E8;
   const Long64_t TS_TO_HR = 1E8*3600;
   const Double_t AMU_TO_MEV = 931.4936148;
   const Double_t c_speed = 299792458; //m/s
   const Double_t MEV_TO_J = 1.60217733e-13;

   //Graphs
   //General
   //bool debug = false;
   //int count46Ar;
   struct Counter{
      long long General;
      long long Ar46;
      long long K47;
      Counter(): General(0), Ar46(0), K47(0){};
   } counter;


   std::string file_name;


   //Reactions


   //Target
   Interpolation *thickness_angle;
   Interpolation *angle_angle;

   //VAMOS/////////////////////////////////////////////////////////////////////////////////
   struct VamosId
   { //Identification struct
      double En;
      double D_En;
      double D_En2;
      double Path;
      double T;
      double V;
      double Beta;
      double Gamma;
      double M_Q;
      double M;
      double Charge;
      bool Identified;
      std::string Nucleus;
      TLorentzVector p4;
      std::string GetNucleus()
      {
         if (Nucleus.compare(""))
            return Nucleus.substr(Nucleus.find_last_of("_"));
         else
            return std::string("ANY");
      }
      std::string GetMass()
      {
         if (Nucleus.compare(""))
            return Nucleus.substr(0, Nucleus.find_first_of("_"));
         else
            return std::string("ANY");
      }
      VamosId() : En(0),
                  D_En(0),
                  D_En2(0),
                  Path(0),
                  T(0),
                  V(0),
                  Beta(0),
                  Gamma(0),
                  M_Q(0),
                  M(0),
                  Charge(0),
                  Identified(false),
                  Nucleus(""),
                  p4(){};
   };

   VamosId *IdentifiedNucleus;

   struct VamosConf{//Configuration spectra
      TH2D* mdE_E;
      std::map<std::string,TH2D*> mQ_MQ;
      std::map<std::string, TH1D*> hAmass;
   };

   struct VamosData{//Data spectra
      std::map<std::string,std::map<std::string, int>> evNr;
      std::map<std::string,std::map<std::string, TH2D*>> mTW_Brho;
   };


   //AGATA////////////////////////////////////////////////////////////////////////////////////
   struct AgataConf{
      TH3D *mmAGATA3D;
      TH1D *hAddTS_LTS;
   };

   std::vector<std::string> masses = {"46", "47"};
   std::vector<std::string> nuclei = {"K", "Ar"};

   std::vector<std::string> AGATAconditions = {"NONE", "MUGASTtap", "MUGASTan"};

   struct AgataData{
      std::map<std::string, std::map<std::string,std::map<std::string,TH1D*>>> hDC; //[mass][nucleus][condition]
      std::map<std::string, std::map<std::string,TH2D*>> mDC;//[mass][nucleus]
      std::map<std::string, std::map<std::string,TH2D*>> mDC_ThetaMUGAST;//[mass][nucleus]
      std::map<std::string, std::map<std::string,std::map<std::string, TH2D*>>> mEx_DC;//[mass][nucleus][particle]
      std::map<std::string, std::map<std::string,std::map<std::string, TH2D*>>> mELab_ThetaLab;//[mass][nucleus][particle]
   };


   //Cats//////////////////////////////////////////////////////////////////////////////////////
   struct CatsConf{//
      TH2D* mCATSpos;
   };

   struct CatsData{//
   };

   //Mugast////////////////////////////////////////////////////////////////////////////////////
   std::vector<std::string> silicons = {"MM1", "MM2", "MM3", "MM4", "MG1", "MG2", "MG3", "MG4", "MG5", "MG6", "MG7", "MG8", "MG9", "MG10", "MG11"};
   std::vector<std::string> siliconsMM = {"MM1", "MM2", "MM3", "MM4"};
   std::vector<std::string> siliconsMG = {"MG1", "MG2", "MG3", "MG4", "MG5", "MG6", "MG7", "MG8", "MG9", "MG10", "MG11"};
   std::vector<std::string> strips;
   std::vector<std::string> alpha_correction = {"70","80", "90", "100", "110","120","130","140","150","160","170"};
   std::vector<std::string> particles={"p", "d", "ANY"};
   std::vector<std::string> gammas={"360 keV", "ANY"};

   struct SiConf{//
      std::map<std::string, TH2D*> mE_TOF;   //[M*#] 
     // std::map<std::string, std::map<std::string,TH2D*>> mE_TOF_Corrected;   //[M*#][alpha] 
      std::map<std::string, TH2D*> mdE_E_Si; //[MM#] 
      std::map<std::string, std::map<std::string,TH2D*>>  mStrip_E; //[M*#][X\Y]
      std::map<std::string, std::map<std::string,TH2D*>>  mStrip_T; //[M*#][X\Y]
   };

   struct SiData{//
      std::map<std::string, std::map<std::string, std::map<std::string,TH2D*>>> mE_TOF;//[M*#][Mass][Nucl]
      std::map<std::string, std::map<std::string, std::map<std::string,TH2D*>>> mdE_E_Si;//[M*#][Mass][Nucl]
      std::map<std::string, std::map<std::string, std::map<std::string,TH1D*>>> hEx;//[Mass][Nucl][Parcle]
      std::map<std::string, std::map<std::string, std::map<std::string,TH2D*>>> mEx_TW;//[Mass][Nucl][Parcle]
      std::map<std::string, std::map<std::string, std::map<std::string,TH2D*>>> mECM_ThetaCM;
      std::map<std::string, std::map<std::string, std::map<std::string,std::map<std::string,TH2D*>>>> mELab_ThetaLab;
   };

   inline double AlignT(int detector, int strip, double time){
      if (detector!=11) return time;
      else{
         double reference=356.391;
         switch (strip)
         {
         case 0: return time;
         case 1: return  time;
         case 2: return  time;
         case 3: return  time;
         case 4: return  time;
         case 5: return  time;
         case 6: return  time;
         case 7: return  time;
         case 8: return  time;
         case 9: return  time;
         case 10: return time;
         case 11: return time;
         case 12: return time;
         case 13: return time;
         case 14: return time;
         case 15: return time;
         case 16: return time;
         default: return 0;
         }
      }
   }

   inline double AlignPunch(int detector, double time){
      double reference = 358.9;
      switch(detector){
         case 1: return time; break;
         case 2: return time;
         case 3: return time-361.1+reference;
         case 4: return time-375.3+reference;
         case 5: return time-354.7+reference;
         case 6: return time;
         case 7: return time-350.3+reference;
         case 8: return time;
         case 9: return time;
         case 10: return time;
         case 11: return time-357.1+reference;
         default: return 0;
      }
   }

   //Conf
   struct Conf{
      VamosConf VAMOS;
      CatsConf CATS;
      AgataConf AGATA;
      SiConf SI; 
   };

   struct Data{
      VamosData VAMOS;
      CatsData CATS;
      AgataData AGATA;
      SiData SI; 
   };

   Conf pConf;
   Data pData;



   //Graphical cuts////////////////////////////////////////////////////////////////////////////////////
   std::map<std::string, std::map< std::string, TCutG*>> cut;
   std::vector<std::string> Qcuts = {"Q14", "Q15", "Q16", "Q17", "Q18"};
   std::vector<std::string> Zcuts = {"dE_E_K", "dE_E_Ar"};
   std::vector<std::string> Mgates= {"45", "46", "47"};

   //std::vector<std::string> QcutsAr = {"Q12Ar","Q13Ar","Q14Ar", "Q15Ar", "Q16Ar", "Q17Ar", "Q18Ar"};
   //std::vector<std::string> QcutsK = {"Q13K","Q14K", "Q15K", "Q16K", "Q17K", "Q18K", "Q19K"};
   TFile* VAMOScuts;
   TFile* MUGASTcuts;

   //long long counter;

   std::map<std::string, double> mass; 

   private:
   inline Double_t CorrectDoppler(const TLorentzVector &, const Double_t &, const Double_t &, const Double_t &, const Double_t &);
   inline Double_t CorrectTOF(const Double_t &, const TVector3 &, const Double_t &, const std::string &);
   const Double_t agata_Zshift=-4;
   bool GetSettings();

   inline bool Fill(TH1* , const Double_t &);
   inline bool Fill(TH2* , const Double_t &, const Double_t &);
   inline bool Fill(TH3* , const Double_t &, const Double_t &, const Double32_t &);

   //ClassDef(Selector,0);

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
