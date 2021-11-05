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

#include <MugastData.h>
#include <ReactionReconstruction.h>
#include <Units.h>
#include <EnergyLoss.h>

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH2D.h>
#include <TFile.h>
#include <TSelector.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>

#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>

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
   TTreeReaderValue<Int_t> ChargeState = {fReader, "ChargeState"};
   TTreeReaderValue<Int_t> VamosAcceptedth = {fReader, "VamosAcceptedth"};
   TTreeReaderValue<Int_t> VamosAcceptedph = {fReader, "VamosAcceptedph"};
   TTreeReaderValue<Int_t> VamosAccepteddelta = {fReader, "VamosAccepteddelta"};
   TTreeReaderArray<Double_t> Ex = {fReader, "Ex"};
   TTreeReaderArray<Double_t> ELab = {fReader, "ELab"};
   TTreeReaderArray<Double_t> EDep = {fReader, "EDep"};
   TTreeReaderArray<Double_t> ThetaLab = {fReader, "ThetaLab"};
   //TTreeReaderArray<Double_t> PhiLab = {fReader, "PhiLab"};
   TTreeReaderArray<Double_t> ThetaCM = {fReader, "ThetaCM"};
   TTreeReaderValue<Int_t> Run = {fReader, "Run"};
   TTreeReaderArray<Double_t> X = {fReader, "X"};
   TTreeReaderArray<Double_t> Y = {fReader, "Y"};
   TTreeReaderArray<Double_t> Z = {fReader, "Z"};
   TTreeReaderArray<Double_t> dE = {fReader, "dE"};
   TTreeReaderArray<Double_t> dTheta = {fReader, "dTheta"};
   TTreeReaderArray<Int_t> DSSD_X = {fReader, "DSSD_X"};
   TTreeReaderArray<Int_t> DSSD_Y = {fReader, "DSSD_Y"};
   TTreeReaderArray<Int_t> TelescopeNr = {fReader, "TelescopeNr"};

   TVector3 targetPos;
   ReactionReconstruction2body<long double>* reaction;

   struct Id{
   public:
    int x;
    int y;
    int mg;
   public:
    Id(int x=0, int y=0, int mg=0): x(x), y(y), mg(mg){};
    bool operator ==  (const Id& rhs)const{
        return x == rhs.x && y == rhs.y && mg == rhs.mg;
    }
    bool operator !=  (const Id& rhs)const{
        return ! (*this == rhs);
    }
    bool operator <  (const Id& rhs)const{
        if(x != rhs.x) return x < rhs.x;
        if(y != rhs.y) return y < rhs.y;
        if(mg != rhs.mg) return mg < rhs.mg;
        return false;
    }
    friend class boost::serialization::access;

    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & x;
        ar & y;
        ar & mg;
    }
   };

    struct Pos{
    public:
      double x;
      double y;
      double z;
    public:
      Pos(double x=0, double y=0, double z=0):x(x), y(y), z(z){};
      bool operator == (const Pos& rhs)const{
          return x == rhs.x && y == rhs.y && z == rhs.z;
      }
        bool operator != (const Pos& rhs)const{
            return ! (*this == rhs);
        }
      bool operator <  (const Pos& rhs)const{
          if(x != rhs.x) return x < rhs.x;
          if(y != rhs.y) return y < rhs.y;
          if(z != rhs.z) return z < rhs.z;
          return false;
      }
      double dist(const Pos& other, const double& scale)const{
          return sqrt(pow(x-other.x*scale, 2)+pow(y-other.y*scale, 2)+pow(z-other.z*scale, 2));
      }
        friend class boost::serialization::access;

        template<class Archive>
        void serialize(Archive & ar, const unsigned int version)
        {
            ar & x;
            ar & y;
            ar & z;
        }
    };

   typedef std::map<Id, std::map<Pos, int>> Transform;

   Transform transform;
   std::map<Id, double> threasholdStripT;
    std::map<Id, double> threasholdStripE;

   Transform GetTsf(){return transform; };
   void SetThreasholds(const std::map<Id, double> thrT, const std::map<Id, double> thrE){threasholdStripE = thrE; threasholdStripT = thrT;};
   void SetOutputName(const std::string& name){outputFileName = name;};
    std::string outputFileName{"out.root"};

   TTree* tree;
   TFile* outfile;
   MugastData mugastData;

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
