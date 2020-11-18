#pragma once

//Cpp headers
#include <map>
#include <string>
#include <vector>
#include <fstream>

//ROOT headers
#include <TFile.h>
#include <TSelector.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TTreeReaderValue.h>
#include <TMath.h>
#include <TH1D.h>
#include <TH2.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TLorentzVector.h>

// External libraries
#include "TCATSPhysics.h"
#include "TMugastPhysics.h"
#include "TMust2Physics.h"

// Internal libraries
#include "Interpolation.h"
#include "VamosIdentification.h"
#include "MugastIdentification.h"
#include "AgataProcessing.h"




class Selector : public TSelector{
public:
    TTreeReader fReader; //!the tree reader
    TTree *fChain = nullptr;   //!pointer to the analyzed TTree or TChain
    long long int total_entries{};

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
    //TTreeReaderArray<UShort_t> T_FPMW__HF =       {fReader, "T_FPMW_HF"};
    //TTreeReaderArray<Float_t> T_FPMW_HF_C =       {fReader, "T_FPMW_HF_C"};
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

    Selector(TTree * /*tree*/ = 0);
    virtual ~Selector();
    virtual Int_t Version() const { return 2; }
    virtual void Begin(TTree *tree);
    virtual void SlaveBegin(TTree *tree);
    virtual void Init(TTree *tree);
    virtual Bool_t Notify();
    virtual Bool_t Process(Long64_t entry);
    virtual Int_t GetEntry(Long64_t entry, Int_t getall = 0) { return fChain ? fChain->GetTree()->GetEntry(entry, getall) : 0; }
    virtual void SetOption(const char *option) { fOption = option; }
    virtual void SetObject(TObject *obj) { fObject = obj; }
    virtual void SetInputList(TList *input) { fInput = input; }
    virtual TList *GetOutputList() const { return fOutput; }
    virtual void SlaveTerminate();
    virtual void Terminate();

    std::vector<std::pair<double, double>> GetTWvsIce();
    std::vector<TObject*> Output;

    //My methods
    inline void LoadVamosData();
    inline void LoadMugastData();
    inline void LoadAgataData();

    inline void PlotVamosGraphs();
    inline void PlotMugastGraphs();
    inline void PlotAgataGraphs();
    inline void PlotCatsGraphs();

    //Constants
    const Long64_t TS_TO_S = 1E8;
    const Long64_t TS_TO_HR = 1E8 * 3600;

    //Jolly pointer to fill histogram
    std::unique_ptr<TH3D> general_histo_ptr{};

    std::string file_name;

    //Analysis classes/////////////////////////////////////////////////////////////////////////////////
    VamosIdentification vamos_fragment;
    MugastIdentification mugast_fragment;
    AgataProcessing agata_gammas;

    struct VamosConf{ //Configuration spectra
        std::unique_ptr<TH2D> mdE_E;
        std::unique_ptr<TH2D> mdE2_E;
        std::unordered_map<int, std::unique_ptr<TH2D>> mQ_MQ; //map index over [Z]
        std::unordered_map<int, std::unique_ptr<TH2D>> Xf_MQ; //map index over [Z]
        //std::unordered_map<int, std::unordered_map<TH1D *>> hAmass;
    };

    struct VamosData{ //Data spectra
        std::unordered_map<int, std::unordered_map<int, std::unique_ptr<TH2D>>> mTW_Brho; //[M][Z]
        std::unordered_map<int, std::unordered_map<int,std::unordered_map<std::string, std::unique_ptr<TH2D>>>> mE_Theta; //[M][Z]
    };

    //AGATA////////////////////////////////////////////////////////////////////////////////////
    struct AgataConf{ //
        std::unique_ptr<TH3D> mmAGATA3D;
        std::unique_ptr<TH1D> hAddTS_LTS;
    };

    std::vector<std::string> AGATAconditions = {"NONE", "MUGASTtap", "MUGASTan"};

    struct AgataData{
        std::unordered_map<int, std::unordered_map<int, std::unordered_map<std::string, std::unique_ptr<TH1D>>>> hDC;            //[mass][nucleus][condition]
        std::unordered_map<int, std::unordered_map<int, std::unique_ptr<TH2D>>> mDC;                                             //[mass][nucleus]
        std::unordered_map<int, std::unordered_map<int, std::unordered_map<std::string, std::unique_ptr<TH2D>>>> mELab_ThetaLab; //[mass][nucleus][particle]
    };

    //Cats//////////////////////////////////////////////////////////////////////////////////////
    struct CatsConf{ //
        std::unique_ptr<TH2D> mCATSpos;
    };

    struct CatsData{};

    //Mugast////////////////////////////////////////////////////////////////////////////////////
    std::vector<std::string> alpha_correction = {"70", "80", "90", "100", "110", "120", "130", "140", "150", "160", "170"};
    std::vector<std::string> particles = {"p", "d", "ANY"};
    std::vector<std::string> gammas = {"360", "1319", "1660", "NOCONDITION"};
    double gamma_gate = 30.;

    struct MGConf{
        std::unique_ptr<TH3D> hit;
        std::unique_ptr<TH2D> hit_XY;
        std::unique_ptr<TH2D> hit_XZ;
        std::unique_ptr<TH2D> hit_YZ;
        std::unique_ptr<TH2D> hit_ThetaPhi;
        std::unordered_map<int, std::unique_ptr<TH2D>> mE_TOF;                                     //[M*#]
        std::unordered_map<int, std::unique_ptr<TH2D>> mE_TOF2;                                    //[M*#]
        std::unordered_map<int, std::unordered_map<std::string, std::unique_ptr<TH2D>>> mStrip_E;  //[M*#][X\Y]
        std::unordered_map<int, std::unordered_map<std::string, std::unique_ptr<TH2D>>> mStrip_T;  //[M*#][X\Y]
        std::unordered_map<int, std::unordered_map<std::string, std::unique_ptr<TH2D>>> mStrip_T2; //[M*#][X\Y]
        std::unordered_map<int, std::unordered_map<int, std::unordered_map<std::string, std::unordered_map<std::string, std::unique_ptr<TH2D>>>>> mELab_ESI;
        std::unordered_map<int, std::unordered_map<int, std::unordered_map<std::string, std::unordered_map<std::string, std::unique_ptr<TH2D>>>>> mThetaLab_ELost;
    };

    struct MGData{
        std::unordered_map<int, std::unordered_map<int, std::unordered_map<int, std::unique_ptr<TH2D>>>> mE_TOF;      //[M][Z][M*#]
        std::unordered_map<int, std::unordered_map<int, std::unordered_map<int, std::unique_ptr<TH2D>>>> mE_TOF2;     //[M][Z][M*#]
        std::unordered_map<int, std::unordered_map<int, std::unordered_map<std::string, std::unique_ptr<TH1D>>>> hEx; //[M][Z][Particle]
        //std::unordered_map<int, std::unordered_map<int, std::unordered_map<std::string, std::unique_ptr<TH2D>>>> mECM_ThetaCM;
        std::unordered_map<int, std::unordered_map<int, std::unordered_map<std::string, std::unordered_map<std::string, std::unique_ptr<TH2D>>>>> mELab_ThetaLab;
        std::unordered_map<int, std::unordered_map<int, std::unordered_map<std::string, std::unordered_map<std::string, std::unique_ptr<TH2D>>>>> mEx_ThetaLab;
        std::unordered_map<int, std::unordered_map<int, std::unordered_map<std::string, std::unordered_map<std::string, std::unique_ptr<TH2D>>>>> mEx_Phi;
        std::unordered_map<int, std::unordered_map<int, std::unordered_map<std::string, std::unordered_map<std::string, std::unique_ptr<TH1D>>>>> hThetaCM;
        std::unordered_map<int, std::unordered_map<int, std::unordered_map<std::string, std::unique_ptr<TH2D>>>> mEx_EDC;
    };

    struct MMConf{
        std::unordered_map<int, std::unique_ptr<TH2D>> mE_TOF;                                    //[M*#]
        std::unordered_map<int, std::unique_ptr<TH2D>> mdE_E_Si;                                  //[MM#]
        std::unordered_map<int, std::unordered_map<std::string, std::unique_ptr<TH2D>>> mStrip_E; //[M*#][X\Y]
        std::unordered_map<int, std::unordered_map<std::string, std::unique_ptr<TH2D>>> mStrip_T; //[M*#][X\Y]
    };

    struct MMData{
        std::unordered_map<int, std::unordered_map<int, std::unordered_map<std::string, std::unique_ptr<TH2D>>>> mE_TOF;   //[Mass][Nucl][M*#]
        std::unordered_map<int, std::unordered_map<int, std::unordered_map<std::string, std::unique_ptr<TH2D>>>> mdE_E_Si; //[Mass][Nucl][M*#]
        std::unordered_map<int, std::unordered_map<int, std::unordered_map<std::string, std::unique_ptr<TH1D>>>> hEx;      //[Mass][Nucl][Parcle]
        std::unordered_map<int, std::unordered_map<int, std::unordered_map<std::string, std::unique_ptr<TH2D>>>> mEx_TW;   //[Mass][Nucl][Parcle]
        std::unordered_map<int, std::unordered_map<int, std::unordered_map<std::string, std::unique_ptr<TH2D>>>> mECM_ThetaCM;
        std::unordered_map<int, std::unordered_map<int, std::unordered_map<std::string, std::unordered_map<std::string, std::unique_ptr<TH2D>>>>> mELab_ThetaLab;
    };

    //Conf
    struct Conf{ //All the conf plots
        VamosConf VAMOS;
        CatsConf CATS;
        AgataConf AGATA;
        MGConf MG;
        MMConf MM;
    };

    struct Data{ //All the data plots
        VamosData VAMOS;
        CatsData CATS;
        AgataData AGATA;
        MGData MG;
        MMData MM;
    };

    Conf pConf;
    Data pData;
    std::unique_ptr<TTree> tree;

private:
    //struct HistogramSettings{
    //    std::vector<std::string> lims;
    //    bool enabled;
    //    HistogramSettings(bool enabled): enabled(enabled){};
    //};
    //std::unordered_map<std::string, HistogramSettings> enabled_histograms;
    std::unordered_map<std::string, bool> enabled_histograms;
    std::ofstream *new_graph_file = nullptr;
    bool GetSettings();

    template <class T>
    bool Istantiate(std::unique_ptr<T>& ptr, T *ptr_value){
        if (enabled_histograms.find(ptr_value->GetName()) != enabled_histograms.end()){
            ptr.reset(ptr_value);
            Output.push_back(ptr.get());
            //auto const *lims = &enabled_histograms[ptr_value->GetName()].lims;
            //if (lims->size()%3){
            //    for (int ii=0; ii<lims->size(); ii+=3){
            //        //ptr->GetXaxis();
            //    }
            //}
        }else if (new_graph_file != nullptr){
            *new_graph_file << ptr_value->GetName()
                            << "\t\t\t\t\t"
                            << "false"
                            << std::endl;
        }else{
            delete ptr_value;
            ptr.reset();
        }
        return true;
    };

    static inline bool Fill(TH1D* , const Double_t &);
    static inline bool Fill(TH2D* , const Double_t &, const Double_t &);
    static inline bool Fill(TH3D* , const Double_t &, const Double_t &, const Double32_t &);
    static inline bool Fill(TTree*);

    //ClassDef(Selector,0);
};