#define Selector_cxx

#include "Selector.h"

Selector::Selector(TTree *) {}

Selector::~Selector() {
    TIter iter(fOutput);
    TObject *obj;
    while ((obj = iter())) {
        delete obj;
    }
}

void Selector::Begin(TTree * /*tree*/) {
    // The Begin() function is called at the start of the query.
    // When running with PROOF Begin() is only called on the client.
    // The tree argument is deprecated (on PROOF 0 is passed).
#ifdef VERBOSE_DEBUG
    std::cout << "------------>Selector::Begin()\n";
#endif

    TString option = GetOption();
}

void Selector::SlaveBegin(TTree * /*tree*/) {
    // The SlaveBegin() function is called after the Begin() function.
    // When running with PROOF SlaveBegin() is called on each slave server.
    // The tree argument is deprecated (on PROOF 0 is passed).
#ifdef VERBOSE_DEBUG
    std::cout << "------------>Selector::SlaveBegin()\n";
#endif

    TString option = GetOption();
    file_name = "./Out/" + option;

    //Initializing Histograms

    //Configurations
    //Target
    thickness_angle = new Interpolation("./Configs/Interpolations/GasThickness.txt");
    angle_angle = new Interpolation("./Configs/Interpolations/EntranceAngleHavar.txt");

    //////////////////////////////////////////////////////////////////////////////////////////////////
    ///Histogram pointer initialization//////////////////////////////////////////////////////////////

    //Config histograms//////////////////////////////////////////////////////////////////////////////
    //VAMOS
    pConf.VAMOS.mdE_E = new TH2D("pConf_VAMOS_mdE_E", "dE E in VAMOS", 4000, 10, 350, 4000, 10, 140);
    fOutput->Add(pConf.VAMOS.mdE_E);

    pConf.VAMOS.mdE2_E = new TH2D("pConf_VAMOS_mdE2_E", "dE2 E in VAMOS", 2000, 10, 350, 2000, 10, 140);
    fOutput->Add(pConf.VAMOS.mdE2_E);

    for (const auto &Z_it : vamos_fragment.cuts_Z) {
        pConf.VAMOS.mQ_MQ[Z_it] = new TH2D(Form("pConf_VAMOS_mQ_MQ_M%i", Z_it), Form("M/Q vs Q with Z%i selection", Z_it), 1000, 2, 4, 1000, 3, 24);
        fOutput->Add(pConf.VAMOS.mQ_MQ[Z_it]);
    }

    //    pConf.VAMOS.hAmass["Ar"] = new TH1D("pConf-VAMOS-hAmass-Ar", "Mass histogram of Ar", 1000, 20, 50);
    //    fOutput->Add(pConf.VAMOS.hAmass["Ar"]);
    //    pConf.VAMOS.hAmass["K"] = new TH1D("pConf-VAMOS-hAmass-K", "Mass histogram of K", 1000, 20, 50);
    //    fOutput->Add(pConf.VAMOS.hAmass["K"]);

    //AGATA
    pConf.AGATA.mmAGATA3D = new TH3D("pConf-AGATA-mmAGATA3D", "Hit patter on AGATA", 50, -300, 300, 50, -300, 300, 50, -300, 300);
    fOutput->Add(pConf.AGATA.mmAGATA3D);
    pConf.AGATA.hAddTS_LTS = new TH1D("pConf-AGATA-hAddTS_LTS", "Difference between AddTS and LTS", 1000, 0, 300);
    fOutput->Add(pConf.AGATA.hAddTS_LTS);

    //CATS
    pConf.CATS.mCATSpos = new TH2D("pConf-CATS-mCATSpos", "Beam position on CATS", 1000, -50, 50, 1000, -50, 50);
    fOutput->Add(pConf.CATS.mCATSpos);

    //MUGAST
    for (int ii = 0; ii < 128; ii++) {
        strips.push_back(to_string(ii));
    }

    for (const auto &MM : siliconsMM) {
        pConf.SI.mdE_E_Si[MM] = new TH2D(Form("pConf-SI-mdE_E_Si-%s", MM.c_str()), Form("dE E of %s", MM.c_str()), 1000, 0, 28, 1000, 0, 28);
        fOutput->Add(pConf.SI.mdE_E_Si[MM]);
    }
    for (const auto &SI : silicons) {
        pConf.SI.mE_TOF[SI] = new TH2D(Form("pConf-SI-mE_TOF-%s", SI.c_str()), Form("E vs TOF of %s", SI.c_str()), 1000, 0, 28, 1000, 260, 380);
        fOutput->Add(pConf.SI.mE_TOF[SI]);
        pConf.SI.mStrip_E[SI]["X"] = new TH2D(Form("pConf-SI-mStrip_E-%s-%s", SI.c_str(), "X"), Form("Strip vs E of %s %s", SI.c_str(), "X"), 128, 0, 128, 1000, 0, 30);
        fOutput->Add(pConf.SI.mStrip_E[SI]["X"]);
        pConf.SI.mStrip_E[SI]["Y"] = new TH2D(Form("pConf-SI-mStrip_E-%s-%s", SI.c_str(), "Y"), Form("Strip vs E of %s %s", SI.c_str(), "Y"), 128, 0, 128, 1000, 0, 30);
        fOutput->Add(pConf.SI.mStrip_E[SI]["Y"]);
        pConf.SI.mStrip_T[SI]["X"] = new TH2D(Form("pConf-SI-mStrip_T-%s-%s", SI.c_str(), "X"), Form("Strip vs T of %s %s", SI.c_str(), "X"), 128, 0, 128, 1000, 0, 1500);
        fOutput->Add(pConf.SI.mStrip_T[SI]["X"]);
        pConf.SI.mStrip_T[SI]["Y"] = new TH2D(Form("pConf-SI-mStrip_T-%s-%s", SI.c_str(), "Y"), Form("Strip vs T of %s %s", SI.c_str(), "Y"), 128, 0, 128, 1000, 0, 1500);
        fOutput->Add(pConf.SI.mStrip_T[SI]["Y"]);
    }

    //Data histograms//////////////////////////////////////////////////////////////////////////////
    for (const auto &it_M : vamos_fragment.cuts_M) {
        for (const auto &it_Z : vamos_fragment.cuts_Z) {
            pData.VAMOS.mTW_Brho[it_M][it_Z] = new TH2D(Form("pData-VAMOS-mTW_Brho-%i-%i", it_M, it_Z), Form("Time vs Brho with %i %i in VAMOS", it_M, it_Z), 5000, 242, 328, 1000, 0.5, 1.5);
            fOutput->Add(pData.VAMOS.mTW_Brho[it_M][it_Z]);
        }
    }
    //VAMOS

    //AGATA
    for (const auto &it_M : vamos_fragment.cuts_M) {
        for (const auto &it_Z : vamos_fragment.cuts_Z) {
            for (const auto &condition : AGATAconditions) {
                pData.AGATA.hDC[it_M][it_Z][condition] = new TH1D(Form("pData_AGATA_hDC_M%i_Z%i_cond%s", it_M, it_Z, condition.c_str()), Form("DC spectrum of M%i Z%i with condition %s", it_M, it_Z, condition.c_str()), 4000, 0, 4000);
                fOutput->Add(pData.AGATA.hDC[it_M][it_Z][condition]);
            }
            pData.AGATA.mDC[it_M][it_Z] = new TH2D(Form("pData_AGATA_mDC_M%i_Z%i", it_M, it_Z), Form("DC gamma gamma of M%i Z%i", it_M, it_Z), 4000, 0, 4000, 4000, 0, 4000);
            fOutput->Add(pData.AGATA.mDC[it_M][it_Z]);
            pData.AGATA.mDC_ThetaMUGAST[it_M][it_Z] = new TH2D(Form("pData_AGATA_mDC_ThetaMUGAST_M%i_Z%i", it_M, it_Z), Form("DC gamma vs Theta on MUGAST of M%i Z%i", it_M, it_Z), 4000, 0, 4000, 180, 0, 180);
            fOutput->Add(pData.AGATA.mDC_ThetaMUGAST[it_M][it_Z]);
            for (const auto &particle : particles) {
                pData.AGATA.mEx_DC[it_M][it_Z][particle] = new TH2D(Form("pData_AGATA_mEx_DC_M%i_Z%i_%s", it_M, it_Z, particle.c_str()), Form("Excitation energy AGATA vs MUGAST M%i Z%i and %s", it_M, it_Z, particle.c_str()), 1000, 0, 10, 1000, 0, 10);
                fOutput->Add(pData.AGATA.mEx_DC[it_M][it_Z][particle]);
                pData.AGATA.mELab_ThetaLab[it_M][it_Z][particle] = new TH2D(Form("pData_AGATA_ELab_ThetaLab_M%i_Z%i_%s", it_M, it_Z, particle.c_str()), Form("Excitation energy AGATA vs MUGAST M%i Z%i and %s", it_M, it_Z, particle.c_str()), 1000, 0, 10, 1000, 0, 10);
                fOutput->Add(pData.AGATA.mELab_ThetaLab[it_M][it_Z][particle]);
            }
        }
    }


#ifdef VERBOSE_DEBUG
    std::cout << "------------>finished: agata initialization\n";
#endif

    //CATS

    //MUGAST
    for (const auto &it_M : vamos_fragment.cuts_M) {
        for (const auto &it_Z : vamos_fragment.cuts_Z) {
            for (const auto &MM : siliconsMM) {
                pData.SI.mdE_E_Si[it_M][it_Z][MM] = new TH2D(Form("pData_SI_mdE_E_Si_M%i_Z%i_%s", it_M, it_Z, MM.c_str()), Form("dE E of %s with M%i Z%i", MM.c_str(), it_M, it_Z), 1000, 0, 28, 1000, 0, 28);
                fOutput->Add(pData.SI.mdE_E_Si[it_M][it_Z][MM]);
            }
            for (const auto &SI : silicons) {
                pData.SI.mE_TOF[it_M][it_Z][SI] = new TH2D(Form("pData_SI_mE_TOF_M%i_Z%i_%s", it_M, it_Z, SI.c_str()), Form("E vs TOF of %s with M%i Z%i", SI.c_str(), it_M, it_Z), 1000, 0, 28, 1000, 260, 380);
                fOutput->Add(pData.SI.mE_TOF[it_M][it_Z][SI]);
            }
            pData.SI.mE_TOF[it_M][it_Z]["MG"] = new TH2D(Form("pData_SI_mE_TOF_MG_M%i_Z%i", it_M, it_Z), Form("E vs TOF of all MG with M%i Z%i", it_M, it_Z), 1000, 0, 28, 1000, 260, 380);
            fOutput->Add(pData.SI.mE_TOF[it_M][it_Z]["MG"]);
            for (const auto &particle : particles) {
                pData.SI.hEx[it_M][it_Z][particle] = new TH1D(Form("pData_SI-hEx_M%i_Z%i_%s", it_M, it_Z, particle.c_str()), Form("Excitation energy with M%i Z%i in VAMOS and %s in MUGAST", it_M, it_Z, particle.c_str()), 1000, -60, 60);
                fOutput->Add(pData.SI.hEx[it_M][it_Z][particle]);
                pData.SI.mEx_TW[it_M][it_Z][particle] = new TH2D(Form("pData_SI_mEx_TW_M%i_Z%i_%s", it_M, it_Z, particle.c_str()), Form("Excitation energy vs Time with M%i Z%i in VAMOS and %s in MUGAST", it_M, it_Z, particle.c_str()), 5000, 242, 328, 1000, -60, 60);
                fOutput->Add(pData.SI.mEx_TW[it_M][it_Z][particle]);
                pData.SI.mECM_ThetaCM[it_M][it_Z][particle] = new TH2D(Form("pData_SI_mECM_ThetaCM_M%i_Z%i_%s", it_M, it_Z, particle.c_str()), Form("E CM vs Theta CM with M%i Z%i in VAMOS and %s in MUGAST", it_M, it_Z, particle.c_str()), 1000, 0, 180, 1000, 0, 60);
                fOutput->Add(pData.SI.mECM_ThetaCM[it_M][it_Z][particle]);
                for (const auto &gamma : gammas) {
                    pData.SI.mELab_ThetaLab[it_M][it_Z][particle][gamma] = new TH2D(Form("pData_SI_mELab_ThetaLab_M%i_Z%i_%s_%s", it_M, it_Z, particle.c_str(), gamma.c_str()), Form("ELab vs Theta Lab with M%i Z%i in VAMOS and %s in MUGAST and %s in AGATA", it_M, it_Z, particle.c_str(), gamma.c_str()), 1000, 0, 180, 1000, 0, 60);
                    fOutput->Add(pData.SI.mELab_ThetaLab[it_M][it_Z][particle][gamma]);
                }
            }
        }
    }
    for (const auto &particle : particles) {
        pData.SI.mELab_ThetaLab[0][0][particle]["ANY"] = new TH2D(Form("pData-SI-mELab_ThetaLab-%s-%s-%s", "ANY", "ANY", particle.c_str()), Form("E Lab vs Theta Lab with %s %s in VAMOS and %s in MUGAST", "ANY", "ANY", particle.c_str()), 1000, 0, 180, 1000, 0, 60);
        fOutput->Add(pData.SI.mELab_ThetaLab[0][0][particle]["ANY"]);
    }

#ifdef VERBOSE_DEBUG
    std::cout << "------------>finished: SI initialization\n";
#endif

    ///////////////////////////////////////////////////////////////////////////////////////////
    ///Loading things from TFiles/////////////////////////////////////////////////////////////
    std::ifstream *input_file = nullptr;

    //GCuts///////////////////////////////////////////////////////////////////////////////////
    std::string VAMOS_cuts_file = "./Configs/Cuts/VAMOS.root";
    std::string MUGAST_cuts_file = "./Configs/Cuts/MUGAST.root";

    //TODO: use constructor to initialize this!!!
    vamos_fragment.LoadCuts(VAMOS_cuts_file);
    vamos_fragment.Initialize();

#ifdef VERBOSE_DEBUG
    std::cout << "------------>finished: vamos_fragment initialization\n";
#endif

    //TODO: update with new identification class
    input_file = new std::ifstream(MUGAST_cuts_file);
    if (input_file) {
        input_file->close();
        TFile *MUGASTcuts = new TFile(MUGAST_cuts_file.c_str(), "READ");
        if (!(MUGASTcuts->IsOpen())) {
            std::cout << "MUGAST file not opened\n";
        } else {
            TIter contents(MUGASTcuts->GetListOfKeys());
            TKey *key;
            TObject *obj;
            while ((key = (TKey *)contents())) {
                obj = MUGASTcuts->Get(key->GetName());
                //         TClass *cl = gROOT->GetClass(key->GetClassName());
                //         if (cl->InheritsFrom("TCutG")){
                if (obj->InheritsFrom("TCutG")) {
                    TCutG *tmp = (TCutG *)obj;
                    cut["MUGAST"][tmp->GetName()] = tmp;
                    //std::cout << "Found cut in MUGAST :" << tmp->GetName() << std::endl;
                }
            }
        }
    } else {
        std::cout << "MUGAST cuts not found\n";
    }

    //Interpolations///////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////////
    ///Other configs/////////////////////////////////////////////////////////////////////////////////
    mass["2_H"] = 2.01410177812 * AMU_TO_MEV;  //In MeV
    mass["1_H"] = 1.00782503223 * AMU_TO_MEV;  //In MeV

    //////////////////////////////////////////////////////////////////////////////////////////////////
    ///Temporary histogram///////////////////////////////////////////////////////////////////////////
    general_histo_ptr = new TH2D("test", "test", 1000, -400, 400, 1000, 0, 100);
    fOutput->Add(general_histo_ptr);

    //////////////////////////////////////////////////////////////////////////////////////////////////
    ///Deciding which graphs to plot/////////////////////////////////////////////////////////////////
    GetSettings();  //Decides which histograms to fill
#ifdef VERBOSE_DEBUG
    std::cout << "------------>finished: Selector::SlaveBegin()\n";
#endif
}

Bool_t Selector::Process(Long64_t entry) {
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

#ifdef VERBOSE_DEBUG
    std::cout << "------------>Selector::Process()\n";
#endif

    //Agata vs Ancillary coincidence gate
    Bool_t AGATA_GOOD = *AddTS - *LTS > 175 && *AddTS - *LTS < 184;

    //Vamos identification////////////////////////////////////////////////////////////////////////////////////
    //Configuration spectra are filled in the identification

    LoadVamosData();
#ifdef VERBOSE_DEBUG
    std::cout << "------------>Finished: Loading Vamos data\n";
#endif
    vamos_fragment.Identify();
#ifdef VERBOSE_DEBUG
    std::cout << "------------>Finished: Identification\n";
#endif
    PlotVamosGraphs();
    if (!vamos_fragment.Identified()) goto mugast_label;

#ifdef VERBOSE_DEBUG
    std::cout << "------------>Finished: Vamos identification, positive exit\n";
#endif


    //AGATA///////////////////////////////////////////////////////////////////////////////////////////////////
    if (AGATA_GOOD) {
        for (long unsigned int ii = 0; ii < AddE.GetSize(); ii++) {
            if (AddE[ii] > 10) {
                Fill(pData.AGATA.hDC[vamos_fragment.Get_id_M()][vamos_fragment.Get_id_Z()]["NONE"], 1E3 * CorrectDoppler(*vamos_fragment.Get_p4(), AddE[ii] / 1E3, AddX[ii], AddY[ii], AddZ[ii]));
                for (long unsigned int kk = 0; kk < (*Mugast).PosX.size(); kk++) {
                    TVector3 vec((*Mugast).PosX[kk], (*Mugast).PosY[kk], (*Mugast).PosZ[kk]);
                    Fill(pData.AGATA.mDC_ThetaMUGAST[vamos_fragment.Get_id_M()][vamos_fragment.Get_id_Z()], 1E3 * CorrectDoppler(*vamos_fragment.Get_p4(), AddE[ii] / 1E3, AddX[ii], AddY[ii], AddZ[ii]), vec.Theta() * TMath::RadToDeg());
                }
                //Gamma Gamma matrices
                for (long unsigned int jj = 0; jj < AddE.GetSize(); jj++) {
                    if (AddE[jj] > 10 && ii != jj) {
                        Fill(pData.AGATA.mDC[vamos_fragment.Get_id_M()][vamos_fragment.Get_id_Z()], 1E3 * CorrectDoppler(*vamos_fragment.Get_p4(), AddE[ii] / 1E3, AddX[ii], AddY[ii], AddZ[ii]), 1E3 * CorrectDoppler(*vamos_fragment.Get_p4(), AddE[jj] / 1E3, AddX[jj], AddY[jj], AddZ[jj]));
                        Fill(pData.AGATA.mDC[vamos_fragment.Get_id_M()][vamos_fragment.Get_id_Z()], 1E3 * CorrectDoppler(*vamos_fragment.Get_p4(), AddE[jj] / 1E3, AddX[jj], AddY[jj], AddZ[jj]), 1E3 * CorrectDoppler(*vamos_fragment.Get_p4(), AddE[ii] / 1E3, AddX[ii], AddY[ii], AddZ[ii]));
                    }
                }
            }
        }
    }

//MUGAST//////////////////////////////////////////////////////////////////////////////////////////////////////
mugast_label:  //Label of goto previous to VAMOS

#ifdef VERBOSE_DEBUG
    std::cout << "------------>Finished: Vamos identification, negative exit\n";
#endif

    FillMugastConfHistograms();

    //SI data loops
    try {
        //E-TOF plots
        if (vamos_fragment.Identified()) {
            for (long unsigned int ii = 0; ii < (*Mugast).DSSD_E.size(); ii++) {
                //TVector3 hitPos((*Mugast).PosX[ii], (*Mugast).PosY[ii], (*Mugast).PosZ[ii]);
                //Fill(pData.SI.mE_TOF[vamos_fragment.Get_id_M()][vamos_fragment.Get_id_Z()][Form("MG%i", (*Mugast).TelescopeNumber[ii])], (*Mugast).DSSD_E[ii], AlignT((*Mugast).TelescopeNumber[ii], (*Mugast).DSSD_Y[ii], (*Mugast).DSSD_T[ii]));
                Fill(pData.SI.mE_TOF[vamos_fragment.Get_id_M()][vamos_fragment.Get_id_Z()][Form("MG%i", (*Mugast).TelescopeNumber[ii])], 
                        (*Mugast).DSSD_E[ii], (*Mugast).DSSD_T[ii]);
                Fill(pData.SI.mE_TOF[vamos_fragment.Get_id_M()][vamos_fragment.Get_id_Z()]["MG"], 
                        (*Mugast).DSSD_E[ii], AlignPunch((*Mugast).TelescopeNumber[ii], (*Mugast).DSSD_T[ii]));
            }
        }

        //Loop on physics data
        int MugastEvents = 0;
        for (long unsigned int ii = 0; ii < DetID.GetSize(); ii++) {
            Fill(pData.SI.mELab_ThetaLab[0][0]["ANY"]["ANY"], ThetaLab[MugastEvents], ELab[MugastEvents]);
            if (DetID[ii] >= 100)  //These events are in Must2
                continue;
            if (DetID[ii] != (*Mugast).TelescopeNumber[MugastEvents])
                std::cout << "Something is wrong matching Mugast events\n";
            //Excitation energy and kinematic lines

            Fill(pData.SI.hEx           [vamos_fragment.Get_id_M()][vamos_fragment.Get_id_Z()]["ANY"], Ex[MugastEvents]);
            Fill(pData.SI.mEx_TW        [vamos_fragment.Get_id_M()][vamos_fragment.Get_id_Z()]["ANY"], *TW, Ex[MugastEvents]);
            Fill(pData.SI.mELab_ThetaLab[vamos_fragment.Get_id_M()][vamos_fragment.Get_id_Z()]["ANY"]["ANY"], ThetaLab[MugastEvents], ELab[MugastEvents]);
            Fill(pData.SI.mECM_ThetaCM  [vamos_fragment.Get_id_M()][vamos_fragment.Get_id_Z()]["ANY"], ThetaCM[MugastEvents], Ecm[MugastEvents]);

            if (AGATA_GOOD) {
                //Loop over gammas
                for (long unsigned int ii = 0; ii < AddE.GetSize(); ii++) {
                    if (AddE[ii] > 10) {
                        Fill(pData.AGATA.mEx_DC     [vamos_fragment.Get_id_M()][vamos_fragment.Get_id_Z()]["ANY"], 
                                Ex[MugastEvents], CorrectDoppler(*vamos_fragment.Get_p4(), AddE[ii] / 1E3, AddX[ii], AddY[ii], AddZ[ii]));
                    }
                    if (AddE[ii] > 320 && AddE[ii] < 390) {
                        Fill(pData.SI.mELab_ThetaLab[vamos_fragment.Get_id_M()][vamos_fragment.Get_id_Z()]["ANY"]["360 keV"], 
                                ThetaLab[MugastEvents], ELab[MugastEvents]);
                    }
                }
            }
            //if (cut.at("MUGAST").at(Form("CUT_ETOF_MG%d", (*Mugast).TelescopeNumber[MugastEvents]))->IsInside((*Mugast).DSSD_E[MugastEvents], AlignT((*Mugast).TelescopeNumber[MugastEvents], (*Mugast).DSSD_Y[MugastEvents], (*Mugast).DSSD_T[MugastEvents])))
            for (const auto &particle : particles) {
                if (particle == "ANY")
                    continue;
                if (cut.at("MUGAST").at(Form("E_TOF_MG%d_%s", (*Mugast).TelescopeNumber[MugastEvents], particle.c_str()))->IsInside((*Mugast).DSSD_E[MugastEvents], (*Mugast).DSSD_T[MugastEvents])) {
                    if (AGATA_GOOD) {
                        //Loop over gammas
                        for (long unsigned int ii = 0; ii < AddE.GetSize(); ii++) {
                            if (AddE[ii] > 10) {
                                Fill(pData.AGATA.mEx_DC[vamos_fragment.Get_id_M()][vamos_fragment.Get_id_Z()][particle], Ex[MugastEvents], CorrectDoppler(*vamos_fragment.Get_p4(), AddE[ii] / 1E3, AddX[ii], AddY[ii], AddZ[ii]));
                            }
                            if (AddE[ii] > 320 && AddE[ii] < 390) {
                                if (particle == "2_H")
                                    //cout << "Filling gammas "
                                    //     << " id->GetMass(): " << IdentifiedNucleus->GetMass() << " i->GetNucleus(): " << IdentifiedNucleus->GetNucleus() << " particle: " << particle << " gamma: "
                                    //     << "360 keV"
                                    //     << " ELab: " << ThetaLab[MugastEvents] << " ELab: " << ELab[MugastEvents] << std::endl;
                                    Fill(pData.SI.mELab_ThetaLab[vamos_fragment.Get_id_M()][vamos_fragment.Get_id_Z()][particle]["360 keV"], ThetaLab[MugastEvents], ELab[MugastEvents]);
                            }
                        }
                    }
                    Fill(pData.SI.hEx           [vamos_fragment.Get_id_M()][vamos_fragment.Get_id_Z()][particle],           Ex[MugastEvents]);
                    Fill(pData.SI.mEx_TW        [vamos_fragment.Get_id_M()][vamos_fragment.Get_id_Z()][particle],           *TW, Ex[MugastEvents]);
                    Fill(pData.SI.mELab_ThetaLab[vamos_fragment.Get_id_M()][vamos_fragment.Get_id_Z()][particle]["ANY"],    ThetaLab[MugastEvents], ELab[MugastEvents]);
                    Fill(pData.SI.mECM_ThetaCM  [vamos_fragment.Get_id_M()][vamos_fragment.Get_id_Z()][particle],           ThetaCM[MugastEvents], Ecm[MugastEvents]);
                    Fill(pData.SI.mELab_ThetaLab[0][0][particle]["ANY"],                                                    ThetaLab[MugastEvents], ELab[MugastEvents]);
                }
            }
            MugastEvents++;
        }
        //Cats///////////////////////////////////////////////////////////////////////////////////////////////////
        for (long unsigned int ii = 0; ii < (*CATS).PositionX.size(); ii++) {
            Fill(pConf.CATS.mCATSpos, (*CATS).PositionX[ii], (*CATS).PositionY[ii]);
        }
    } catch (std::out_of_range &e) {
        std::cerr << "Silicon Physics loop :" << e.what() << std::endl;
    }

    //MUST2
    for (long unsigned int ii = 0; ii < (*MUST2).Si_E.size(); ii++) {
        Fill(pConf.SI.mE_TOF[Form("MM%i", (*MUST2).TelescopeNumber[ii])], (*MUST2).Si_E[ii], (*MUST2).Si_T[ii]);
        if (vamos_fragment.Identified()) {
            Fill(pData.SI.mE_TOF[vamos_fragment.Get_id_M()][vamos_fragment.Get_id_Z()][Form("MM%i",(*MUST2).TelescopeNumber[ii])],
                 (*MUST2).Si_E[ii], (*MUST2).Si_T[ii]);
        }
    }

    //AGATA
    // Fill(pConf.AGATA.hAddTS_LTS,*AddTS-*LTS);
    for (int ii = 0; ii < *nbAdd; ii++) {
        Fill(pConf.AGATA.mmAGATA3D, AddX[ii], AddY[ii], AddZ[ii]);
    }
    return kTRUE;
}

void Selector::SlaveTerminate() {
    // The SlaveTerminate() function is called after all entries or objects
    // have been processed. When running with PROOF SlaveTerminate() is called
    // on each slave server.
#ifdef VERBOSE_DEBUG
    std::cout << "------------>Selector::SlaveTerminate()\n";
#endif
    TFile *top = new TFile(file_name.c_str(), "recreate");
    std::cout << "Output file : " << file_name << "\n";
    TIter iter(fOutput);
    TObject *obj;
    while ((obj = iter())) {
        obj->Write();
    }
    top->Close();
}

void Selector::Terminate() {
    // The Terminate() function is the last function to be called during
    // a query. It always runs on the client, it can be used to present
    // the results graphically or save the results to file.
#ifdef VERBOSE_DEBUG
    std::cout << "------------>Selector::Terminate()\n";
#endif
}

//Identification of the VAMOS fragments
//inline void Selector::IdentifyFragment() {

inline void Selector::LoadVamosData() {
    vamos_fragment.SetData(new VamosIdentification::Data(&IC, &Path, &Brho, &Xf,
                                                         &ThetaL, &PhiL, &AGAVA_VAMOSTS,
                                                         &T_FPMW_CATS2_C));
}

inline void Selector::PlotVamosGraphs() {
    //dE-E plot, no conditions
    Fill(pConf.VAMOS.mdE_E, vamos_fragment.Get_En(), vamos_fragment.Get_D_En());
    Fill(pConf.VAMOS.mdE2_E, vamos_fragment.Get_En(), vamos_fragment.Get_D_En2());
    if (!vamos_fragment.Identified()) return;

    Fill(pConf.VAMOS.mQ_MQ[vamos_fragment.Get_id_Z()], vamos_fragment.Get_M_Q(), vamos_fragment.Get_Charge());
    Fill(pData.VAMOS.mTW_Brho[vamos_fragment.Get_id_M()][vamos_fragment.Get_id_Z()], *TW, *Brho);
    //Fill(pConf.VAMOS.hAmass[Nucleus_tmp], vamos_fragment.Get_M_Q * stoi(Qcut.substr(1, 2)));
}

inline void Selector::FillMugastConfHistograms() {
    for (long unsigned int ii = 0; ii < (*Mugast).DSSD_E.size(); ii++) {
        //(XE)
        Fill(pConf.SI.mStrip_E[Form("MG%i", (*Mugast).TelescopeNumber[ii])]["X"], (*Mugast).DSSD_X[ii], (*Mugast).DSSD_E[ii]);
        //(YE)
        Fill(pConf.SI.mStrip_E[Form("MG%i", (*Mugast).TelescopeNumber[ii])]["Y"], (*Mugast).DSSD_Y[ii], (*Mugast).DSSD_E[ii]);
        //(E TOF)
        Fill(pConf.SI.mE_TOF[Form("MG%i", (*Mugast).TelescopeNumber[ii])], (*Mugast).DSSD_T[ii]);
        //Fill(pConf.SI.mE_TOF[Form("MG%i", (*Mugast).TelescopeNumber[ii])], (*Mugast).DSSD_E[ii], AlignT((*Mugast).TelescopeNumber[ii], (*Mugast).DSSD_Y[ii], (*Mugast).DSSD_T[ii]));
    }
    for (long unsigned int ii = 0; ii < (*Mugast).DSSD_T.size(); ii++) {
        //(TE)
        Fill(pConf.SI.mStrip_T[Form("MG%i", (*Mugast).TelescopeNumber[ii])]["X"], (*Mugast).DSSD_X[ii], (*Mugast).DSSD_T[ii]);
        //(TE)
        Fill(pConf.SI.mStrip_T[Form("MG%i", (*Mugast).TelescopeNumber[ii])]["Y"], (*Mugast).DSSD_Y[ii], (*Mugast).DSSD_T[ii]);
    }
}

inline Double_t Selector::CorrectDoppler(const TLorentzVector &p4, const Double_t &Egamma,
                                         const Double_t &X, const Double_t &Y, const Double_t &Z) {
    TLorentzVector pgamma(Egamma, 0, 0, Egamma);
    TVector3 PosGamma(X, Y, Z + agata_Zshift);
    pgamma.SetPhi(PosGamma.Phi());
    pgamma.SetTheta(PosGamma.Theta());
    pgamma.SetE(Egamma);
    pgamma.Boost(-p4.BoostVector());
    return pgamma.Energy();
}

inline Double_t Selector::CorrectTOF(const Double_t &tof, const TVector3 &pos,
                                     const Double_t &Ek, const std::string &coeff) {
    //return ((tof)-2.99792* pos.Mag() *std::stof(coeff)*(Ek+mass["1H"])/(sqrt((Ek+mass["1H"])*(Ek+mass["1H"])-mass["1H"]*mass["1H"])));
    return (tof - std::stof(coeff)) / pos.Mag();
    //return tof - 3.*pos.Mag()*sqrt(mass["1H"]/(2*Ek));
}

bool Selector::GetSettings() {
    std::string config_file = "./Configs/GraphsEnabled.txt";
    std::ifstream ifile(config_file);
    if (ifile) {
        std::ifstream file(config_file);
        std::string Line;
        while (std::getline(file, Line)) {
            std::istringstream str(Line);
            std::string Graph;
            str >> Graph;
            std::string enabled;
            str >> enabled;
            if (!enabled.compare("false")) {
                TObject *to_delete = fOutput->FindObject(Graph.c_str());
                if (to_delete) {
                    fOutput->Remove(to_delete);
                    //std::cout<<"Deleted :"<< Graph <<std::endl;
                }
            }
        }
    } else {
        std::ofstream file(config_file);
        TIter iter(fOutput);
        TObject *obj;
        while ((obj = iter())) {
            file << obj->GetName() << "\t\t\t\t"
                 << "false" << std::endl;
        }
        file.close();
    }
    return true;
}

inline bool Selector::Fill(TH1 *histo, const Double_t &data1) {
    if (histo == nullptr || !fOutput->FindObject(histo)) {
        return false;
    } else {
        histo->Fill(data1);
        //  std::cout << "filling :"<<data1 <<std::endl;
    }
    return true;
}

inline bool Selector::Fill(TH2 *histo, const Double_t &data1, const Double_t &data2) {
    if (histo == nullptr || !fOutput->FindObject(histo)) {
        return false;
    } else {
        histo->Fill(data1, data2);
    }
    return true;
}

inline bool Selector::Fill(TH3 *histo, const Double_t &data1,
                           const Double_t &data2, const Double_t &data3) {
    if (histo == nullptr || !fOutput->FindObject(histo)) {
        return false;
    } else {
        histo->Fill(data1, data2, data3);
    }
    return true;
}