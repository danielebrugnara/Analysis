#define Selector_cxx

#include "Selector.h"

Selector::Selector(TTree *){}

Selector::~Selector(){
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

    TString option = GetOption();
}

void Selector::SlaveBegin(TTree * /*tree*/) {
    // The SlaveBegin() function is called after the Begin() function.
    // When running with PROOF SlaveBegin() is called on each slave server.
    // The tree argument is deprecated (on PROOF 0 is passed).

    TString option = GetOption();
    file_name = "./Out/" + option;

    //Initializing Histograms

    //Configurations
    //Target
    thickness_angle = new Interpolation("./Configs/Interpolations/GasThickness.txt");
    angle_angle     = new Interpolation("./Configs/Interpolations/EntranceAngleHavar.txt");

    //////////////////////////////////////////////////////////////////////////////////////////////////
    ///Histogram pointer initialization//////////////////////////////////////////////////////////////

    //Config histograms//////////////////////////////////////////////////////////////////////////////
    //VAMOS
    pConf.VAMOS.mdE_E = new TH2D("pConf-VAMOS-mdE_E", "dE E in VAMOS", 4000, 10, 350, 4000, 10, 140);
    fOutput->Add(pConf.VAMOS.mdE_E);

    pConf.VAMOS.mdE2_E = new TH2D("pConf-VAMOS-mdE2_E", "dE2 E in VAMOS", 2000, 10, 350, 2000, 10, 140);
    fOutput->Add(pConf.VAMOS.mdE2_E);

    pConf.VAMOS.mQ_MQ["Ar"] = new TH2D("pConf-VAMOS-mQ_MQ-Ar", "M/Q vs Q with Ar selection", 1000, 2, 4, 1000, 3, 24);
    fOutput->Add(pConf.VAMOS.mQ_MQ["Ar"]);
    pConf.VAMOS.mQ_MQ["K"] = new TH2D("pConf-VAMOS-mQ_MQ-K", "M/Q vs Q with K selection", 1000, 2, 4, 1000, 3, 24);
    fOutput->Add(pConf.VAMOS.mQ_MQ["K"]);

    pConf.VAMOS.hAmass["Ar"] = new TH1D("pConf-VAMOS-hAmass-Ar", "Mass histogram of Ar", 1000, 20, 50);
    fOutput->Add(pConf.VAMOS.hAmass["Ar"]);
    pConf.VAMOS.hAmass["K"] = new TH1D("pConf-VAMOS-hAmass-K", "Mass histogram of K", 1000, 20, 50);
    fOutput->Add(pConf.VAMOS.hAmass["K"]);

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
    for (const auto &nucleus : nuclei) {
        for (const auto &mass : masses) {
            pData.VAMOS.mTW_Brho[mass][nucleus] = new TH2D(Form("pData-VAMOS-mTW_Brho-%s-%s", mass.c_str(), nucleus.c_str()), Form("Time vs Brho with %s %s in VAMOS", mass.c_str(), nucleus.c_str()), 5000, 242, 328, 1000, 0.5, 1.5);
            fOutput->Add(pData.VAMOS.mTW_Brho[mass][nucleus]);
        }
    }
    //VAMOS

    //AGATA
    for (const auto &nucleus : nuclei) {
        for (const auto &mass : masses) {
            for (const auto &condition : AGATAconditions) {
                pData.AGATA.hDC[mass][nucleus][condition] = new TH1D(Form("pData-AGATA-hDC-%s-%s-cond%s", mass.c_str(), nucleus.c_str(), condition.c_str()), Form("DC spectrum of %s %s with condition %s", mass.c_str(), nucleus.c_str(), condition.c_str()), 4000, 0, 4000);
                fOutput->Add(pData.AGATA.hDC[mass][nucleus][condition]);
            }
            pData.AGATA.mDC[mass][nucleus] = new TH2D(Form("pData-AGATA-mDC-%s-%s", mass.c_str(), nucleus.c_str()), Form("DC gamma gamma of %s %s", mass.c_str(), nucleus.c_str()), 4000, 0, 4000, 4000, 0, 4000);
            fOutput->Add(pData.AGATA.mDC[mass][nucleus]);
            pData.AGATA.mDC_ThetaMUGAST[mass][nucleus] = new TH2D(Form("pData-AGATA-mDC_ThetaMUGAST-%s-%s", mass.c_str(), nucleus.c_str()), Form("DC gamma vs Theta on MUGAST of %s %s", mass.c_str(), nucleus.c_str()), 4000, 0, 4000, 180, 0, 180);
            fOutput->Add(pData.AGATA.mDC_ThetaMUGAST[mass][nucleus]);
            for (const auto &particle : particles) {
                pData.AGATA.mEx_DC[mass][nucleus][particle] = new TH2D(Form("pData-AGATA-mEx_DC-%s-%s-%s", mass.c_str(), nucleus.c_str(), particle.c_str()), Form("Excitation energy AGATA vs MUGAST %s %s and %s", mass.c_str(), nucleus.c_str(), particle.c_str()), 1000, 0, 10, 1000, 0, 10);
                fOutput->Add(pData.AGATA.mEx_DC[mass][nucleus][particle]);
                pData.AGATA.mELab_ThetaLab[mass][nucleus][particle] = new TH2D(Form("pData-AGATA-mELab_ThetaLab-%s-%s-%s", mass.c_str(), nucleus.c_str(), particle.c_str()), Form("Excitation energy AGATA vs MUGAST %s %s and %s", mass.c_str(), nucleus.c_str(), particle.c_str()), 1000, 0, 10, 1000, 0, 10);
                fOutput->Add(pData.AGATA.mELab_ThetaLab[mass][nucleus][particle]);
            }
        }
    }

    //CATS

    //MUGAST
    for (const auto &nucleus : nuclei) {
        for (const auto &mass : masses) {
            for (const auto &MM : siliconsMM) {
                pData.SI.mdE_E_Si[MM][mass][nucleus] = new TH2D(Form("pData-SI-mdE_E_Si-%s-%s-%s", MM.c_str(), mass.c_str(), nucleus.c_str()), Form("dE E of %s with %s %s", MM.c_str(), mass.c_str(), nucleus.c_str()), 1000, 0, 28, 1000, 0, 28);
                fOutput->Add(pData.SI.mdE_E_Si[MM][mass][nucleus]);
            }
            for (const auto &SI : silicons) {
                pData.SI.mE_TOF[SI][mass][nucleus] = new TH2D(Form("pData-SI-mE_TOF-%s-%s-%s", SI.c_str(), mass.c_str(), nucleus.c_str()), Form("E vs TOF of %s with %s %s", SI.c_str(), mass.c_str(), nucleus.c_str()), 1000, 0, 28, 1000, 260, 380);
                fOutput->Add(pData.SI.mE_TOF[SI][mass][nucleus]);
            }
            pData.SI.mE_TOF["MG"][mass][nucleus] = new TH2D(Form("pData-SI-mE_TOF-MG-%s-%s", mass.c_str(), nucleus.c_str()), Form("E vs TOF of all MG with %s %s", mass.c_str(), nucleus.c_str()), 1000, 0, 28, 1000, 260, 380);
            fOutput->Add(pData.SI.mE_TOF["MG"][mass][nucleus]);
            for (const auto &particle : particles) {
                pData.SI.hEx[mass][nucleus][particle] = new TH1D(Form("pData-SI-hEx-%s-%s-%s", mass.c_str(), nucleus.c_str(), particle.c_str()), Form("Excitation energy with %s %s in VAMOS and %s in MUGAST", mass.c_str(), nucleus.c_str(), particle.c_str()), 1000, -60, 60);
                fOutput->Add(pData.SI.hEx[mass][nucleus][particle]);
                pData.SI.mEx_TW[mass][nucleus][particle] = new TH2D(Form("pData-SI-mEx_TW-%s-%s-%s", mass.c_str(), nucleus.c_str(), particle.c_str()), Form("Excitation energy vs Time with %s %s in VAMOS and %s in MUGAST", mass.c_str(), nucleus.c_str(), particle.c_str()), 5000, 242, 328, 1000, -60, 60);
                fOutput->Add(pData.SI.mEx_TW[mass][nucleus][particle]);
                pData.SI.mECM_ThetaCM[mass][nucleus][particle] = new TH2D(Form("pData-SI-mECM_ThetaCM-%s-%s-%s", mass.c_str(), nucleus.c_str(), particle.c_str()), Form("E CM vs Theta CM with %s %s in VAMOS and %s in MUGAST", mass.c_str(), nucleus.c_str(), particle.c_str()), 1000, 0, 180, 1000, 0, 60);
                fOutput->Add(pData.SI.mECM_ThetaCM[mass][nucleus][particle]);
                for (const auto &gamma : gammas) {
                    pData.SI.mELab_ThetaLab[mass][nucleus][particle][gamma] = new TH2D(Form("pData-SI-mELab_ThetaLab-%s-%s-%s-%s", mass.c_str(), nucleus.c_str(), particle.c_str(), gamma.c_str()), Form("ELab vs Theta Lab with %s %s in VAMOS and %s in MUGAST and %s in AGATA", mass.c_str(), nucleus.c_str(), particle.c_str(), gamma.c_str()), 1000, 0, 180, 1000, 0, 60);
                    fOutput->Add(pData.SI.mELab_ThetaLab[mass][nucleus][particle][gamma]);
                }
            }
        }
    }
    for (const auto &particle : particles) {
        pData.SI.mELab_ThetaLab["ANY"]["ANY"][particle]["ANY"] = new TH2D(Form("pData-SI-mELab_ThetaLab-%s-%s-%s", "ANY", "ANY", particle.c_str()), Form("E Lab vs Theta Lab with %s %s in VAMOS and %s in MUGAST", "ANY", "ANY", particle.c_str()), 1000, 0, 180, 1000, 0, 60);
        fOutput->Add(pData.SI.mELab_ThetaLab["ANY"]["ANY"][particle]["ANY"]);
    }

    ///////////////////////////////////////////////////////////////////////////////////////////
    ///Loading things from TFiles/////////////////////////////////////////////////////////////
    std::ifstream *input_file = nullptr;  

    //GCuts///////////////////////////////////////////////////////////////////////////////////
    std::string VAMOS_cuts_file = "./Configs/Cuts/VAMOS.root";
    std::string MUGAST_cuts_file = "./Configs/Cuts/MUGAST.root";

    input_file = new std::ifstream(VAMOS_cuts_file);
    if (input_file) {
        input_file->close();
        VamosFragment.LoadCuts(VAMOS_cuts_file);
        /*
        TFile *VAMOScuts = new TFile(VAMOS_cuts_file.c_str(), "READ");
        if (!(VAMOScuts->IsOpen())) {
            std::cout << "VAMOS file not opened\n";
        } else {
            TIter contents(VAMOScuts->GetListOfKeys());
            TKey *key;
            TObject *obj;
            while ((key = (TKey *)contents())) {
                obj = VAMOScuts->Get(key->GetName());
                if (obj->InheritsFrom("TCutG")) {
                    TCutG *tmp = (TCutG *)obj;
                    cut["VAMOS"][tmp->GetName()] = tmp;
                    //std::cout << "Found cut in VAMOS :" << tmp->GetName() << std::endl;
                }
            }
        }
        */
    } else {
        std::cout << "VAMUS cuts not found\n";
    }

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
    CheckCutsPresence();

    //Interpolations///////////////////////////////////////////////////////////////////////////////

    //FP-time
    std::string FP_interpolation_file = "./Configs/Interpolations/InterpolationTimeFP.root";
    input_file = new std::ifstream(FP_interpolation_file);
    if (input_file) {
        input_file->close();
        TFile *tmp = new TFile(FP_interpolation_file.c_str(), "READ");
        if (!(tmp->IsOpen())) throw std::runtime_error("Interpolation Time FP ROOT file not opened\n"); 
        FP_time_interpolation = new Interpolation(tmp);
    } else {
        std::cout << "FP interpolation file not found\n";
    }


    //////////////////////////////////////////////////////////////////////////////////////////////////
    ///Other configs/////////////////////////////////////////////////////////////////////////////////
    mass["46_Ar"] = 45.968082712 * AMU_TO_MEV;  //in MeV
    mass["47_Ar"] = 46.972934865 * AMU_TO_MEV;  //in MeV
    mass["47_K"] = 46.961661614 * AMU_TO_MEV;   //in MeV
    mass["46_K"] = 45.961981586 * AMU_TO_MEV;   //In MeV
    mass["2_H"] = 2.01410177812 * AMU_TO_MEV;   //In MeV
    mass["1_H"] = 1.00782503223 * AMU_TO_MEV;   //In MeV


    //////////////////////////////////////////////////////////////////////////////////////////////////
    ///Temporary histogram///////////////////////////////////////////////////////////////////////////
    general_histo_ptr = new TH2D("test", "test", 1000, -400, 400, 1000, 0, 100);
    fOutput->Add(general_histo_ptr);


    //////////////////////////////////////////////////////////////////////////////////////////////////
    ///Deciding which graphs to plot/////////////////////////////////////////////////////////////////
    GetSettings();  //Decides which histograms to fill
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

    ++counter.General;

    //Agata vs Ancillary coincidence gate
    Bool_t AGATA_GOOD = *AddTS - *LTS > 175 && *AddTS - *LTS < 184;

    //Vamos identification////////////////////////////////////////////////////////////////////////////////////
    //Configuration spectra are filled in the identification
    if (!IdentifiedNucleus)
        delete IdentifiedNucleus;
    IdentifiedNucleus = new VamosId();
    if (IC.GetSize() == 0 || !(MW_N.GetSize() >  0)) {  //VAMOS is not OK, only look at the SI detectors
    //if (IC.GetSize() == 0 ) {  //VAMOS is not OK, only look at the SI detectors
        goto mugast_label;
    }
    IdentifyFragment();

    if (IdentifiedNucleus->Identified) {
        Fill(pData.VAMOS.mTW_Brho[IdentifiedNucleus->GetMass()][IdentifiedNucleus->GetNucleus()], *TW, *Brho);
    }

    //AGATA///////////////////////////////////////////////////////////////////////////////////////////////////
    if (AGATA_GOOD) {
        for (long unsigned int ii = 0; ii < AddE.GetSize(); ii++) {
            if (AddE[ii] > 10) {
                Fill(pData.AGATA.hDC[IdentifiedNucleus->GetMass()][IdentifiedNucleus->GetNucleus()]["NONE"], 1E3 * CorrectDoppler(IdentifiedNucleus->p4, AddE[ii] / 1E3, AddX[ii], AddY[ii], AddZ[ii]));
                for (long unsigned int kk = 0; kk < (*Mugast).PosX.size(); kk++) {
                    TVector3 vec((*Mugast).PosX[kk], (*Mugast).PosY[kk], (*Mugast).PosZ[kk]);
                    Fill(pData.AGATA.mDC_ThetaMUGAST[IdentifiedNucleus->GetMass()][IdentifiedNucleus->GetNucleus()], 1E3 * CorrectDoppler(IdentifiedNucleus->p4, AddE[ii] / 1E3, AddX[ii], AddY[ii], AddZ[ii]), vec.Theta() * TMath::RadToDeg());
                }
                //Gamma Gamma matrices
                for (long unsigned int jj = 0; jj < AddE.GetSize(); jj++) {
                    if (AddE[jj] > 10 && ii != jj) {
                        Fill(pData.AGATA.mDC[IdentifiedNucleus->GetMass()][IdentifiedNucleus->GetNucleus()], 1E3 * CorrectDoppler(IdentifiedNucleus->p4, AddE[ii] / 1E3, AddX[ii], AddY[ii], AddZ[ii]), 1E3 * CorrectDoppler(IdentifiedNucleus->p4, AddE[jj] / 1E3, AddX[jj], AddY[jj], AddZ[jj]));
                        Fill(pData.AGATA.mDC[IdentifiedNucleus->GetMass()][IdentifiedNucleus->GetNucleus()], 1E3 * CorrectDoppler(IdentifiedNucleus->p4, AddE[jj] / 1E3, AddX[jj], AddY[jj], AddZ[jj]), 1E3 * CorrectDoppler(IdentifiedNucleus->p4, AddE[ii] / 1E3, AddX[ii], AddY[ii], AddZ[ii]));
                    }
                }
            }
        }
    }

//MUGAST//////////////////////////////////////////////////////////////////////////////////////////////////////
mugast_label:  //Label of goto previous to VAMOS
    FillMugastConfHistograms();

    //SI data loops
    try {
        //E-TOF plots
        if (IdentifiedNucleus->Identified) {
            for (long unsigned int ii = 0; ii < (*Mugast).DSSD_E.size(); ii++) {
                //TVector3 hitPos((*Mugast).PosX[ii], (*Mugast).PosY[ii], (*Mugast).PosZ[ii]);
                Fill(pData.SI.mE_TOF[Form("MG%i", (*Mugast).TelescopeNumber[ii])][IdentifiedNucleus->GetMass()][IdentifiedNucleus->GetNucleus()], (*Mugast).DSSD_E[ii], AlignT((*Mugast).TelescopeNumber[ii], (*Mugast).DSSD_Y[ii], (*Mugast).DSSD_T[ii]));
                Fill(pData.SI.mE_TOF["MG"][IdentifiedNucleus->GetMass()][IdentifiedNucleus->GetNucleus()], (*Mugast).DSSD_E[ii], AlignPunch((*Mugast).TelescopeNumber[ii], (*Mugast).DSSD_T[ii]));
            }
        }

        //Loop on physics data
        int MugastEvents = 0;
        for (long unsigned int ii = 0; ii < DetID.GetSize(); ii++) {
            Fill(pData.SI.mELab_ThetaLab["ANY"]["ANY"]["ANY"]["ANY"], ThetaLab[MugastEvents], ELab[MugastEvents]);
            if (DetID[ii] >= 100)  //These events are in Must2
                continue;
            if (DetID[ii] != (*Mugast).TelescopeNumber[MugastEvents])
                std::cout << "Something is wrong matching Mugast events\n";
            //Excitation energy and kinematic lines

            Fill(pData.SI.hEx[IdentifiedNucleus->GetMass()][IdentifiedNucleus->GetNucleus()]["ANY"], Ex[MugastEvents]);
            Fill(pData.SI.mEx_TW[IdentifiedNucleus->GetMass()][IdentifiedNucleus->GetNucleus()]["ANY"], *TW, Ex[MugastEvents]);
            Fill(pData.SI.mELab_ThetaLab[IdentifiedNucleus->GetMass()][IdentifiedNucleus->GetNucleus()]["ANY"]["ANY"], ThetaLab[MugastEvents], ELab[MugastEvents]);
            Fill(pData.SI.mECM_ThetaCM[IdentifiedNucleus->GetMass()][IdentifiedNucleus->GetNucleus()]["ANY"], ThetaCM[MugastEvents], Ecm[MugastEvents]);

            if (AGATA_GOOD) {
                //Loop over gammas
                for (long unsigned int ii = 0; ii < AddE.GetSize(); ii++) {
                    if (AddE[ii] > 10) {
                        Fill(pData.AGATA.mEx_DC[IdentifiedNucleus->GetMass()][IdentifiedNucleus->GetNucleus()]["ANY"], Ex[MugastEvents], CorrectDoppler(IdentifiedNucleus->p4, AddE[ii] / 1E3, AddX[ii], AddY[ii], AddZ[ii]));
                    }
                    if (AddE[ii] > 320 && AddE[ii] < 390) {
                        Fill(pData.SI.mELab_ThetaLab[IdentifiedNucleus->GetMass()][IdentifiedNucleus->GetNucleus()]["ANY"]["360 keV"], ThetaLab[MugastEvents], ELab[MugastEvents]);
                    }
                }
            }
            //if (cut.at("MUGAST").at(Form("CUT_ETOF_MG%d", (*Mugast).TelescopeNumber[MugastEvents]))->IsInside((*Mugast).DSSD_E[MugastEvents], AlignT((*Mugast).TelescopeNumber[MugastEvents], (*Mugast).DSSD_Y[MugastEvents], (*Mugast).DSSD_T[MugastEvents])))
            for (const auto &particle : particles) {
                if (particle == "ANY")
                    continue;
                if (cut.at("MUGAST").at(Form("E_TOF_MG%d_%s", (*Mugast).TelescopeNumber[MugastEvents], particle.c_str()))->IsInside((*Mugast).DSSD_E[MugastEvents], AlignT((*Mugast).TelescopeNumber[MugastEvents], (*Mugast).DSSD_Y[MugastEvents], (*Mugast).DSSD_T[MugastEvents]))) {
                    if (AGATA_GOOD) {
                        //Loop over gammas
                        for (long unsigned int ii = 0; ii < AddE.GetSize(); ii++) {
                            if (AddE[ii] > 10) {
                                Fill(pData.AGATA.mEx_DC[IdentifiedNucleus->GetMass()][IdentifiedNucleus->GetNucleus()][particle], Ex[MugastEvents], CorrectDoppler(IdentifiedNucleus->p4, AddE[ii] / 1E3, AddX[ii], AddY[ii], AddZ[ii]));
                            }
                            if (AddE[ii] > 320 && AddE[ii] < 390) {
                                if (particle == "2_H")
                                    cout << "Filling gammas "
                                         << " id->GetMass(): " << IdentifiedNucleus->GetMass() << " i->GetNucleus(): " << IdentifiedNucleus->GetNucleus() << " particle: " << particle << " gamma: "
                                         << "360 keV"
                                         << " ELab: " << ThetaLab[MugastEvents] << " ELab: " << ELab[MugastEvents] << std::endl;
                                Fill(pData.SI.mELab_ThetaLab[IdentifiedNucleus->GetMass()][IdentifiedNucleus->GetNucleus()][particle]["360 keV"], ThetaLab[MugastEvents], ELab[MugastEvents]);
                            }
                        }
                    }
                    Fill(pData.SI.hEx[IdentifiedNucleus->GetMass()][IdentifiedNucleus->GetNucleus()][particle], Ex[MugastEvents]);
                    Fill(pData.SI.mEx_TW[IdentifiedNucleus->GetMass()][IdentifiedNucleus->GetNucleus()][particle], *TW, Ex[MugastEvents]);
                    Fill(pData.SI.mELab_ThetaLab[IdentifiedNucleus->GetMass()][IdentifiedNucleus->GetNucleus()][particle]["ANY"], ThetaLab[MugastEvents], ELab[MugastEvents]);
                    Fill(pData.SI.mELab_ThetaLab["ANY"]["ANY"][particle]["ANY"], ThetaLab[MugastEvents], ELab[MugastEvents]);
                    Fill(pData.SI.mECM_ThetaCM[IdentifiedNucleus->GetMass()][IdentifiedNucleus->GetNucleus()][particle], ThetaCM[MugastEvents], Ecm[MugastEvents]);
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
        if (IdentifiedNucleus->Identified) {
            Fill(pData.SI.mE_TOF[Form("MM%i", (*MUST2).TelescopeNumber[ii])][IdentifiedNucleus->GetMass()][IdentifiedNucleus->GetNucleus()], (*MUST2).Si_E[ii], (*MUST2).Si_T[ii]);
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
}

//Identification of the VAMOS fragments
//inline void Selector::IdentifyFragment() {

inline void Selector::LoadVamosData(){
    delete VamosFragment.data;
    VamosFragment.data = new VamosIdentification::Data(IC, *Path, *Brho, *Xf, *AGAVA_VAMOSTS, *T_FPMW_CATS2_C);
}

inline void Selector::PlotVamosGraphs() {
    //dE-E plot, no conditions
    Fill(pConf.VAMOS.mdE_E, VamosFragment.Get_En(), VamosFragment.Get_D_En());
    Fill(pConf.VAMOS.mdE2_E, VamosFragment.Get_En(), VamosFragment.Get_D_En2());
    if (VamosFragment.Get_En() > 70 &&VamosFragment.Get_En() < 95 &&VamosFragment.Get_D_En2() > 30 &&VamosFragment.Get_D_En2() < 45  )
    {
       general_histo_ptr->Fill(*Xf, VamosFragment.Get_En);
    }
    for (const auto &Zcut : Zcuts) {
        //if (cut.at("VAMOS").at(Zcut)->IsInside(VamosFragment.Get_En, VamosFragment.Get_D_En)) {
        if (cut.at("VAMOS").at(Zcut)->IsInside(VamosFragment.Get_En(), VamosFragment.Get_D_En2())) {
            std::string Nucleus_tmp = Zcut.substr(Zcut.find_last_of("_") + 1);
            Fill(pConf.VAMOS.mQ_MQ[Nucleus_tmp], VamosFragment.Get_M_Q(), VamosFragment.Get_Charge());
            for (const auto &Qcut : Qcuts) {
                if (cut.at("VAMOS").at(Qcut)->IsInside(VamosFragment.Get_M_Q(), VamosFragment.Get_Charge())) {
                    Fill(pConf.VAMOS.hAmass[Nucleus_tmp], VamosFragment.Get_M_Q * stoi(Qcut.substr(1, 2)));
                    for (const auto &Mgate : Mgates) {
                        if (VamosFragment.Get_M_Q * stoi(Qcut.substr(1, 2)) < stod(Mgate) + 0.5 && VamosFragment.Get_M_Q * stoi(Qcut.substr(1, 2)) > stod(Mgate) - 0.5) {
                            if (VamosFragment.Get_Nucleus().compare("") == 0) throw std::runtime_error("Fragment already identified, overlapping gates??");
                            VamosFragment.Get_Nucleus() = Mgate + "_" + Nucleus_tmp;  //Format is 46_Ar or 47_K
                            VamosFragment.Get_p4().SetT(mass[VamosFragment.Get_Nucleus()]);
                            TVector3 b4(0, 0, VamosFragment.Get_Beta());
                            b4.SetMagThetaPhi(VamosFragment.Get_Beta(), *ThetaL, *PhiL);
                            VamosFragment.Get_p4().Boost(b4);
                            if (VamosFragment.Get_Nucleus() == "46_Ar") ++counter.Ar46;
                            if (VamosFragment.Get_Nucleus() == "47_K") ++counter.K47;
                        }
                    }
                }
            }
        }
    }
}

inline double Selector::ComputeFPTimeCorrection(const double & FP_X){
    return 1.;
}

inline void Selector::FillMugastConfHistograms() {
    for (long unsigned int ii = 0; ii < (*Mugast).DSSD_E.size(); ii++) {
        //(XE)
        Fill(pConf.SI.mStrip_E[Form("MG%i", (*Mugast).TelescopeNumber[ii])]["X"], (*Mugast).DSSD_X[ii], (*Mugast).DSSD_E[ii]);
        //(YE)
        Fill(pConf.SI.mStrip_E[Form("MG%i", (*Mugast).TelescopeNumber[ii])]["Y"], (*Mugast).DSSD_Y[ii], (*Mugast).DSSD_E[ii]);
        //(E TOF)
        Fill(pConf.SI.mE_TOF[Form("MG%i", (*Mugast).TelescopeNumber[ii])], (*Mugast).DSSD_E[ii], AlignT((*Mugast).TelescopeNumber[ii], (*Mugast).DSSD_Y[ii], (*Mugast).DSSD_T[ii]));
    }
    for (long unsigned int ii = 0; ii < (*Mugast).DSSD_T.size(); ii++) {
        //(TE)
        Fill(pConf.SI.mStrip_T[Form("MG%i", (*Mugast).TelescopeNumber[ii])]["X"], (*Mugast).DSSD_X[ii], (*Mugast).DSSD_T[ii]);
        //(TE)
        Fill(pConf.SI.mStrip_T[Form("MG%i", (*Mugast).TelescopeNumber[ii])]["Y"], (*Mugast).DSSD_Y[ii], (*Mugast).DSSD_T[ii]);
    }
}

inline Double_t Selector::CorrectDoppler(const TLorentzVector &p4, const Double_t &Egamma, const Double_t &X, const Double_t &Y, const Double_t &Z) {
    TLorentzVector pgamma(Egamma, 0, 0, Egamma);
    TVector3 PosGamma(X, Y, Z + agata_Zshift);
    pgamma.SetPhi(PosGamma.Phi());
    pgamma.SetTheta(PosGamma.Theta());
    pgamma.SetE(Egamma);
    pgamma.Boost(-p4.BoostVector());
    return pgamma.Energy();
}

inline Double_t Selector::CorrectTOF(const Double_t &tof, const TVector3 &pos, const Double_t &Ek, const std::string &coeff) {
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

void Selector::CheckCutsPresence() {
    auto *tmp_ptr = &cut["VAMOS"];
    for (const auto &Zcut : Zcuts) {
        if (tmp_ptr->find(Zcut) == tmp_ptr->end()) {
            throw std::runtime_error("Unable to find cut :" + Zcut + "\n");
        }
    }
    for (const auto &Qcut : Qcuts) {
        if (tmp_ptr->find(Qcut) == tmp_ptr->end()) {
            throw std::runtime_error("Unable to find cut :" + Qcut + "\n");
        }
    }
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

inline bool Selector::Fill(TH3 *histo, const Double_t &data1, const Double_t &data2, const Double_t &data3) {
    if (histo == nullptr || !fOutput->FindObject(histo)) {
        return false;
    } else {
        histo->Fill(data1, data2, data3);
    }
    return true;
}