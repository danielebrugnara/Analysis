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
#ifdef VERBOSE_DEBUG
    std::cout << "------------>Selector::Begin()\n";
#endif

    TString option = GetOption();
}

void Selector::SlaveBegin(TTree * /*tree*/) {
#ifdef VERBOSE_DEBUG
    std::cout << "------------>Selector::SlaveBegin()\n";
#endif
    total_entries = fReader.GetEntries(false);

    TString option = GetOption();
    file_name = "./Out/" + option;

    //////////////////////////////////////////////////////////////////////////////////////////////////
    ///Deciding which graphs to plot/////////////////////////////////////////////////////////////////
    GetSettings();  //Decides which histograms to fill
    ///Loading Identification classes/////////////////////////////////////////////////////////////

    //TODO: use constructor to initialize this!!!
    std::string VAMOS_cuts_file = "./Configs/Cuts/VAMOS.root";
    vamos_fragment.LoadCuts(VAMOS_cuts_file);
    vamos_fragment.Initialize();

    std::string MUGAST_cuts_file = "./Configs/Cuts/MUGAST.root";
    mugast_fragment.LoadCuts(MUGAST_cuts_file);

    //Passing beam energy in MeV, target position mm
    //mugast_fragment.Initialize(379.04, TVector3(0, 0, 25.));
    mugast_fragment.Initialize(453.22, TVector3(0, 0, 25.));

#ifdef VERBOSE_DEBUG
    std::cout << "------------>finished: vamos_fragment initialization\n";
#endif

    //Initializing Histograms

    //Configurations
    //Target
    //thickness_angle = new Interpolation("./Configs/Interpolations/GasThickness.txt");
    //angle_angle = new Interpolation("./Configs/Interpolations/EntranceAngleHavar.txt");

    //////////////////////////////////////////////////////////////////////////////////////////////////
    ///Histogram pointer initialization//////////////////////////////////////////////////////////////

    //Config histograms//////////////////////////////////////////////////////////////////////////////
    //VAMOS
    Istantiate(pConf.VAMOS.mdE_E,
               new TH2D("pConf_VAMOS_mdE_E",
                        "dE E in VAMOS",
                        4000, 10, 350, 4000, 10, 140));

    Istantiate(pConf.VAMOS.mdE2_E,
               new TH2D("pConf_VAMOS_mdE2_E",
                        "dE2 E in VAMOS",
                        2000, 10, 350, 2000, 10, 140));

    for (const auto &Z_it : vamos_fragment.cuts_Z) {
        Istantiate(pConf.VAMOS.mQ_MQ[Z_it],
                   new TH2D(Form("pConf_VAMOS_mQ_MQ_Z%i", Z_it),
                            Form("M/Q vs Q with Z%i selection", Z_it),
                            1000, 2, 4, 1000, 3, 24));

        Istantiate(pConf.VAMOS.Xf_MQ[Z_it],
                   new TH2D(Form("pConf_VAMOS_Xf_MQ_Z%i", Z_it),
                            Form("M/Q vs Xf with Z%i selection", Z_it),
                            1000, 2, 4, 1000, -400, 400));
    }

    //AGATA
    Istantiate(pConf.AGATA.mmAGATA3D,
               new TH3D("pConf-AGATA-mmAGATA3D",
                        "Hit patter on AGATA",
                        50, -300, 300, 50, -300, 300, 50, -300, 300));

    Istantiate(pConf.AGATA.hAddTS_LTS,
               new TH1D("pConf-AGATA-hAddTS_LTS",
                        "Difference between AddTS and LTS",
                        1000, 0, 300));

    //CATS
    Istantiate(pConf.CATS.mCATSpos,
               new TH2D("pConf-CATS-mCATSpos",
                        "Beam position on CATS",
                        1000, -50, 50, 1000, -50, 50));

    //MUGAST
    for (int ii = 0; ii < 128; ii++) {
        strips.push_back(to_string(ii));
    }

    //    for (const auto &MM : siliconsMM) {
    //        Istantiate(pConf.SI.mdE_E_Si[MM],
    //                    new TH2D(Form("pConf-SI-mdE_E_Si-%s", MM.c_str()),
    //                                Form("dE E of %s", MM.c_str()),
    //                                1000, 0, 28, 1000, 0, 28));
    //    }
    //    for (const auto &SI : silicons) {
    //
    //        Istantiate(pConf.SI.mStrip_E[SI]["X"],
    //                    new TH2D(Form("pConf-SI-mStrip_E-%s-%s", SI.c_str(), "X"),
    //                                Form("Strip vs E of %s %s", SI.c_str(), "X"),
    //                                128, 0, 128, 1000, 0, 30));
    //
    //        Istantiate(pConf.SI.mStrip_E[SI]["Y"],
    //                    new TH2D(Form("pConf-SI-mStrip_E-%s-%s", SI.c_str(), "Y"),
    //                                Form("Strip vs E of %s %s", SI.c_str(), "Y"),
    //                                128, 0, 128, 1000, 0, 30));
    //
    //        Istantiate(pConf.SI.mStrip_T[SI]["X"],
    //                    new TH2D(Form("pConf-SI-mStrip_T-%s-%s", SI.c_str(), "X"),
    //                                Form("Strip vs T of %s %s", SI.c_str(), "X"),
    //                                128, 0, 128, 1000, 0, 1500));
    //
    //        Istantiate(pConf.SI.mStrip_T[SI]["Y"],
    //                    new TH2D(Form("pConf-SI-mStrip_T-%s-%s", SI.c_str(), "Y"),
    //                                Form("Strip vs T of %s %s", SI.c_str(), "Y"),
    //                                128, 0, 128, 1000, 0, 1500));
    //    }

    //Data histograms//////////////////////////////////////////////////////////////////////////////
    for (const auto &it_M : vamos_fragment.cuts_M) {
        for (const auto &it_Z : vamos_fragment.cuts_Z) {
            Istantiate(pData.VAMOS.mTW_Brho[it_M][it_Z],
                       new TH2D(Form("pData_VAMOS_mTW_Brho_%i_%i", it_M, it_Z),
                                Form("Time vs Brho with M%i Z%i in VAMOS", it_M, it_Z),
                                5000, 242, 328, 1000, 0.5, 1.5));
        }
    }
    //VAMOS

    //AGATA
    for (const auto &it_M : vamos_fragment.cuts_M) {
        for (const auto &it_Z : vamos_fragment.cuts_Z) {
            for (const auto &condition : AGATAconditions) {
                Istantiate(pData.AGATA.hDC[it_M][it_Z][condition],
                           new TH1D(Form("pData_AGATA_hDC_M%i_Z%i_cond%s", it_M, it_Z, condition.c_str()),
                                    Form("DC spectrum of M%i Z%i with condition %s", it_M, it_Z, condition.c_str()),
                                    4000, 0, 4000));
            }
            Istantiate(pData.AGATA.mDC[it_M][it_Z],
                       new TH2D(Form("pData_AGATA_mDC_M%i_Z%i", it_M, it_Z), Form("DC gamma gamma of M%i Z%i", it_M, it_Z),
                                4000, 0, 4000, 4000, 0, 4000));

            Istantiate(pData.AGATA.mDC_ThetaMUGAST[it_M][it_Z],
                       new TH2D(Form("pData_AGATA_mDC_ThetaMUGAST_M%i_Z%i", it_M, it_Z),
                                Form("DC gamma vs Theta on MUGAST of M%i Z%i", it_M, it_Z),
                                4000, 0, 4000, 180, 0, 180));

            for (const auto &particle : particles) {
                Istantiate(pData.AGATA.mEx_DC[it_M][it_Z][particle],
                           new TH2D(Form("pData_AGATA_mEx_DC_M%i_Z%i_%s", it_M, it_Z, particle.c_str()),
                                    Form("Excitation energy AGATA vs MUGAST M%i Z%i and %s", it_M, it_Z, particle.c_str()),
                                    1000, -10, 10, 4000, 0, 4));

                Istantiate(pData.AGATA.mELab_ThetaLab[it_M][it_Z][particle],
                           new TH2D(Form("pData_AGATA_ELab_ThetaLab_M%i_Z%i_%s", it_M, it_Z, particle.c_str()),
                                    Form("Excitation energy AGATA vs MUGAST M%i Z%i and %s", it_M, it_Z, particle.c_str()),
                                    1000, 0, 10, 1000, 0, 10));
            }
        }
    }

#ifdef VERBOSE_DEBUG
    std::cout << "------------>finished: agata initialization\n";
#endif

    //CATS

    //MUGAST
    for (const auto &it_MG : mugast_fragment.cuts_MG) {
        //E TOF////////////////////////////////////////
        Istantiate(pConf.MG.mE_TOF[it_MG],
                   new TH2D(Form("pConf_MG_mE_TOF_MG%i", it_MG),
                            Form("E vs TOF of MG%i", it_MG),
                            1000, 0, 28, 1000, 260, 380));

        Istantiate(pConf.MG.mE_TOF2[it_MG],
                   new TH2D(Form("pConf_MG_mE_TOF2_MG%i", it_MG),
                            Form("E vs TOF2 of MG%i", it_MG),
                            1000, 0, 28, 1000, 260, 380));

        for (const auto &strip : mugast_fragment.strips) {
            //Strip E////////////////////////////////////////
            Istantiate(pConf.MG.mStrip_E[it_MG][strip],
                       new TH2D(Form("pConf_SI_mStrip_E_MG%i_%s", it_MG, strip.c_str()),
                                Form("Strip vs E of MG%i %s", it_MG, strip.c_str()),
                                129, 0, 121, 1000, 0, 30));

            //Strip T////////////////////////////////////////
            Istantiate(pConf.MG.mStrip_T[it_MG][strip],
                       new TH2D(Form("pConf_SI_mStrip_T_MG%i_%s", it_MG, strip.c_str()),
                                Form("Strip vs T of MG%i %s", it_MG, strip.c_str()),
                                129, 0, 129, 1000, 0, 1500));

            Istantiate(pConf.MG.mStrip_T2[it_MG][strip],
                       new TH2D(Form("pConf_SI_mStrip_T2_MG%i_%s", it_MG, strip.c_str()),
                                Form("Strip vs T2 of MG%i %s", it_MG, strip.c_str()),
                                129, 0, 129, 1000, 0, 1500));
        }

        for (const auto &it_M : vamos_fragment.cuts_M) {
            for (const auto &it_Z : vamos_fragment.cuts_Z) {
                Istantiate(pData.MG.mE_TOF[it_M][it_Z][it_MG],
                           new TH2D(Form("pData_MG_mE_TOF_M%i_Z%i_MG%i", it_M, it_Z, it_MG),
                                    Form("E vs TOF of MG%i with M%i Z%i", it_MG, it_M, it_Z),
                                    1000, 0, 28, 1000, 260, 380));

                Istantiate(pData.MG.mE_TOF2[it_M][it_Z][it_MG],
                           new TH2D(Form("pData_MG_mE_TOF2_M%i_Z%i_MG%i", it_M, it_Z, it_MG),
                                    Form("E vs TOF2 of MG%i with M%i Z%i", it_MG, it_M, it_Z),
                                    1000, 0, 28, 1000, 260, 380));
            }
        }
    }
    for (const auto &it_M : vamos_fragment.cuts_M) {
        for (const auto &it_Z : vamos_fragment.cuts_Z) {
            for (const auto &particle : mugast_fragment.particles) {
                Istantiate(pData.MG.hEx[it_M][it_Z][particle],
                           new TH1D(Form("pData_MG_hEx_M%i_Z%i_%s", it_M, it_Z, particle.c_str()),
                                    Form("Excitation energy with M%i Z%i in VAMOS and %s in MUGAST", it_M, it_Z, particle.c_str()),
                                    1000, -60, 60));
                for (const auto &it_gamma : gammas) {
                    Istantiate(pData.MG.mELab_ThetaLab[it_M][it_Z][particle][it_gamma],
                               new TH2D(Form("pData_MG_mELab_ThetaLab_M%i_Z%i_%s_%s", it_M, it_Z, particle.c_str(), it_gamma.c_str()),
                                        Form("ELab vs Theta Lab with M%i Z%i in VAMOS and %s in MUGAST and %s in AGATA", it_M, it_Z, particle.c_str(), it_gamma.c_str()),
                                        1000, 0, 3.1416, 1000, 0, 20));
                }
            }
            Istantiate(pData.MG.hEx[it_M][it_Z]["NONE"],
                       new TH1D(Form("pData_SI_hEx_M%i_Z%i_%s", it_M, it_Z, "NONE"),
                                Form("Excitation energy with M%i Z%i in VAMOS and %s in MUGAST", it_M, it_Z, "NONE"),
                                1000, -60, 60));
        }
    }
    //TODO: Must2
    //for (const auto &it_MM : must2_fragment.cuts_MM) {
    //    Istantiate(pConf.MM.mE_TOF[it_MM],
    //                new TH2D(Form("pConf-MM-mE_TOF-MM%i", it_MM),
    //                            Form("E vs TOF of MM%i", it_MM),
    //                            1000, 0, 28, 1000, 260, 380));
    //    for (const auto &it_M : vamos_fragment.cuts_M) {
    //        for (const auto &it_Z : vamos_fragment.cuts_Z) {
    //            Istantiate(pData.MM.mdE_E_Si[it_M][it_Z][it_MM],
    //                        new TH2D(Form("pData_SI_mdE_E_Si_M%i_Z%i_%s", it_M, it_Z, it_MM),
    //                                    Form("dE E of %s with M%i Z%i", MM.c_str(), it_M, it_Z),
    //                                    1000, 0, 28, 1000, 0, 28));
    //
    //            Istantiate(pData.MM.mE_TOF[it_M][it_Z][it_MM],
    //                        new TH2D(Form("pData_MM_mE_TOF_M%i_Z%i_MM%i", it_M, it_Z, it_MM),
    //                                    Form("E vs TOF of MM%i with M%i Z%i", it_MM, it_M, it_Z),
    //                                    1000, 0, 28, 1000, 260, 380));
    //        }
    //    }
    //}

    //TODO: Add MUST2

    //Istantiate(pData.SI.mE_TOF[it_M][it_Z]["MG"],
    //            new TH2D(Form("pData_SI_mE_TOF_MG_M%i_Z%i", it_M, it_Z),
    //                        Form("E vs TOF of all MG with M%i Z%i", it_M, it_Z),
    //                        1000, 0, 28, 1000, 260, 380));
    //            for (const auto &particle : particles) {
    //                Istantiate(pData.SI.hEx[it_M][it_Z][particle],
    //                            new TH1D(Form("pData_SI-hEx_M%i_Z%i_%s", it_M, it_Z, particle.c_str()),
    //                                        Form("Excitation energy with M%i Z%i in VAMOS and %s in MUGAST", it_M, it_Z, particle.c_str()),
    //                                        1000, -60, 60));
    //
    //                Istantiate(pData.SI.mEx_TW[it_M][it_Z][particle],
    //                            new TH2D(Form("pData_SI_mEx_TW_M%i_Z%i_%s", it_M, it_Z, particle.c_str()),
    //                                        Form("Excitation energy vs Time with M%i Z%i in VAMOS and %s in MUGAST", it_M, it_Z, particle.c_str()),
    //                                        5000, 242, 328, 1000, -60, 60));
    //
    //                Istantiate(pData.SI.mECM_ThetaCM[it_M][it_Z][particle],
    //                            new TH2D(Form("pData_SI_mECM_ThetaCM_M%i_Z%i_%s", it_M, it_Z, particle.c_str()),
    //                                        Form("E CM vs Theta CM with M%i Z%i in VAMOS and %s in MUGAST", it_M, it_Z, particle.c_str()),
    //                                        1000, 0, 180, 1000, 0, 60));
    //
    //                for (const auto &gamma : gammas) {
    //                    Istantiate(pData.SI.mELab_ThetaLab[it_M][it_Z][particle][gamma],
    //                                new TH2D(Form("pData_SI_mELab_ThetaLab_M%i_Z%i_%s_%s", it_M, it_Z, particle.c_str(), gamma.c_str()),
    //                                            Form("ELab vs Theta Lab with M%i Z%i in VAMOS and %s in MUGAST and %s in AGATA", it_M, it_Z, particle.c_str(), gamma.c_str()),
    //                                            1000, 0, 180, 1000, 0, 60));
    //                }
    //            }
    //        }
    //    }
    //    for (const auto &particle : particles) {
    //        Istantiate(pData.SI.mELab_ThetaLab[0][0][particle]["ANY"],
    //                    new TH2D(Form("pData-SI-mELab_ThetaLab-%s-%s-%s", "ANY", "ANY", particle.c_str()),
    //                                Form("E Lab vs Theta Lab with %s %s in VAMOS and %s in MUGAST", "ANY", "ANY", particle.c_str()),
    //                                1000, 0, 180, 1000, 0, 60));
    //     }

    if (new_graph_file != nullptr) {
        new_graph_file->close();
        //TODO: kill all threads here
    }

#ifdef VERBOSE_DEBUG
    std::cout << "------------>finished: SI initialization\n";
#endif

    ///Temporary histogram///////////////////////////////////////////////////////////////////////////
    general_histo_ptr = new TH2D("test", "test", 1000, -400, 400, 1000, 0, 100);
    fOutput->Add(general_histo_ptr);

#ifdef VERBOSE_DEBUG
    std::cout << "------------>finished: Selector::SlaveBegin()\n";
#endif
}

Bool_t Selector::Process(Long64_t entry) {
    fReader.SetLocalEntry(entry);
   if (entry%5000 == 0) std::cout << "Processed entries : " <<entry <<" of  "<<total_entries<<"\n"; 

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
#ifdef VERBOSE_DEBUG
    std::cout << "------------>Finished: Plotting VAMOS graphs\n";
#endif
    if (!vamos_fragment.Identified()) goto mugast_label;

#ifdef VERBOSE_DEBUG
    std::cout << "------------>Finished: Vamos identification, positive exit\n";
#endif

    //AGATA///////////////////////////////////////////////////////////////////////////////////////////////////
    if (AGATA_GOOD) {
        for (long unsigned int ii = 0; ii < AddE.GetSize(); ii++) {
            //if (AddE[ii] > 10 && *GATCONF_MASTER==1) {
            if (AddE[ii] > 10) {
                Fill(pData.AGATA.hDC[vamos_fragment.Get_id_M()][vamos_fragment.Get_id_Z()]["NONE"],
                     1E3 * CorrectDoppler(*vamos_fragment.Get_p4(),
                                          AddE[ii] / 1E3, AddX[ii], AddY[ii], AddZ[ii]));

                for (long unsigned int kk = 0; kk < (*Mugast).PosX.size(); kk++) {
                    TVector3 vec((*Mugast).PosX[kk], (*Mugast).PosY[kk], (*Mugast).PosZ[kk]);

                    Fill(pData.AGATA.mDC_ThetaMUGAST[vamos_fragment.Get_id_M()][vamos_fragment.Get_id_Z()],
                         1E3 * CorrectDoppler(*vamos_fragment.Get_p4(),
                                              AddE[ii] / 1E3, AddX[ii], AddY[ii], AddZ[ii]),
                         vec.Theta() * TMath::RadToDeg());
                }
                //Gamma Gamma matrices
                for (long unsigned int jj = 0; jj < AddE.GetSize(); jj++) {
                    if (AddE[jj] > 10 && ii != jj) {
                        Fill(pData.AGATA.mDC[vamos_fragment.Get_id_M()][vamos_fragment.Get_id_Z()],
                             1E3 * CorrectDoppler(*vamos_fragment.Get_p4(), AddE[ii] / 1E3,
                                                  AddX[ii], AddY[ii], AddZ[ii]),
                             1E3 * CorrectDoppler(*vamos_fragment.Get_p4(), AddE[jj] / 1E3,
                                                  AddX[jj], AddY[jj], AddZ[jj]));

                        Fill(pData.AGATA.mDC[vamos_fragment.Get_id_M()][vamos_fragment.Get_id_Z()],
                             1E3 * CorrectDoppler(*vamos_fragment.Get_p4(), AddE[jj] / 1E3,
                                                  AddX[jj], AddY[jj], AddZ[jj]),
                             1E3 * CorrectDoppler(*vamos_fragment.Get_p4(), AddE[ii] / 1E3,
                                                  AddX[ii], AddY[ii], AddZ[ii]));
                    }
                }
            }
        }
    }

//MUGAST//////////////////////////////////////////////////////////////////////////////////////////////////////
mugast_label:  //Label of goto previous to VAMOS
    LoadMugastData();
#ifdef VERBOSE_DEBUG
    std::cout << "------------>Finished: Loading Mugast data\n";
#endif
    mugast_fragment.Identify();
#ifdef VERBOSE_DEBUG
    std::cout << "------------>Finished: Mugast identification, negative exit\n";
#endif

    PlotMugastGraphs();

    //SI data loops
//    try {
//E-TOF plots
#ifdef VERBOSE_DEBUG
    std::cout << "------------>Finished: E TOF histos filled\n";
#endif

    //Loop on physics data
//        int MugastEvents = 0;
//        for (long unsigned int ii = 0; ii < DetID.GetSize(); ii++) {
//            Fill(pData.SI.mELab_ThetaLab[0][0]["ANY"]["ANY"], ThetaLab[MugastEvents], ELab[MugastEvents]);
//            if (DetID[ii] >= 100)  //These events are in Must2
//                continue;
//            if (DetID[ii] != (*Mugast).TelescopeNumber[MugastEvents])
//                std::cout << "Something is wrong matching Mugast events\n";
//Excitation energy and kinematic lines

//            Fill(pData.SI.hEx           [vamos_fragment.Get_id_M()][vamos_fragment.Get_id_Z()]["ANY"],
//                    Ex[MugastEvents]);
//
//            Fill(pData.SI.mEx_TW        [vamos_fragment.Get_id_M()][vamos_fragment.Get_id_Z()]["ANY"],
//                    *TW, Ex[MugastEvents]);
//
//            Fill(pData.SI.mELab_ThetaLab[vamos_fragment.Get_id_M()][vamos_fragment.Get_id_Z()]["ANY"]["ANY"],
//                    ThetaLab[MugastEvents], ELab[MugastEvents]);
//
//            Fill(pData.SI.mECM_ThetaCM  [vamos_fragment.Get_id_M()][vamos_fragment.Get_id_Z()]["ANY"],
//                    ThetaCM[MugastEvents], Ecm[MugastEvents]);

//            if (AGATA_GOOD) {
//                //Loop over gammas
//                for (long unsigned int ii = 0; ii < AddE.GetSize(); ii++) {
//                    if (AddE[ii] > 10) {
//                        Fill(pData.AGATA.mEx_DC     [vamos_fragment.Get_id_M()][vamos_fragment.Get_id_Z()]["ANY"],
//                                Ex[MugastEvents],
//                                CorrectDoppler(*vamos_fragment.Get_p4(), AddE[ii] / 1E3, AddX[ii], AddY[ii], AddZ[ii]));
//                    }
//                    if (AddE[ii] > 320 && AddE[ii] < 390) {
//                        Fill(pData.SI.mELab_ThetaLab[vamos_fragment.Get_id_M()][vamos_fragment.Get_id_Z()]["ANY"]["360 keV"],
//                                ThetaLab[MugastEvents], ELab[MugastEvents]);
//                    }
//                }
//            }
#ifdef VERBOSE_DEBUG
    std::cout << "------------>Finished: kinematic histos filled\n";
#endif
    //if (cut.at("MUGAST").at(Form("CUT_ETOF_MG%d", (*Mugast).TelescopeNumber[MugastEvents]))->IsInside((*Mugast).DSSD_E[MugastEvents], AlignT((*Mugast).TelescopeNumber[MugastEvents], (*Mugast).DSSD_Y[MugastEvents], (*Mugast).DSSD_T[MugastEvents])))
    //            for (const auto &particle : particles) {
    //                if (particle == "ANY")
    //                    continue;
    //                if (cut.at("MUGAST").at(Form("E_TOF_MG%d_%s", (*Mugast).TelescopeNumber[MugastEvents], particle.c_str()))->IsInside((*Mugast).DSSD_E[MugastEvents], (*Mugast).DSSD_T[MugastEvents])) {
    //                    if (AGATA_GOOD) {
    //                        //Loop over gammas
    //                        for (long unsigned int ii = 0; ii < AddE.GetSize(); ii++) {
    //                            if (AddE[ii] > 10) {
    //                                Fill(pData.AGATA.mEx_DC[vamos_fragment.Get_id_M()][vamos_fragment.Get_id_Z()][particle],
    //                                        Ex[MugastEvents], CorrectDoppler(*vamos_fragment.Get_p4(),
    //                                        AddE[ii] / 1E3, AddX[ii], AddY[ii], AddZ[ii]));
    //                            }
    //                            if (AddE[ii] > 320 && AddE[ii] < 390) {
    //                                if (particle == "2_H")
    //                                    Fill(pData.SI.mELab_ThetaLab[vamos_fragment.Get_id_M()][vamos_fragment.Get_id_Z()][particle]["360 keV"],
    //                                            ThetaLab[MugastEvents], ELab[MugastEvents]);
    //                            }
    //                        }
    //                    }
    //                    Fill(pData.SI.hEx           [vamos_fragment.Get_id_M()][vamos_fragment.Get_id_Z()][particle],
    //                            Ex[MugastEvents]);
    //
    //                    Fill(pData.SI.mEx_TW        [vamos_fragment.Get_id_M()][vamos_fragment.Get_id_Z()][particle],
    //                            *TW, Ex[MugastEvents]);
    //
    //                    Fill(pData.SI.mELab_ThetaLab[vamos_fragment.Get_id_M()][vamos_fragment.Get_id_Z()][particle]["ANY"],
    //                            ThetaLab[MugastEvents], ELab[MugastEvents]);
    //
    //                    Fill(pData.SI.mECM_ThetaCM  [vamos_fragment.Get_id_M()][vamos_fragment.Get_id_Z()][particle],
    //                            ThetaCM[MugastEvents], Ecm[MugastEvents]);
    //
    //                    Fill(pData.SI.mELab_ThetaLab[0][0][particle]["ANY"],
    //                            ThetaLab[MugastEvents], ELab[MugastEvents]);
    //                }
    //            }
    //            MugastEvents++;
    //        }
    //        //Cats///////////////////////////////////////////////////////////////////////////////////////////////////
    //        for (long unsigned int ii = 0; ii < (*CATS).PositionX.size(); ii++) {
    //            Fill(pConf.CATS.mCATSpos, (*CATS).PositionX[ii], (*CATS).PositionY[ii]);
    //        }
    //    } catch (std::out_of_range &e) {
    //        std::cerr << "Silicon Physics loop :" << e.what() << std::endl;
    //    }

    //MUST2
//    for (long unsigned int ii = 0; ii < (*MUST2).Si_E.size(); ii++) {
//        Fill(pConf.SI.mE_TOF[Form("MM%i", (*MUST2).TelescopeNumber[ii])],
//                (*MUST2).Si_E[ii], (*MUST2).Si_T[ii]);
//        if (vamos_fragment.Identified()) {
//            Fill(pData.SI.mE_TOF[vamos_fragment.Get_id_M()][vamos_fragment.Get_id_Z()][Form("MM%i",(*MUST2).TelescopeNumber[ii])],
//                 (*MUST2).Si_E[ii], (*MUST2).Si_T[ii]);
//        }
//    }
#ifdef VERBOSE_DEBUG
    std::cout << "------------>Finished: Si histos filled\n";
#endif

    //AGATA
    // Fill(pConf.AGATA.hAddTS_LTS,*AddTS-*LTS);
    for (int ii = 0; ii < *nbAdd; ii++) {
        Fill(pConf.AGATA.mmAGATA3D, AddX[ii], AddY[ii], AddZ[ii]);
    }
    return kTRUE;
}

void Selector::SlaveTerminate() {
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
#ifdef VERBOSE_DEBUG
    std::cout << "------------>Selector::Terminate()\n";
#endif
}

inline void Selector::LoadVamosData() {
    vamos_fragment.SetData(new VamosIdentification::Data(&IC, &Path, &Brho, &Xf,
                                                         &ThetaL, &PhiL, &AGAVA_VAMOSTS,
                                                         &T_FPMW_CATS2_C));
}

inline void Selector::LoadMugastData() {
    mugast_fragment.SetData(new MugastIdentification::Data(&Mugast, 
                                                            &TW,
                                                            vamos_fragment.Get_id_M(),
                                                            vamos_fragment.Get_id_Z()));
}

inline void Selector::PlotVamosGraphs() {
#ifdef VERBOSE_DEBUG
    std::cout << "------------>Selector::PlotVamosGraphs()\n";
#endif
    //dE-E plot, no conditions
    Fill(pConf.VAMOS.mdE_E,
         vamos_fragment.Get_En(), vamos_fragment.Get_D_En());
    Fill(pConf.VAMOS.mdE2_E,
         vamos_fragment.Get_En(), vamos_fragment.Get_D_En2());

    if (vamos_fragment.Get_id_Z() == 0) return;
    Fill(pConf.VAMOS.mQ_MQ[vamos_fragment.Get_id_Z()],
         vamos_fragment.Get_M_Q(), vamos_fragment.Get_Charge());
    Fill(pConf.VAMOS.Xf_MQ[vamos_fragment.Get_id_Z()],
         vamos_fragment.Get_M_Q(), *Xf);
    if (!vamos_fragment.Identified()) return;
    Fill(pData.VAMOS.mTW_Brho[vamos_fragment.Get_id_M()][vamos_fragment.Get_id_Z()],
         *TW, *Brho);
}

inline void Selector::PlotMugastGraphs() {
#ifdef VERBOSE_DEBUG
    std::cout << "------------>Selector::PlotMugastGraphs()\n";
#endif
    for (int ii = 0; ii < mugast_fragment.Get_Mult(); ++ii) {
        //Ex
        Fill(pData.MG.hEx[vamos_fragment.Get_id_M()]
                         [vamos_fragment.Get_id_Z()]
                         [mugast_fragment.Get_Particle(ii)],
             mugast_fragment.Get_Ex(ii));

        //ELab Theta Lab
        Fill(pData.MG.mELab_ThetaLab[vamos_fragment.Get_id_M()]
                                    [vamos_fragment.Get_id_Z()]
                                    [mugast_fragment.Get_Particle(ii)]
                                    ["NOCONDITION"],
             mugast_fragment.Get_ThetaLab(ii),
             mugast_fragment.Get_E(ii));

        //E TOF
        Fill(pConf.MG.mE_TOF[mugast_fragment.Get_MG(ii)],
             mugast_fragment.Get_SI_E(ii),
             mugast_fragment.Get_T(ii));

        Fill(pConf.MG.mE_TOF2[mugast_fragment.Get_MG(ii)],
             mugast_fragment.Get_SI_E(ii),
             mugast_fragment.Get_T2(ii));
        //Strip
        //E
        Fill(pConf.MG.mStrip_E[mugast_fragment.Get_MG(ii)]["X"],
             mugast_fragment.Get_SI_X(ii),
             mugast_fragment.Get_SI_E(ii));

        Fill(pConf.MG.mStrip_E[mugast_fragment.Get_MG(ii)]["Y"],
             mugast_fragment.Get_SI_Y(ii),
             mugast_fragment.Get_SI_E(ii));
        //T
        Fill(pConf.MG.mStrip_T[mugast_fragment.Get_MG(ii)]["X"],
             mugast_fragment.Get_SI_X(ii),
             mugast_fragment.Get_T(ii));

        Fill(pConf.MG.mStrip_T[mugast_fragment.Get_MG(ii)]["Y"],
             mugast_fragment.Get_SI_Y(ii),
             mugast_fragment.Get_T(ii));
        //T2
        Fill(pConf.MG.mStrip_T2[mugast_fragment.Get_MG(ii)]["X"],
             mugast_fragment.Get_SI_X(ii),
             mugast_fragment.Get_T2(ii));

        Fill(pConf.MG.mStrip_T2[mugast_fragment.Get_MG(ii)]["Y"],
             mugast_fragment.Get_SI_Y(ii),
             mugast_fragment.Get_T2(ii));

        if (vamos_fragment.Identified()) {
            Fill(pData.MG.mE_TOF[vamos_fragment.Get_id_M()][vamos_fragment.Get_id_Z()][mugast_fragment.Get_MG(ii)],
                 mugast_fragment.Get_SI_E(ii),
                 mugast_fragment.Get_T(ii));

            Fill(pData.MG.mE_TOF2[vamos_fragment.Get_id_M()][vamos_fragment.Get_id_Z()][mugast_fragment.Get_MG(ii)],
                 mugast_fragment.Get_SI_E(ii),
                 mugast_fragment.Get_T2(ii));
        }
    }

    //        if (vamos_fragment.Identified()) {
    //            for (long unsigned int ii = 0; ii < (*Mugast).DSSD_E.size(); ii++) {
    //TVector3 hitPos((*Mugast).PosX[ii], (*Mugast).PosY[ii], (*Mugast).PosZ[ii]);
    //Fill(pData.SI.mE_TOF[vamos_fragment.Get_id_M()][vamos_fragment.Get_id_Z()][Form("MG%i", (*Mugast).TelescopeNumber[ii])], (*Mugast).DSSD_E[ii], AlignT((*Mugast).TelescopeNumber[ii], (*Mugast).DSSD_Y[ii], (*Mugast).DSSD_T[ii]));

    //                Fill(pData.SI.mE_TOF[vamos_fragment.Get_id_M()][vamos_fragment.Get_id_Z()][Form("MG%i", (*Mugast).TelescopeNumber[ii])],
    //                        mugast_fragment.Get_SI_E(ii), mugast_fragment.Get_T(ii));

    //Fill(pData.SI.mE_TOF[vamos_fragment.Get_id_M()][vamos_fragment.Get_id_Z()]["MG"],
    //        (*Mugast).DSSD_E[ii], AlignPunch((*Mugast).TelescopeNumber[ii],
    //        (*Mugast).DSSD_T[ii]));
    //            }
    //        }
}

inline double Selector::CorrectDoppler(const TLorentzVector &p4, const double &Egamma,
                                       const double &X, const double &Y, const double &Z) {
    TLorentzVector pgamma(Egamma, 0, 0, Egamma);
    TVector3 PosGamma(X, Y, Z + agata_Zshift);
    pgamma.SetPhi(PosGamma.Phi());
    pgamma.SetTheta(PosGamma.Theta());
    pgamma.SetE(Egamma);
    pgamma.Boost(-p4.BoostVector());
    return pgamma.Energy();
}

bool Selector::GetSettings() {
    std::string config_file = "./Configs/GraphsEnabled.txt";
    std::ifstream ifile(config_file);
    if (ifile) {
        //       std::ifstream file(config_file);
        std::string line;
        while (std::getline(ifile, line)) {
            std::istringstream str(line);
            std::string Graph;
            str >> Graph;
            std::string enabled;
            str >> enabled;
            if (enabled.compare("false")) {  //IS ENABLED
                enabled_histograms[Graph] = true;
            }
        }
    } else {
        ifile.close();
        new_graph_file = new std::ofstream(config_file);
    }
    return true;
}

inline bool Selector::Fill(TH1 *histo, const double &data1) {
    if (histo == nullptr) {
        return false;
    } else {
#ifndef VERBOSE_DEBUG
        histo->Fill(data1);
#endif
#ifdef VERBOSE_DEBUG
        int i = 0;
        i = histo->Fill(data1);
        if (1) {
            std::cout << "Name : " << histo->GetName() << std::endl;
            std::cout << "Values : " << data1 << std::endl;
        }
#endif
    }
    return true;
}

inline bool Selector::Fill(TH2 *histo,
                           const double &data1, const double &data2) {
    if (histo == nullptr) {
        return false;
    } else {
#ifndef VERBOSE_DEBUG
        histo->Fill(data1, data2);
#endif
#ifdef VERBOSE_DEBUG
        int i = 0;
        i = histo->Fill(data1, data2);
        if (1) {
            std::cout << "Name : " << histo->GetName() << std::endl;
            std::cout << "Values : " << data1 << "  " << data2 << std::endl;
        }
#endif
    }
    return true;
}

inline bool Selector::Fill(TH3 *histo, const double &data1,
                           const double &data2, const double &data3) {
    if (histo == nullptr) {
        return false;
    } else {
#ifndef VERBOSE_DEBUG
        histo->Fill(data1, data2, data3);
#endif
#ifdef VERBOSE_DEBUG
        int i = 0;
        i = histo->Fill(data1, data2, data3);
        if (1) {
            std::cout << "Name : " << histo->GetName() << std::endl;
            std::cout << "Values : " << data1 << "  " << data2 << "   " << data3 << std::endl;
        }
#endif
    }
    return true;
}
