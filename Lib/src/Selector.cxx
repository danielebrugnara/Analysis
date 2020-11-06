#include "Selector.h"

#ifdef VERBOSE_DEBUG
#  define DEBUG(x, y) std::cout << (x) << (y) << std::endl;
#else
#  define DEBUG(x, y) 
#endif

//Constructor///////////////////////////////////////////////////////////////////////////
Selector::Selector(TTree *) {}

//Distructor////////////////////////////////////////////////////////////////////////////
Selector::~Selector(){
//    fOutput->Clear();
//    TIter iter(fOutput);
//    TObject *obj;
//    while (iter != iter.End()){
//        obj = iter();
//        delete obj; //No need as we are using smart pointers
//    }
}

//Begn//////////////////////////////////////////////////////////////////////////////////
void Selector::Begin(TTree * /*tree*/){
    DEBUG("------------>Selector::Begin()", "");
    TString option = GetOption();
}

//Slave Begin///////////////////////////////////////////////////////////////////////////
void Selector::SlaveBegin(TTree * /*tree*/){
    DEBUG("------------>Selector::SlaveBegin()", "");
    total_entries = fReader.GetEntries(false);

    TString option = GetOption();
    file_name = "./Out/" + option;

    ///Deciding which graphs to plot/////////////////////////////////////////////////////////////////
    GetSettings(); //Decides which histograms to fill
    ///Loading Identification classes/////////////////////////////////////////////////////////////

    //TODO: use constructor to initialize this!!!
    std::string VAMOS_cuts_file = "./Configs/Cuts/VAMOS.root";
    vamos_fragment.LoadCuts(VAMOS_cuts_file);
    vamos_fragment.Initialize();

    std::string MUGAST_cuts_file = "./Configs/Cuts/MUGAST.root";
    mugast_fragment.LoadCuts(MUGAST_cuts_file);

    //Passing beam energy in MeV, target position mm
    //mugast_fragment.Initialize(379.04, TVector3(0, 0, 25.));
    mugast_fragment.Initialize(450.00*UNITS::MeV, TVector3(0, 0, 25.));

    DEBUG("------------>finished: vamos_fragment initialization", "");

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

    for (const auto &Z_it : vamos_fragment.cuts_Z){
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

    //Data histograms//////////////////////////////////////////////////////////////////////////////
    for (const auto &it_M : vamos_fragment.cuts_M){
        for (const auto &it_Z : vamos_fragment.cuts_Z){
            Istantiate(pData.VAMOS.mTW_Brho[it_M][it_Z],
                       new TH2D(Form("pData_VAMOS_mTW_Brho_%i_%i", it_M, it_Z),
                                Form("Time vs Brho with M%i Z%i in VAMOS", it_M, it_Z),
                                //5000, 242, 328, 1000, 0.5, 1.5));
                                5000, 0, 450, 1000, 0.5, 1.5));

            for (const auto &particle : mugast_fragment.light_particles){
                Istantiate(pData.VAMOS.mE_Theta[it_M][it_Z][particle],
                       new TH2D(Form("pData_VAMOS_mE_Theta__%i_%i_%s", it_M, it_Z, particle.c_str()),
                                Form("Energy (Brho) vs Theta with M%i Z%i in VAMOS and %s in MUGAST", it_M, it_Z, particle.c_str()),
                                1000, 0, 1, 1000, 1, 450));
            }
        }
    }
    //VAMOS

    //AGATA
    for (const auto &it_M : vamos_fragment.cuts_M){
        for (const auto &it_Z : vamos_fragment.cuts_Z){
            for (const auto &condition : AGATAconditions){
                Istantiate(pData.AGATA.hDC[it_M][it_Z][condition],
                           new TH1D(Form("pData_AGATA_hDC_M%i_Z%i_cond%s", it_M, it_Z, condition.c_str()),
                                    Form("DC spectrum of M%i Z%i with condition %s", it_M, it_Z, condition.c_str()),
                                    4000, 0, 4000));
            }
            Istantiate(pData.AGATA.mDC[it_M][it_Z],
                       new TH2D(Form("pData_AGATA_mDC_M%i_Z%i", it_M, it_Z), Form("DC gamma gamma of M%i Z%i", it_M, it_Z),
                                4000, 0, 4000, 4000, 0, 4000));
        }
    }

    DEBUG("------------>finished: agata initialization\n", "");

    //CATS

    //MUGAST

    Istantiate(pConf.MG.hit , 
                new TH3D("pConf_MG_hit",
                            "Hits on Mugast",
                            100, -150, 150,
                            100, -150, 150, 
                            100, -150, 150));

    Istantiate(pConf.MG.hit_XY , 
                new TH2D("pConf_MG_hit_XY",
                            "Hits on Mugast XY",
                            1000, -150, 150,
                            1000, -150, 150));

    Istantiate(pConf.MG.hit_XZ , 
                new TH2D("pConf_MG_hit_XZ",
                            "Hits on Mugast XZ",
                            1000, -150, 150,
                            1000, -150, 150));

    Istantiate(pConf.MG.hit_YZ , 
                new TH2D("pConf_MG_hit_YZ",
                            "Hits on Mugast YZ",
                            1000, -150, 150,
                            1000, -150, 150));

    for (const auto &it_MG : mugast_fragment.cuts_MG){
        //E TOF////////////////////////////////////////
        Istantiate(pConf.MG.mE_TOF[it_MG],
                   new TH2D(Form("pConf_MG_mE_TOF_MG%i", it_MG),
                            Form("E vs TOF of MG%i", it_MG),
                            1000, 0, 28, 1000, 260, 380));

        Istantiate(pConf.MG.mE_TOF2[it_MG],
                   new TH2D(Form("pConf_MG_mE_TOF2_MG%i", it_MG),
                            Form("E vs TOF2 of MG%i", it_MG),
                            1000, 0, 28, 1000, 260, 380));

        for (const auto &strip : mugast_fragment.strips){
            //Strip E////////////////////////////////////////
            Istantiate(pConf.MG.mStrip_E[it_MG][strip],
                       new TH2D(Form("pConf_SI_mStrip_E_MG%i_%s", it_MG, strip.c_str()),
                                Form("Strip vs E of MG%i %s", it_MG, strip.c_str()),
                                129, 0, 129, 1000, 0, 30));

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

        for (const auto &it_M : vamos_fragment.cuts_M){
            for (const auto &it_Z : vamos_fragment.cuts_Z){
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

    for (const auto &it_M : vamos_fragment.cuts_M){
        for (const auto &it_Z : vamos_fragment.cuts_Z){
            for (const auto &particle : mugast_fragment.light_particles){
                Istantiate(pData.MG.hEx[it_M][it_Z][particle],
                           new TH1D(Form("pData_MG_hEx_M%i_Z%i_%s", it_M, it_Z, particle.c_str()),
                                    Form("Excitation energy with M%i Z%i in VAMOS and %s in MUGAST", it_M, it_Z, particle.c_str()),
                                    1000, -60, 60));

                Istantiate(pData.MG.mEx_EDC[it_M][it_Z][particle],
                           new TH2D(Form("pData_MG_mEx_EDC_M%i_Z%i_%s", it_M, it_Z, particle.c_str()),
                                    Form("Ex vs EDC with M%i Z%i in VAMOS and %s in MUGAST", it_M, it_Z, particle.c_str()),
                                    1000, -20, 20, 2500, 0, 2500));

                for (const auto &it_gamma : gammas){
                    Istantiate(pData.MG.mELab_ThetaLab[it_M][it_Z][particle][it_gamma],
                               new TH2D(Form("pData_MG_mELab_ThetaLab_M%i_Z%i_%s_%s", it_M, it_Z, particle.c_str(), it_gamma.c_str()),
                                        Form("ELab vs Theta Lab with M%i Z%i in VAMOS and %s in MUGAST and %s in AGATA", it_M, it_Z, particle.c_str(), it_gamma.c_str()),
                                        1000, 0, 3.1416, 1000, 0, 20));

                    Istantiate(pData.MG.mEx_ThetaLab[it_M][it_Z][particle][it_gamma],
                               new TH2D(Form("pData_MG_mEx_ThetaLab_M%i_Z%i_%s_%s", it_M, it_Z, particle.c_str(), it_gamma.c_str()),
                                        Form("Ex vs Theta Lab with M%i Z%i in VAMOS and %s in MUGAST and %s in AGATA", it_M, it_Z, particle.c_str(), it_gamma.c_str()),
                                        1000, 0, 3.1416, 1000, -10, 10));

                    Istantiate(pData.MG.hThetaCM[it_M][it_Z][particle][it_gamma],
                               new TH1D(Form("pData_MG_hThetaCM_M%i_Z%i_%s_%s", it_M, it_Z, particle.c_str(), it_gamma.c_str()),
                                        Form("Theta CM with M%i Z%i in VAMOS and %s in MUGAST and %s in AGATA", it_M, it_Z, particle.c_str(), it_gamma.c_str()),
                                        1000, 0, 3.1415));
                }
            }
            Istantiate(pData.MG.hEx[it_M][it_Z]["NONE"],
                       new TH1D(Form("pData_SI_hEx_M%i_Z%i_%s", it_M, it_Z, "NONE"),
                                Form("Excitation energy with M%i Z%i in VAMOS and %s in MUGAST", it_M, it_Z, "NONE"),
                                1000, -60, 60));
        }
    }

    //TODO: Add MUST2

    if (new_graph_file != nullptr){
        new_graph_file->close();
        exit(1);
    }

    DEBUG("------------>finished: SI initialization", "");

    ///Temporary histogram///////////////////////////////////////////////////////////////////////////
    general_histo_ptr.reset(new TH2D("test", "Ex vs theta Lab", 1000, TMath::Pi()/2., TMath::Pi(), 1000, -10, 10));
    //fOutput->Add(general_histo_ptr.get());
    Output.push_back(general_histo_ptr.get());

    DEBUG("------------>finished: Selector::SlaveBegin()", "");
}

//Processing entry//////////////////////////////////////////////////////////////////////
Bool_t Selector::Process(Long64_t entry){
    DEBUG("------------>Selector::Process()", "");
    fReader.SetLocalEntry(entry);
    if (entry % 5000 == 0)
        std::cout << "\rProcessed entries : " << entry ;

    //Vamos/////////////////////////////////////////////
    LoadVamosData();
    DEBUG("------------>Finished: Loading Vamos data", "");
    vamos_fragment.Identify();
    DEBUG("------------>Finished: Identification", "");

    //Agata/////////////////////////////////////////////
    LoadAgataData();
    DEBUG("------------>Finished: Loading Agata data, positive exit", "");
    agata_gammas.Process();
    DEBUG("------------>Finished: Processing Agata data, positive exit", "");

    //MUGAST////////////////////////////////////////////
    LoadMugastData();
    DEBUG("------------>Finished: Loading Mugast data", "");
    mugast_fragment.Identify();
    DEBUG("------------>Finished: Mugast identification, negative exit", "");

    //Plotting histograms///////////////////////////////
    PlotVamosGraphs();
    DEBUG("------------>Finished: Plotting VAMOS graphs", "");
    PlotAgataGraphs();
    DEBUG("------------>Finished: Plotting agata graphs, positive exit", "");
    PlotMugastGraphs();
    DEBUG("------------>Finished: Analysis of one event", "");

    if (entry % 10000 == 0)
        mugast_fragment.StoreTWvsIce();

    return kTRUE;
}

void Selector::SlaveTerminate(){
    DEBUG("------------>Selector::SlaveTerminate()", "");
    auto *top = new TFile(file_name.c_str(), "recreate");
    std::cout << "Output file : " << file_name << "\n";
    for(const auto& it: Output){
        it->Write();
    }
    top->Close();
}

void Selector::Terminate(){
    DEBUG("------------>Selector::Terminate()", "");
}

inline void Selector::LoadVamosData(){
    vamos_fragment.SetData(new VamosIdentification::Data(&IC, &Path, &Brho, &Xf,
                                                         &ThetaL, &PhiL, &AGAVA_VAMOSTS,
                                                         &T_FPMW_CATS2_C));
}

inline void Selector::LoadMugastData(){
    mugast_fragment.SetData(new MugastIdentification::Data(&Mugast,
                                                           &CATS,
                                                           &TW,
                                                           vamos_fragment.Get_id_M(),
                                                           vamos_fragment.Get_id_Z()));
}

inline void Selector::LoadAgataData(){
    agata_gammas.SetData(new AgataProcessing::Data(&nbAdd,
                                                   &TSHit,
                                                   &AddTS,
                                                   &LTS,
                                                   &AddE,
                                                   &AddX,
                                                   &AddY,
                                                   &AddZ,
                                                   vamos_fragment.Get_p4()));
}

inline void Selector::PlotVamosGraphs(){
    DEBUG("------------>Selector::PlotVamosGraphs()", "");
    //dE-E plot, no conditions
    Fill(pConf.VAMOS.mdE_E.get(),vamos_fragment.Get_En(), vamos_fragment.Get_D_En());
    Fill(pConf.VAMOS.mdE2_E.get(),vamos_fragment.Get_En(), vamos_fragment.Get_D_En2());

    if (vamos_fragment.Get_id_Z() == 0)
        return;
    Fill(pConf.VAMOS.mQ_MQ[vamos_fragment.Get_id_Z()].get(),
         vamos_fragment.Get_M_Q(),
         vamos_fragment.Get_Charge());
    Fill(pConf.VAMOS.Xf_MQ[vamos_fragment.Get_id_Z()].get(),
         vamos_fragment.Get_M_Q(),
         *Xf);
    if (!vamos_fragment.Identified())
        return;
    Fill(pData.VAMOS.mTW_Brho[vamos_fragment.Get_id_M()][vamos_fragment.Get_id_Z()].get(),
         *TW, *Brho);

    for (unsigned int ii = 0; ii < mugast_fragment.Get_Mult(); ++ii){
        Fill(pData.VAMOS.mE_Theta[vamos_fragment.Get_id_M()][vamos_fragment.Get_id_Z()][mugast_fragment.Get_Particle(ii)].get(),
             vamos_fragment.Get_p4()->Angle(TVector3(0, 0, 1)),
            vamos_fragment.Get_EnFromBrho());
    }
}

inline void Selector::PlotMugastGraphs(){
    DEBUG("------------>Selector::PlotMugastGraphs()", "");
    for (unsigned int ii = 0; ii < mugast_fragment.Get_Mult(); ++ii){
        if (mugast_fragment.Get_T(ii)>260 && mugast_fragment.Get_T(ii)<380){
            //Only events with tof
            Fill(pConf.MG.hit.get(),
                    mugast_fragment.Get_Pos(ii)->X(),
                    mugast_fragment.Get_Pos(ii)->Y(),
                    mugast_fragment.Get_Pos(ii)->Z());

            Fill(pConf.MG.hit_XY.get(),
                    mugast_fragment.Get_Pos(ii)->X(),
                    mugast_fragment.Get_Pos(ii)->Y());

            Fill(pConf.MG.hit_XZ.get(),
                    mugast_fragment.Get_Pos(ii)->X(),
                    mugast_fragment.Get_Pos(ii)->Z());

            Fill(pConf.MG.hit_YZ.get(),
                    mugast_fragment.Get_Pos(ii)->Y(),
                    mugast_fragment.Get_Pos(ii)->Z());

        }
        //Ex
        Fill(pData.MG.hEx[vamos_fragment.Get_id_M()]
                         [vamos_fragment.Get_id_Z()]
                         [mugast_fragment.Get_Particle(ii)].get(),
             mugast_fragment.Get_Ex(ii));

        Fill(pData.MG.mELab_ThetaLab[vamos_fragment.Get_id_M()]
                                    [vamos_fragment.Get_id_Z()]
                                    [mugast_fragment.Get_Particle(ii)]
                                    ["NOCONDITION"].get(),
             mugast_fragment.Get_ThetaLab(ii),
             mugast_fragment.Get_E(ii));

        Fill(pData.MG.mEx_ThetaLab[vamos_fragment.Get_id_M()]
                                    [vamos_fragment.Get_id_Z()]
                                    [mugast_fragment.Get_Particle(ii)]
                                    ["NOCONDITION"].get(),
             mugast_fragment.Get_ThetaLab(ii),
             mugast_fragment.Get_Ex(ii));

        Fill(pData.MG.hThetaCM[vamos_fragment.Get_id_M()]
                                    [vamos_fragment.Get_id_Z()]
                                    [mugast_fragment.Get_Particle(ii)]
                                    ["NOCONDITION"].get(),
             mugast_fragment.Get_ThetaCM(ii));

        if (vamos_fragment.Get_id_M()==47 && vamos_fragment.Get_id_Z()==19 && mugast_fragment.Get_Particle(ii)=="m2_z1"){
        general_histo_ptr->Fill(
             mugast_fragment.Get_ThetaLab(ii),
             mugast_fragment.Get_Ex(ii));
        }
        if (agata_gammas.In_Coincidence()){
            //ELab Theta Lab
            for (long unsigned int jj = 0; jj < agata_gammas.Get_Mult(); ++jj){
                Fill(pData.MG.mEx_EDC[vamos_fragment.Get_id_M()]
                                     [vamos_fragment.Get_id_Z()]
                                     [mugast_fragment.Get_Particle(ii)].get(),
                     mugast_fragment.Get_Ex(ii),
                     agata_gammas.Get_EDC(jj));

                for (const auto &it_gamma : gammas){
                    if (it_gamma == "NOCONDITION")
                        continue;
                    if (abs(std::stod(it_gamma) - agata_gammas.Get_EDC(jj)) < gamma_gate){
                        //ELab ThetaLab with gamma condition
                        Fill(pData.MG.mELab_ThetaLab[vamos_fragment.Get_id_M()]
                                                    [vamos_fragment.Get_id_Z()]
                                                    [mugast_fragment.Get_Particle(ii)]
                                                    [it_gamma].get(),
                             mugast_fragment.Get_ThetaLab(ii),
                             mugast_fragment.Get_E(ii));

                        Fill(pData.MG.mEx_ThetaLab[vamos_fragment.Get_id_M()]
                                                    [vamos_fragment.Get_id_Z()]
                                                    [mugast_fragment.Get_Particle(ii)]
                                                    [it_gamma].get(),
                             mugast_fragment.Get_ThetaLab(ii),
                             mugast_fragment.Get_Ex(ii));

                        Fill(pData.MG.hThetaCM[vamos_fragment.Get_id_M()]
                                                    [vamos_fragment.Get_id_Z()]
                                                    [mugast_fragment.Get_Particle(ii)]
                                                    [it_gamma].get(),
                             mugast_fragment.Get_ThetaCM(ii));
                    }
                }
            }
        }

        //E TOF
        Fill(pConf.MG.mE_TOF[mugast_fragment.Get_MG(ii)].get(),
             mugast_fragment.Get_SI_E(ii),
             mugast_fragment.Get_T(ii));

        Fill(pConf.MG.mE_TOF2[mugast_fragment.Get_MG(ii)].get(),
             mugast_fragment.Get_SI_E(ii),
             mugast_fragment.Get_T2(ii));

            Fill(pConf.MG.mStrip_E[mugast_fragment.Get_MG(ii)]["X"].get(),
                 mugast_fragment.Get_SI_X(ii),
                 mugast_fragment.Get_SI_E(ii));

            Fill(pConf.MG.mStrip_T[mugast_fragment.Get_MG(ii)]["X"].get(),
                 mugast_fragment.Get_SI_X(ii),
                 mugast_fragment.Get_T(ii));

            Fill(pConf.MG.mStrip_T2[mugast_fragment.Get_MG(ii)]["X"].get(),
                 mugast_fragment.Get_SI_X(ii),
                 mugast_fragment.Get_T2(ii));

            Fill(pConf.MG.mStrip_E[mugast_fragment.Get_MG(ii)]["Y"].get(),
                 mugast_fragment.Get_SI_Y(ii),
                 mugast_fragment.Get_SI_E(ii));

            Fill(pConf.MG.mStrip_T[mugast_fragment.Get_MG(ii)]["Y"].get(),
                 mugast_fragment.Get_SI_Y(ii),
                 mugast_fragment.Get_T(ii));

            Fill(pConf.MG.mStrip_T2[mugast_fragment.Get_MG(ii)]["Y"].get(),
                 mugast_fragment.Get_SI_Y(ii),
                 mugast_fragment.Get_T2(ii));

        if (vamos_fragment.Identified())
        {
            Fill(pData.MG.mE_TOF[vamos_fragment.Get_id_M()][vamos_fragment.Get_id_Z()][mugast_fragment.Get_MG(ii)].get(),
                 mugast_fragment.Get_SI_E(ii),
                 mugast_fragment.Get_T(ii));

            Fill(pData.MG.mE_TOF2[vamos_fragment.Get_id_M()][vamos_fragment.Get_id_Z()][mugast_fragment.Get_MG(ii)].get(),
                 mugast_fragment.Get_SI_E(ii),
                 mugast_fragment.Get_T2(ii));
        }
    }
}

inline void Selector::PlotAgataGraphs(){
    for (int ii = 0; ii < *nbAdd; ii++){
        Fill(pConf.AGATA.mmAGATA3D.get(), AddX[ii], AddY[ii], AddZ[ii]);
    }
    if (agata_gammas.In_Coincidence()){
        for (unsigned int ii = 0; ii < agata_gammas.Get_Mult(); ++ii)
        {
            //Only agata
            Fill(pData.AGATA.hDC[vamos_fragment.Get_id_M()][vamos_fragment.Get_id_Z()]["NONE"].get(),
                 agata_gammas.Get_EDC(ii));
            for (unsigned int jj = 0; jj < agata_gammas.Get_Mult(); ++jj){
                if (ii == jj)
                    continue;

                Fill(pData.AGATA.mDC[vamos_fragment.Get_id_M()][vamos_fragment.Get_id_Z()].get(),
                     agata_gammas.Get_EDC(ii),
                     agata_gammas.Get_EDC(jj));
            }
        }
    }
}

inline void Selector::PlotCatsGraphs() {
    Fill(pConf.CATS.mCATSpos.get(),
         CATS->PositionOnTargetX,
         CATS->PositionOnTargetY);
}

bool Selector::GetSettings(){
    std::string config_file = "./Configs/GraphsEnabled.txt";
    std::ifstream ifile(config_file);
    if (ifile){
        std::string line;
        while (std::getline(ifile, line)){
            std::istringstream str(line);
            std::string Graph;
            str >> Graph;
            std::string enabled;
            str >> enabled;
            if (enabled !="false"){ //IS ENABLED
                enabled_histograms.emplace(Graph, true);
                //enabled_histograms[Graph] = true;

                std::string tmp_str;
                std::vector<std::string> settings;
                while(str>>tmp_str){
                    settings.push_back(tmp_str);
                }
                //if (settings.size()==0) continue;
                //if (settings.size()%3){
                //    enabled_histograms[Graph].lims = settings;
                //}else{
                //    std::cout   <<"Settings for : " 
                //                << Graph 
                //                <<" incorrect, keeping default ones\n";
                //}

            }
        }
    }else{
        ifile.close();
        new_graph_file = new std::ofstream(config_file);
    }
    return true;
}

std::vector<std::pair<double, double>> Selector::GetTWvsIce(){
    return mugast_fragment.GetTWvsIce();
}

inline bool Selector::Fill(TH1D* histo, const double &data1){
    if (histo == nullptr){
        return false;
    }else{
#ifndef VERBOSE_DEBUG
        histo->Fill(data1);
#endif
#ifdef VERBOSE_DEBUG
        int i = 0;
        i = histo->Fill(data1);
        if (1)
        {
            std::cout << "Name : " << histo->GetName() << std::endl;
            std::cout << "Values : " << data1 << std::endl;
        }
#endif
    }
    return true;
}

inline bool Selector::Fill(TH2D* histo,
                           const double &data1, const double &data2){
    if (histo == nullptr){
        return false;
    }else{
#ifndef VERBOSE_DEBUG
        histo->Fill(data1, data2);
#endif
#ifdef VERBOSE_DEBUG
        int i = 0;
        i = histo->Fill(data1, data2);
        if (1)
        {
            std::cout << "Name : " << histo->GetName() << std::endl;
            std::cout << "Values : " << data1 << "  " << data2 << std::endl;
        }
#endif
    }
    return true;
}

inline bool Selector::Fill(TH3D* histo, const double &data1,
                           const double &data2, const double &data3){
    if (histo == nullptr){
        return false;
    }else{
#ifndef VERBOSE_DEBUG
        histo->Fill(data1, data2, data3);
#endif
#ifdef VERBOSE_DEBUG
        int i = 0;
        i = histo->Fill(data1, data2, data3);
        if (1)
        {
            std::cout << "Name : " << histo->GetName() << std::endl;
            std::cout << "Values : " << data1 << "  " << data2 << "   " << data3 << std::endl;
        }
#endif
    }
    return true;
}
void Selector::Init(TTree *tree){
    fReader.SetTree(tree);
}

Bool_t Selector::Notify(){
    return kTRUE;
}