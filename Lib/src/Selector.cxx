#define Selector_cxx

#include "Selector.h"

//Constructor///////////////////////////////////////////////////////////////////////////
Selector::Selector(TTree *) {}

//Distructor////////////////////////////////////////////////////////////////////////////
Selector::~Selector()
{
    TIter iter(fOutput);
    TObject *obj;
    while ((obj = iter()))
    {
        delete obj;
    }
}

//Begn//////////////////////////////////////////////////////////////////////////////////
void Selector::Begin(TTree * /*tree*/)
{
#ifdef VERBOSE_DEBUG
    std::cout << "------------>Selector::Begin()\n";
#endif

    TString option = GetOption();
}

//Slave Begin///////////////////////////////////////////////////////////////////////////
void Selector::SlaveBegin(TTree * /*tree*/)
{
#ifdef VERBOSE_DEBUG
    std::cout << "------------>Selector::SlaveBegin()\n";
#endif
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
    mugast_fragment.Initialize(450.00, TVector3(0, 0, 25.));

#ifdef VERBOSE_DEBUG
    std::cout << "------------>finished: vamos_fragment initialization\n";
#endif

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

    for (const auto &Z_it : vamos_fragment.cuts_Z)
    {
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
    for (const auto &it_M : vamos_fragment.cuts_M)
    {
        for (const auto &it_Z : vamos_fragment.cuts_Z)
        {
            Istantiate(pData.VAMOS.mTW_Brho[it_M][it_Z],
                       new TH2D(Form("pData_VAMOS_mTW_Brho_%i_%i", it_M, it_Z),
                                Form("Time vs Brho with M%i Z%i in VAMOS", it_M, it_Z),
                                5000, 242, 328, 1000, 0.5, 1.5));
        }
    }
    //VAMOS

    //AGATA
    for (const auto &it_M : vamos_fragment.cuts_M)
    {
        for (const auto &it_Z : vamos_fragment.cuts_Z)
        {
            for (const auto &condition : AGATAconditions)
            {
                Istantiate(pData.AGATA.hDC[it_M][it_Z][condition],
                           new TH1D(Form("pData_AGATA_hDC_M%i_Z%i_cond%s", it_M, it_Z, condition.c_str()),
                                    Form("DC spectrum of M%i Z%i with condition %s", it_M, it_Z, condition.c_str()),
                                    4000, 0, 4000));
            }
            Istantiate(pData.AGATA.mDC[it_M][it_Z],
                       new TH2D(Form("pData_AGATA_mDC_M%i_Z%i", it_M, it_Z), Form("DC gamma gamma of M%i Z%i", it_M, it_Z),
                                4000, 0, 4000, 4000, 0, 4000));

            for (const auto &particle : particles)
            {
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
    for (const auto &it_MG : mugast_fragment.cuts_MG)
    {
        //E TOF////////////////////////////////////////
        Istantiate(pConf.MG.mE_TOF[it_MG],
                   new TH2D(Form("pConf_MG_mE_TOF_MG%i", it_MG),
                            Form("E vs TOF of MG%i", it_MG),
                            1000, 0, 28, 1000, 260, 380));

        Istantiate(pConf.MG.mE_TOF2[it_MG],
                   new TH2D(Form("pConf_MG_mE_TOF2_MG%i", it_MG),
                            Form("E vs TOF2 of MG%i", it_MG),
                            1000, 0, 28, 1000, 260, 380));

        for (const auto &strip : mugast_fragment.strips)
        {
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

        for (const auto &it_M : vamos_fragment.cuts_M)
        {
            for (const auto &it_Z : vamos_fragment.cuts_Z)
            {
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
    for (const auto &it_M : vamos_fragment.cuts_M)
    {
        for (const auto &it_Z : vamos_fragment.cuts_Z)
        {
            for (const auto &particle : mugast_fragment.light_particles)
            {
                Istantiate(pData.MG.hEx[it_M][it_Z][particle],
                           new TH1D(Form("pData_MG_hEx_M%i_Z%i_%s", it_M, it_Z, particle.c_str()),
                                    Form("Excitation energy with M%i Z%i in VAMOS and %s in MUGAST", it_M, it_Z, particle.c_str()),
                                    1000, -60, 60));
                Istantiate(pData.MG.mEx_EDC[it_M][it_Z][particle],
                           new TH2D(Form("pData_MG_mEx_EDC_M%i_Z%i_%s", it_M, it_Z, particle.c_str()),
                                    Form("Ex vs EDC with M%i Z%i in VAMOS and %s in MUGAST", it_M, it_Z, particle.c_str()),
                                    1000, -10, 10, 2000, 0, 2000));
                for (const auto &it_gamma : gammas)
                {
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

    //TODO: Add MUST2

    if (new_graph_file != nullptr)
    {
        new_graph_file->close();
        //TODO: kill all threads here
        exit(1);
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

//Processing entry//////////////////////////////////////////////////////////////////////
Bool_t Selector::Process(Long64_t entry)
{
#ifdef VERBOSE_DEBUG
    std::cout << "------------>Selector::Process()\n";
#endif
    fReader.SetLocalEntry(entry);
    if (entry % 5000 == 0)
        std::cout << "Processed entries : " << entry << " of  " << total_entries << "\n";


    //Vamos/////////////////////////////////////////////
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

    //Agata/////////////////////////////////////////////
    LoadAgataData();
#ifdef VERBOSE_DEBUG
    std::cout << "------------>Finished: Loading Agata data, positive exit\n";
#endif

    PlotAgataGraphs();
#ifdef VERBOSE_DEBUG
    std::cout << "------------>Finished: Plotting agata graphs, positive exit\n";
#endif

    //MUGAST////////////////////////////////////////////
    LoadMugastData();
#ifdef VERBOSE_DEBUG
    std::cout << "------------>Finished: Loading Mugast data\n";
#endif

    mugast_fragment.Identify();
#ifdef VERBOSE_DEBUG
    std::cout << "------------>Finished: Mugast identification, negative exit\n";
#endif

    PlotMugastGraphs();
#ifdef VERBOSE_DEBUG
    std::cout << "------------>Finished: Analysis of one event\n";
#endif

    if (entry % 10000 == 0)
        mugast_fragment.StoreTWvsIce();

    return kTRUE;
}

void Selector::SlaveTerminate()
{
#ifdef VERBOSE_DEBUG
    std::cout << "------------>Selector::SlaveTerminate()\n";
#endif
    TFile *top = new TFile(file_name.c_str(), "recreate");
    std::cout << "Output file : " << file_name << "\n";
    TIter iter(fOutput);
    TObject *obj;
    while ((obj = iter()))
    {
        obj->Write();
    }
    top->Close();
}

void Selector::Terminate()
{
#ifdef VERBOSE_DEBUG
    std::cout << "------------>Selector::Terminate()\n";
#endif
}

inline void Selector::LoadVamosData()
{
    vamos_fragment.SetData(new VamosIdentification::Data(&IC, &Path, &Brho, &Xf,
                                                         &ThetaL, &PhiL, &AGAVA_VAMOSTS,
                                                         &T_FPMW_CATS2_C));
}

inline void Selector::LoadMugastData()
{
    mugast_fragment.SetData(new MugastIdentification::Data(&Mugast,
                                                           &TW,
                                                           vamos_fragment.Get_id_M(),
                                                           vamos_fragment.Get_id_Z()));
}

inline void Selector::LoadAgataData()
{
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

inline void Selector::PlotVamosGraphs()
{
#ifdef VERBOSE_DEBUG
    std::cout << "------------>Selector::PlotVamosGraphs()\n";
#endif
    //dE-E plot, no conditions
    Fill(pConf.VAMOS.mdE_E,
         vamos_fragment.Get_En(), vamos_fragment.Get_D_En());
    Fill(pConf.VAMOS.mdE2_E,
         vamos_fragment.Get_En(), vamos_fragment.Get_D_En2());

    if (vamos_fragment.Get_id_Z() == 0)
        return;
    Fill(pConf.VAMOS.mQ_MQ[vamos_fragment.Get_id_Z()],
         vamos_fragment.Get_M_Q(), vamos_fragment.Get_Charge());
    Fill(pConf.VAMOS.Xf_MQ[vamos_fragment.Get_id_Z()],
         vamos_fragment.Get_M_Q(), *Xf);
    if (!vamos_fragment.Identified())
        return;
    Fill(pData.VAMOS.mTW_Brho[vamos_fragment.Get_id_M()][vamos_fragment.Get_id_Z()],
         *TW, *Brho);
}

inline void Selector::PlotMugastGraphs()
{
#ifdef VERBOSE_DEBUG
    std::cout << "------------>Selector::PlotMugastGraphs()\n";
#endif
    for (int ii = 0; ii < mugast_fragment.Get_Mult(); ++ii)
    {
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

        if (agata_gammas.In_Coincidence())
        {
            for (long unsigned int jj = 0; jj < agata_gammas.Get_Mult(); ++jj)
            {
                Fill(pData.MG.mEx_EDC[vamos_fragment.Get_id_M()]
                                     [vamos_fragment.Get_id_Z()]
                                     [mugast_fragment.Get_Particle(ii)],
                     mugast_fragment.Get_Ex(ii),
                     agata_gammas.Get_EDC(jj));
                for (const auto &gamma : gammas)
                {
                    if (gamma == "NOCONDITION")
                        continue;
                    if (abs(std::stod(gamma) - agata_gammas.Get_EDC(jj)) > 20)
                        continue;

                    Fill(pData.MG.mELab_ThetaLab[vamos_fragment.Get_id_M()]
                                                [vamos_fragment.Get_id_Z()]
                                                [mugast_fragment.Get_Particle(ii)]
                                                [gamma],
                         mugast_fragment.Get_ThetaLab(ii),
                         mugast_fragment.Get_E(ii));
                }
            }
        }

        //E TOF
        Fill(pConf.MG.mE_TOF[mugast_fragment.Get_MG(ii)],
             mugast_fragment.Get_SI_E(ii),
             mugast_fragment.Get_T(ii));

        Fill(pConf.MG.mE_TOF2[mugast_fragment.Get_MG(ii)],
             mugast_fragment.Get_SI_E(ii),
             mugast_fragment.Get_T2(ii));
        //Strip
        //E

        for (const auto &dimension : mugast_fragment.strips)
        {
            Fill(pConf.MG.mStrip_E[mugast_fragment.Get_MG(ii)][dimension],
                 mugast_fragment.Get_SI_X(ii),
                 mugast_fragment.Get_SI_E(ii));

            Fill(pConf.MG.mStrip_T[mugast_fragment.Get_MG(ii)][dimension],
                 mugast_fragment.Get_SI_X(ii),
                 mugast_fragment.Get_T(ii));

            Fill(pConf.MG.mStrip_T2[mugast_fragment.Get_MG(ii)][dimension],
                 mugast_fragment.Get_SI_X(ii),
                 mugast_fragment.Get_T2(ii));
        }

        if (vamos_fragment.Identified())
        {
            Fill(pData.MG.mE_TOF[vamos_fragment.Get_id_M()][vamos_fragment.Get_id_Z()][mugast_fragment.Get_MG(ii)],
                 mugast_fragment.Get_SI_E(ii),
                 mugast_fragment.Get_T(ii));

            Fill(pData.MG.mE_TOF2[vamos_fragment.Get_id_M()][vamos_fragment.Get_id_Z()][mugast_fragment.Get_MG(ii)],
                 mugast_fragment.Get_SI_E(ii),
                 mugast_fragment.Get_T2(ii));
        }
    }
}

inline void Selector::PlotAgataGraphs()
{
    for (int ii = 0; ii < *nbAdd; ii++)
    {
        Fill(pConf.AGATA.mmAGATA3D, AddX[ii], AddY[ii], AddZ[ii]);
    }
    if (agata_gammas.In_Coincidence())
    {
        for (unsigned int ii = 0; ii < agata_gammas.Get_Mult(); ++ii)
        {
            Fill(pData.AGATA.hDC[vamos_fragment.Get_id_M()][vamos_fragment.Get_id_Z()]["NONE"],
                 agata_gammas.Get_EDC(ii));
            for (unsigned int jj = 0; jj < agata_gammas.Get_Mult(); ++jj)
            {
                if (ii == jj)
                    continue;

                Fill(pData.AGATA.mDC[vamos_fragment.Get_id_M()][vamos_fragment.Get_id_Z()],
                     agata_gammas.Get_EDC(ii),
                     agata_gammas.Get_EDC(jj));
            }
        }
    }
}

bool Selector::GetSettings()
{
    std::string config_file = "./Configs/GraphsEnabled.txt";
    std::ifstream ifile(config_file);
    if (ifile)
    {
        std::string line;
        while (std::getline(ifile, line))
        {
            std::istringstream str(line);
            std::string Graph;
            str >> Graph;
            std::string enabled;
            str >> enabled;
            if (enabled.compare("false"))
            { //IS ENABLED
                enabled_histograms[Graph] = true;
            }
        }
    }
    else
    {
        ifile.close();
        new_graph_file = new std::ofstream(config_file);
    }
    return true;
}

std::vector<std::pair<double, double>> Selector::GetTWvsIce()
{
    return mugast_fragment.GetTWvsIce();
}

inline bool Selector::Fill(TH1 *histo, const double &data1)
{
    if (histo == nullptr)
    {
        return false;
    }
    else
    {
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

inline bool Selector::Fill(TH2 *histo,
                           const double &data1, const double &data2)
{
    if (histo == nullptr)
    {
        return false;
    }
    else
    {
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

inline bool Selector::Fill(TH3 *histo, const double &data1,
                           const double &data2, const double &data3)
{
    if (histo == nullptr)
    {
        return false;
    }
    else
    {
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
