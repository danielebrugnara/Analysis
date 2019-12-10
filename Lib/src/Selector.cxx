#define Selector_cxx
// The class definition in Selector.h has been generated automatically
// by the ROOT utility TTree::MakeSelector(). This class is derived
// from the ROOT class TSelector. For more information on the TSelector
// framework see $ROOTSYS/README/README.SELECTOR or the ROOT User Manual.


// The following methods are defined in this file:
//    Begin():        called every time a loop on the tree starts,
//                    a convenient place to create your histograms.
//    SlaveBegin():   called after Begin(), when on PROOF called only on the
//                    slave servers.
//    Process():      called for each event, in this function you decide what
//                    to read and fill your histograms.
//    SlaveTerminate: called at the end of the loop on the tree, when on PROOF
//                    called only on the slave servers.
//    Terminate():    called at the end of the loop on the tree,
//                    a convenient place to draw/fit your histograms.
//
// To use this file, try the following session on your Tree T:
//
// root> T->Process("Selector.C")
// root> T->Process("Selector.C","some options")
// root> T->Process("Selector.C+")
//


#include "Selector.h"

void Selector::Begin(TTree * /*tree*/)
{
   // The Begin() function is called at the start of the query.
   // When running with PROOF Begin() is only called on the client.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();
}

void Selector::SlaveBegin(TTree * /*tree*/)
{
   // The SlaveBegin() function is called after the Begin() function.
   // When running with PROOF SlaveBegin() is called on each slave server.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();

   //   Long64_t Start = 1E10;
   //   Long64_t Stop = 5E13;

   //Initializing Histograms

   //Configurations
   //Target
   thickness_angle = new Interpolation("thickess.txt");
   angle_angle = new Interpolation("angle.txt");

   //VAMOS
   pConf.VAMOS.mdE_E = new TH2D("pConf-VAMOS-mdE_E", "dE E in VAMOS", 400, 0, 400, 400, 0, 400);
   fOutput->Add(pConf.VAMOS.mdE_E);

   pConf.VAMOS.mQ_MQ["Ar"] = new TH2D("pConf-VAMOS-mQ_MQ-Ar", "M/Q vs Q with Ar selection", 1000, 2, 4, 1000, 10, 24);
   fOutput->Add(pConf.VAMOS.mQ_MQ["Ar"]);
   pConf.VAMOS.mQ_MQ["K"] = new TH2D("pConf-VAMOS-mQ_MQ-K", "M/Q vs Q with K selection", 1000, 2, 4, 1000, 10, 24);
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
   for (int ii = 0; ii < 128; ii++)
   {
      strips.push_back(to_string(ii));
   }

   for (const auto &MM : siliconsMM)
   {
      pConf.SI.mdE_E_Si[MM] = new TH2D(Form("pConf-SI-mdE_E_Si-%s", MM.c_str()), Form("dE E of %s", MM.c_str()), 1000, 0, 28, 1000, 0, 28);
      fOutput->Add(pConf.SI.mdE_E_Si[MM]);
   }
   for (const auto &SI : silicons)
   {
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

   //Data
   for (const auto &nucleus : nuclei)
   {
      for (const auto &mass : masses)
      {
         pData.VAMOS.mTW_Brho[mass][nucleus] = new TH2D(Form("pData-VAMOS-mTW_Brho-%s-%s", mass.c_str(), nucleus.c_str()), Form("Time vs Brho with %s %s in VAMOS", mass.c_str(), nucleus.c_str()), 5000, 242, 328, 1000, 0.5, 1.5);
         fOutput->Add(pData.VAMOS.mTW_Brho[mass][nucleus]);
      }
   }
   //VAMOS

   //AGATA
   for (const auto &nucleus : nuclei)
   {
      for (const auto &mass : masses)
      {
         for (const auto &condition : AGATAconditions)
         {
            pData.AGATA.hDC[mass][nucleus][condition] = new TH1D(Form("pData-AGATA-hDC-%s-%s-cond%s", mass.c_str(), nucleus.c_str(), condition.c_str()), Form("DC spectrum of %s %s with condition %s", mass.c_str(), nucleus.c_str(), condition.c_str()), 4000, 0, 4000);
            fOutput->Add(pData.AGATA.hDC[mass][nucleus][condition]);
         }
         pData.AGATA.mDC[mass][nucleus] = new TH2D(Form("pData-AGATA-mDC-%s-%s", mass.c_str(), nucleus.c_str()), Form("DC gamma gamma of %s %s", mass.c_str(), nucleus.c_str()), 4000, 0, 4000, 4000, 0, 4000);
         fOutput->Add(pData.AGATA.mDC[mass][nucleus]);
         pData.AGATA.mDC_ThetaMUGAST[mass][nucleus] = new TH2D(Form("pData-AGATA-mDC_ThetaMUGAST-%s-%s", mass.c_str(), nucleus.c_str()), Form("DC gamma vs Theta on MUGAST of %s %s", mass.c_str(), nucleus.c_str()), 4000, 0, 4000, 180, 0, 180);
         fOutput->Add(pData.AGATA.mDC_ThetaMUGAST[mass][nucleus]);
         for(const auto & particle: particles){
            pData.AGATA.mEx_DC[mass][nucleus][particle] = new TH2D(Form("pData-AGATA-mEx_DC-%s-%s-%s", mass.c_str(), nucleus.c_str(), particle.c_str()), Form("Excitation energy AGATA vs MUGAST %s %s and %s", mass.c_str(), nucleus.c_str(), particle.c_str()), 1000, 0, 10, 1000, 0, 10);
            fOutput->Add(pData.AGATA.mEx_DC[mass][nucleus][particle]);
            pData.AGATA.mELab_ThetaLab[mass][nucleus][particle] = new TH2D(Form("pData-AGATA-mELab_ThetaLab-%s-%s-%s", mass.c_str(), nucleus.c_str(), particle.c_str()), Form("Excitation energy AGATA vs MUGAST %s %s and %s", mass.c_str(), nucleus.c_str(), particle.c_str()), 1000, 0, 10, 1000, 0, 10);
            fOutput->Add(pData.AGATA.mELab_ThetaLab[mass][nucleus][particle]);
         }
      }
   }

   //CATS

   //MUGAST
   for (const auto &nucleus : nuclei)
   {
      for (const auto &mass : masses)
      {
         for (const auto &MM : siliconsMM)
         {
            pData.SI.mdE_E_Si[MM][mass][nucleus] = new TH2D(Form("pData-SI-mdE_E_Si-%s-%s-%s", MM.c_str(), mass.c_str(), nucleus.c_str()), Form("dE E of %s with %s %s", MM.c_str(), mass.c_str(), nucleus.c_str()), 1000, 0, 28, 1000, 0, 28);
            fOutput->Add(pData.SI.mdE_E_Si[MM][mass][nucleus]);
         }
         for (const auto &SI : silicons)
         {
            pData.SI.mE_TOF[SI][mass][nucleus] = new TH2D(Form("pData-SI-mE_TOF-%s-%s-%s", SI.c_str(), mass.c_str(), nucleus.c_str()), Form("E vs TOF of %s with %s %s", SI.c_str(), mass.c_str(), nucleus.c_str()), 1000, 0, 28, 1000, 260, 380);
            fOutput->Add(pData.SI.mE_TOF[SI][mass][nucleus]);
         }
         pData.SI.mE_TOF["MG"][mass][nucleus] = new TH2D(Form("pData-SI-mE_TOF-MG-%s-%s", mass.c_str(), nucleus.c_str()), Form("E vs TOF of all MG with %s %s", mass.c_str(), nucleus.c_str()), 1000, 0, 28, 1000, 260, 380);
         fOutput->Add(pData.SI.mE_TOF["MG"][mass][nucleus]);
         for (const auto &particle : particles)
         {
            pData.SI.hEx[mass][nucleus][particle] = new TH1D(Form("pData-SI-hEx-%s-%s-%s", mass.c_str(), nucleus.c_str(), particle.c_str()), Form("Excitation energy with %s %s in VAMOS and %s in MUGAST", mass.c_str(), nucleus.c_str(), particle.c_str()), 1000, -60, 60);
            fOutput->Add(pData.SI.hEx[mass][nucleus][particle]);
            pData.SI.mEx_TW[mass][nucleus][particle] = new TH2D(Form("pData-SI-mEx_TW-%s-%s-%s", mass.c_str(), nucleus.c_str(), particle.c_str()), Form("Excitation energy vs Time with %s %s in VAMOS and %s in MUGAST", mass.c_str(), nucleus.c_str(), particle.c_str()), 5000,  242, 328, 1000, -60, 60);
            fOutput->Add(pData.SI.mEx_TW[mass][nucleus][particle]);
            pData.SI.mECM_ThetaCM[mass][nucleus][particle] = new TH2D(Form("pData-SI-mECM_ThetaCM-%s-%s-%s", mass.c_str(), nucleus.c_str(), particle.c_str()), Form("E CM vs Theta CM with %s %s in VAMOS and %s in MUGAST", mass.c_str(), nucleus.c_str(), particle.c_str()), 1000, 0, 180, 1000, 0, 60);
            fOutput->Add(pData.SI.mECM_ThetaCM[mass][nucleus][particle]);
            for (const auto &gamma: gammas)
            {
               pData.SI.mELab_ThetaLab[mass][nucleus][particle][gamma] = new TH2D(Form("pData-SI-mELab_ThetaLab-%s-%s-%s-%s", mass.c_str(), nucleus.c_str(), particle.c_str(), gamma.c_str()), Form("ELab vs Theta Lab with %s %s in VAMOS and %s in MUGAST and %s in AGATA", mass.c_str(), nucleus.c_str(), particle.c_str(), gamma.c_str()), 1000, 0, 180, 1000, 0, 60);
               fOutput->Add(pData.SI.mELab_ThetaLab[mass][nucleus][particle][gamma]);
            }
         }
      }
   }
for (const auto & particle: particles){
   pData.SI.mELab_ThetaLab["ANY"]["ANY"][particle]["ANY"] = new TH2D(Form("pData-SI-mELab_ThetaLab-%s-%s-%s", "ANY", "ANY", particle.c_str()), Form("E Lab vs Theta Lab with %s %s in VAMOS and %s in MUGAST", "ANY", "ANY", particle.c_str()), 1000, 0, 180, 1000, 0, 60);
   fOutput->Add(pData.SI.mELab_ThetaLab["ANY"]["ANY"][particle]["ANY"]);
}



   //Loading graphical cuts
//   struct stat buffer;
//   if (stat("./cuts/VAMOS.root", &buffer) == 0)
  // if (1)
//   std::;

 //  if(std::filesystem::exists(vamos_cuts))

   std::ifstream ifile("./cuts/VAMOS.root");
   if(ifile)
   {
      ifile.close();
      VAMOScuts = new TFile("./cuts/VAMOS.root", "READ");
      if (!(VAMOScuts->IsOpen()))
      {
         std::cout << "VAMOS file not opened\n";
      }
      else
      {
         std::cout << "VAMOS file opened\n";
         TIter contents(VAMOScuts->GetListOfKeys());
         TKey *key;
         TObject *obj;
         while ((key = (TKey *)contents()))
         {
            obj = VAMOScuts->Get(key->GetName());
            //         TClass *cl = gROOT->GetClass(key->GetClassName());
            //         if (cl->InheritsFrom("TCutG")){
            if (obj->InheritsFrom("TCutG"))
            {
               TCutG *tmp = (TCutG *)obj;
               cut["VAMOS"][tmp->GetName()] = tmp;
               std::cout << "Found cut in VAMOS :" << tmp->GetName() << std::endl;
            }
         }
      }
   }

   ifile.close();
   ifile.open("./cuts/MUGAST_tmp.root");
   if (ifile)
   {
      ifile.close();
      MUGASTcuts = new TFile("./cuts/MUGAST_tmp.root", "READ");
      if (!(MUGASTcuts->IsOpen()))
      {
         std::cout << "MUGAST file not opened\n";
      }
      else
      {
         std::cout << "MUGAST file opened\n";
         TIter contents(MUGASTcuts->GetListOfKeys());
         TKey *key;
         TObject *obj;
         while ((key = (TKey *)contents()))
         {
            obj = MUGASTcuts->Get(key->GetName());
            //         TClass *cl = gROOT->GetClass(key->GetClassName());
            //         if (cl->InheritsFrom("TCutG")){
            if (obj->InheritsFrom("TCutG"))
            {
               TCutG *tmp = (TCutG *)obj;
               cut["MUGAST"][tmp->GetName()] = tmp;
               std::cout << "Found cut in MUGAST :" << tmp->GetName() << std::endl;
            }
         }
      }
   }

   mass["46Ar"] = 45.968082712 * AMU_TO_MEV; //in MeV
   mass["47Ar"] = 46.972934865 * AMU_TO_MEV; //in MeV
   mass["47K"] = 46.961661614 * AMU_TO_MEV;  //in MeV
   mass["46K"] = 45.961981586 * AMU_TO_MEV;
   mass["d"] = 2.01410177812 * AMU_TO_MEV;
   mass["1H"] = 1.00782503223 * AMU_TO_MEV;
   //   if (stat("./cuts/MUGAST.root", &buffer)==0);
   //   if (stat("./cuts/AGATA.root", &buffer)==0);
   //   if (stat("./cuts/CATS.root", &buffer)==0);

   counter = 0;
   count46Ar = 0;

   GetSettings();
}

Bool_t Selector::Process(Long64_t entry)
{
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

   counter++;

   Double_t mE2 = 0;
   Double_t mdE2 = 0;
   Double_t mdE3 = 0;
   Double_t mT2 = 0;
   Double_t mD2 = 0;
   Double_t mV2 = 0;
   Double_t mBeta2 = 0;
   Double_t mGamma2 = 0;
   Double_t mM2 = 0;
   Double_t mM_Q2 = 0;
   Double_t mQ2 = 0;

   Bool_t VAMOS_GOOD = true;
   Bool_t AGATA_GOOD = *AddTS - *LTS > 175 && *AddTS - *LTS < 184;
   IdentifiedNucleus.found = false;

   //Initialization of VAMOS variables

   //BASIC FOM DATA
   if (IC.GetSize() > 0)
   {
      //mE2 = 0.85*(pow(*Xf+450*tan(*Tf/1000.),1)*0.0000015*IC[0] + 2.2*IC[0]+1.1*IC[1]*(IC[0]>0.1)+1.34*1.95*IC[2]*(IC[0]>0.1)*(IC[1]>0.1)+3.1*IC[3]*(IC[0]>0.1)*(IC[1]>0.1)*(IC[2]>0.1)+1.05*IC[4]*(IC[0]>0.1)*(IC[1]>0.1)*(IC[2]>0.1)*(IC[3]>0.1)+IC[5]*(IC[0]>0.1)*(IC[1]>0.1)*(IC[2]>0.1)*(IC[3]>0.1)*(IC[4]>0.1));
      //mE2 = 0.85*(pow(*Xf+351*tan(*Tf/1000.),1)*0.0000015*IC[0] +  2.2*IC[0]+1.1*IC[1]*(IC[0]>0.1)+1.34*1.95*IC[2]*(IC[0]>0.1)*(IC[1]>0.1)+3.1*IC[3]*(IC[0]>0.1)*(IC[1]>0.1)*(IC[2]>0.1)+1.05*IC[4]*(IC[0]>0.1)*(IC[1]>0.1)*(IC[2]>0.1)*(IC[3]>0.1)+IC[5]*(IC[0]>0.1)*(IC[1]>0.1)*(IC[2]>0.1)*(IC[3]>0.1)*(IC[4]>0.1));
      mE2 = 2.2 * IC[0] + 1.1 * IC[1] * (IC[0] > 0.1) + 1.34 * 1.95 * IC[2] * (IC[0] > 0.1) * (IC[1] > 0.1) + 3.1 * IC[3] * (IC[0] > 0.1) * (IC[1] > 0.1) * (IC[2] > 0.1) + 1.05 * IC[4] * (IC[0] > 0.1) * (IC[1] > 0.1) * (IC[2] > 0.1) * (IC[3] > 0.1) + IC[5] * (IC[0] > 0.1) * (IC[1] > 0.1) * (IC[2] > 0.1) * (IC[3] > 0.1) * (IC[4] > 0.1);
      //mdE2 = IC[0]+IC[1]*(IC[0]>0.1);
      //mdE2 = pow(*Xf+450*tan(*Tf/1000.),1)*0.0000015*IC[0] + 2.2*IC[0];
      mdE2 = 2.2 * IC[0] + 1.1 * IC[1] * (IC[0] > 0.1);
      mdE3 = IC[0];
   }
   else
      VAMOS_GOOD = false;
   if (MW_N.GetSize() == 1)
      mT2 = 540.5 * (*AGAVA_VAMOSTS < 104753375647998) + 537.9 * (*AGAVA_VAMOSTS >= 104753375647998) - 2. * (*T_FPMW_CATS2_C) + 2.7 * (MW_N[0] == 16) + 2.7 * (MW_N[0] == 15) + 2.9 * (MW_N[0] == 14) + 2.9 * (MW_N[0] == 13) + 2.4 * (MW_N[0] == 12) + 1.3 * (MW_N[0] == 11) + 1.5 * (MW_N[0] == 10) + 1.6 * (MW_N[0] == 9) - 0.6 * (MW_N[0] == 8) + 2.5 * (MW_N[0] == 7) + 2. * (MW_N[0] == 6) + 1.6 * (MW_N[0] == 5) + 1.1 * (MW_N[0] == 4) - 0.6 * (MW_N[0] == 3) - 1.2 * (MW_N[0] == 2) - 4.0 * (MW_N[0] == 1);
   else
      VAMOS_GOOD = false;
   mD2 = *Path + 5;

   if (VAMOS_GOOD)
   {

      //STRUCTURED
      mV2 = mD2 / mT2;
      mBeta2 = mV2 / 29.9792;
      mGamma2 = 1. / sqrt(1.0 - mBeta2 * mBeta2);
      mM2 = (mE2) / 931.5016 / (mGamma2 - 1.);
      //mM2 = 18./20.8*(mE2)/931.5016/(mGamma2-1.);
      mM_Q2 = *Brho / 3.105 / mBeta2 / mGamma2;
      mQ2 = mM2 / mM_Q2;
   }

   //Start of analysis
   //Filling historgrams
   if (VAMOS_GOOD)
   {
      if (abs(*Xf) < 300 && *EWIRE_1_2 >8000 && *EWIRE_1_2 < 16000 && *EWIRE_1_1 > 8000 && *EWIRE_1_1< 16000 && *EWIRE_2_1>10000 && *EWIRE_2_1<18000 && *EWIRE_2_2>4000 && *EWIRE_2_2<7000)
         Fill(pConf.VAMOS.mdE_E, mE2, mdE2);
      try
      {
         //46Ar
         TLorentzVector p4(0, 0, 0, mass["46Ar"]);
         TVector3 b4(0, 0, mBeta2);
         b4.SetMagThetaPhi(mBeta2, *ThetaL, *PhiL);
         //b4.SetMagThetaPhi(mBeta2, 0, 0);
         //b4.SetMagThetaPhi(0, 0, 0);
         p4.Boost(b4);
         for (int ii = 0; ii < (*CATS).PositionX.size(); ii++)
         {
            Fill(pConf.CATS.mCATSpos, (*CATS).PositionX[ii], (*CATS).PositionY[ii]);
         }
         if (cut.at("VAMOS").at("dE_E_Ar")->IsInside(mE2, mdE2))
         {
            Fill(pConf.VAMOS.mQ_MQ["Ar"], mM_Q2, mQ2);
            for (const auto &Qcut : QcutsAr)
            {
               //for (const auto &Qcut : Qcuts)
               //{
               if (cut.at("VAMOS").at(Qcut)->IsInside(mM_Q2, mQ2))
               {
                  count46Ar++;
                  Fill(pConf.VAMOS.hAmass["Ar"], mM_Q2 * stoi(Qcut.substr(1, 2)));
                  //if (mM_Q2 * stoi(Qcut.substr(1, 2)) < 46.5 && mM_Q2 * stoi(Qcut.substr(1, 2)) > 45.5)
                  if (1)
                  {
                     IdentifiedNucleus.mass = "46";
                     IdentifiedNucleus.name = "Ar";
                     IdentifiedNucleus.found = true;
                     //if (*AddTS - *LTS > 175 && *AddTS - *LTS < 184)
                     if(AGATA_GOOD)
                     {
                        for (int ii = 0; ii < AddE.GetSize(); ii++)
                        {
                           if (AddE[ii] > 10)
                           {
                              Fill(pData.AGATA.hDC["46"]["Ar"]["NONE"], 1E3 * CorrectDoppler(p4, AddE[ii] / 1E3, AddX[ii], AddY[ii], AddZ[ii]));
                              for (int kk = 0; kk < (*Mugast).PosX.size(); kk++)
                              {
                                 TVector3 vec((*Mugast).PosX[kk], (*Mugast).PosY[kk], (*Mugast).PosZ[kk]);
                                 Fill(pData.AGATA.mDC_ThetaMUGAST["46"]["Ar"], 1E3 * CorrectDoppler(p4, AddE[ii] / 1E3, AddX[ii], AddY[ii], AddZ[ii]), vec.Theta() * TMath::RadToDeg());
                              }
                              //Gamma Gamma matrices
                              for (int jj = 0; jj < AddE.GetSize(); jj++)
                              {
                                 if (AddE[jj] > 10 && ii != jj)
                                 {
                                    Fill(pData.AGATA.mDC["46"]["Ar"], 1E3 * CorrectDoppler(p4, AddE[ii] / 1E3, AddX[ii], AddY[ii], AddZ[ii]), 1E3 * CorrectDoppler(p4, AddE[jj] / 1E3, AddX[jj], AddY[jj], AddZ[jj]));
                                    Fill(pData.AGATA.mDC["46"]["Ar"], 1E3 * CorrectDoppler(p4, AddE[jj] / 1E3, AddX[jj], AddY[jj], AddZ[jj]), 1E3 * CorrectDoppler(p4, AddE[ii] / 1E3, AddX[ii], AddY[ii], AddZ[ii]));
                                 }
                              }
                           }
                        }
                     }
                  }
               }
            }
         }
      }
      catch (std::out_of_range &e)
      {
         std::cerr << e.what() << std::endl;
      }

      try
      {
         //47K
         TLorentzVector p4(0, 0, 0, mass["47K"]);
         TVector3 b4(0, 0, mBeta2);
         b4.SetMagThetaPhi(mBeta2, *ThetaL, *PhiL);
         //b4.SetMagThetaPhi(mBeta2, 0, 0);
         p4.Boost(b4);
         if (cut.at("VAMOS").at("dE_E_K")->IsInside(mE2, mdE2))
         {
            Fill(pConf.VAMOS.mQ_MQ["K"], mM_Q2, mQ2);
            //for (const auto & Qcut: Qcuts){
            for (const auto &Qcut : QcutsK)
            {
               if (cut.at("VAMOS").at(Qcut)->IsInside(mM_Q2, mQ2))
               {
                  Fill(pConf.VAMOS.hAmass["K"], mM_Q2 * stoi(Qcut.substr(1, 2)));
                  IdentifiedNucleus.mass = "47";
                  IdentifiedNucleus.name = "K";
                  IdentifiedNucleus.found = true;
                  if (AGATA_GOOD && *GATCONF_MASTER == 1)//GATCONF_MASTER==1 =>MUGAST trigger
                  {
                     for (int ii = 0; ii < AddE.GetSize(); ii++)
                     {
                        if (AddE[ii] > 10)
                        {
                           Fill(pData.AGATA.hDC["47"]["K"]["NONE"], 1E3 * CorrectDoppler(p4, AddE[ii] / 1E3, AddX[ii], AddY[ii], AddZ[ii]));
                           for (int kk = 0; kk < (*Mugast).PosX.size(); kk++)
                           {
                              TVector3 vec((*Mugast).PosX[kk], (*Mugast).PosY[kk], (*Mugast).PosZ[kk]);
                              Fill(pData.AGATA.mDC_ThetaMUGAST["47"]["K"], 1E3 * CorrectDoppler(p4, AddE[ii] / 1E3, AddX[ii], AddY[ii], AddZ[ii]), vec.Theta() * TMath::RadToDeg());
                           }
                           for (int jj = 0; jj < AddE.GetSize(); jj++)
                           {
                              if (AddE[jj] > 10 && ii != jj)
                              {
                                 Fill(pData.AGATA.mDC["47"]["K"], 1E3 * CorrectDoppler(p4, AddE[ii] / 1E3, AddX[ii], AddY[ii], AddZ[ii]), 1E3 * CorrectDoppler(p4, AddE[jj] / 1E3, AddX[jj], AddY[jj], AddZ[jj]));
                                 Fill(pData.AGATA.mDC["47"]["K"], 1E3 * CorrectDoppler(p4, AddE[jj] / 1E3, AddX[jj], AddY[jj], AddZ[jj]), 1E3 * CorrectDoppler(p4, AddE[ii] / 1E3, AddX[ii], AddY[ii], AddZ[ii]));
                              }
                           }
                        }
                     }
                  }
                  //}
               }
            }
         }
      }
      catch (std::out_of_range &e)
      {
         std::cerr << e.what() << std::endl;
      }
   }
   //VAMOS Data
   if (IdentifiedNucleus.found)
   {
      Fill(pData.VAMOS.mTW_Brho[IdentifiedNucleus.mass][IdentifiedNucleus.name], *TW, *Brho);
   }

   //MUGAST

   for (int ii = 0; ii < (*Mugast).DSSD_E.size(); ii++)
   {
      //(XE)
      Fill(pConf.SI.mStrip_E[Form("MG%i", (*Mugast).TelescopeNumber[ii])]["X"], (*Mugast).DSSD_X[ii], (*Mugast).DSSD_E[ii]);
      //(YE)
      Fill(pConf.SI.mStrip_E[Form("MG%i", (*Mugast).TelescopeNumber[ii])]["Y"], (*Mugast).DSSD_Y[ii], (*Mugast).DSSD_E[ii]);
   }
   for (int ii = 0; ii < (*Mugast).DSSD_T.size(); ii++)
   {
      //(TE)
      Fill(pConf.SI.mStrip_T[Form("MG%i", (*Mugast).TelescopeNumber[ii])]["X"], (*Mugast).DSSD_X[ii], (*Mugast).DSSD_T[ii]);
      //(TE)
      Fill(pConf.SI.mStrip_T[Form("MG%i", (*Mugast).TelescopeNumber[ii])]["Y"], (*Mugast).DSSD_Y[ii], (*Mugast).DSSD_T[ii]);
   }


   //SI data loops
   try
   {
      //E-TOF plots
      for (int ii = 0; ii < (*Mugast).DSSD_E.size(); ii++)
      {
         TVector3 hitPos((*Mugast).PosX[ii], (*Mugast).PosY[ii], (*Mugast).PosZ[ii]);
         Fill(pConf.SI.mE_TOF[Form("MG%i", (*Mugast).TelescopeNumber[ii])], (*Mugast).DSSD_E[ii], AlignT((*Mugast).TelescopeNumber[ii], (*Mugast).DSSD_Y[ii], (*Mugast).DSSD_T[ii]));
         if (IdentifiedNucleus.found)
         {
            Fill(pData.SI.mE_TOF[Form("MG%i", (*Mugast).TelescopeNumber[ii])][IdentifiedNucleus.mass][IdentifiedNucleus.name], (*Mugast).DSSD_E[ii], AlignT((*Mugast).TelescopeNumber[ii], (*Mugast).DSSD_Y[ii], (*Mugast).DSSD_T[ii]));
            Fill(pData.SI.mE_TOF["MG"][IdentifiedNucleus.mass][IdentifiedNucleus.name], (*Mugast).DSSD_E[ii], AlignPunch((*Mugast).TelescopeNumber[ii], (*Mugast).DSSD_T[ii]));
         }
      }

      //Loop on physics data
      int MugastEvents = 0;
      for (int ii = 0; ii < DetID.GetSize(); ii++)
      {
         Fill(pData.SI.mELab_ThetaLab["ANY"]["ANY"]["ANY"]["ANY"], ThetaLab[MugastEvents], ELab[MugastEvents]);
         if (DetID[ii] >= 100)//These events are in Must2
            continue;
         if (DetID[ii] != (*Mugast).TelescopeNumber[MugastEvents])
            std::cout << "Something is wrong matching Mugast events\n";
         //Excitation energy and kinematic lines

         Fill(pData.SI.hEx[IdentifiedNucleus.mass][IdentifiedNucleus.name]["ANY"], Ex[MugastEvents]);
         Fill(pData.SI.mEx_TW[IdentifiedNucleus.mass][IdentifiedNucleus.name]["ANY"], *TW, Ex[MugastEvents]);
         Fill(pData.SI.mELab_ThetaLab[IdentifiedNucleus.mass][IdentifiedNucleus.name]["ANY"]["ANY"], ThetaLab[MugastEvents], ELab[MugastEvents]);
         Fill(pData.SI.mECM_ThetaCM[IdentifiedNucleus.mass][IdentifiedNucleus.name]["ANY"], ThetaCM[MugastEvents], Ecm[MugastEvents]);

         //Beta reconstruction from VAMOS,  p4 is the heavy fragment
         TLorentzVector p4(0, 0, 0, mass["46Ar"]);
         TVector3 b4(0, 0, mBeta2);
         b4.SetMagThetaPhi(mBeta2, *ThetaL, *PhiL);
         p4.Boost(b4);


         if (AGATA_GOOD)
         {
            //Loop over gammas
            for (int ii = 0; ii < AddE.GetSize(); ii++)
            {
               if (AddE[ii] > 10)
               {
                  Fill(pData.AGATA.mEx_DC[IdentifiedNucleus.mass][IdentifiedNucleus.name]["ANY"], Ex[MugastEvents], CorrectDoppler(p4, AddE[ii] / 1E3, AddX[ii], AddY[ii], AddZ[ii]));
               }
               if (AddE[ii] > 320 && AddE[ii]<390){
                  Fill(pData.SI.mELab_ThetaLab[IdentifiedNucleus.mass][IdentifiedNucleus.name]["ANY"]["360 keV"], ThetaLab[MugastEvents], ELab[MugastEvents]);
               }

            }
         }
         //if (cut.at("MUGAST").at(Form("CUT_ETOF_MG%d", (*Mugast).TelescopeNumber[MugastEvents]))->IsInside((*Mugast).DSSD_E[MugastEvents], AlignT((*Mugast).TelescopeNumber[MugastEvents], (*Mugast).DSSD_Y[MugastEvents], (*Mugast).DSSD_T[MugastEvents])))
         for (const auto & particle : particles)
         {
            if (particle =="ANY") continue;
            if (cut.at("MUGAST").at(Form("E_TOF_MG%d_%s", (*Mugast).TelescopeNumber[MugastEvents], particle.c_str()))->IsInside((*Mugast).DSSD_E[MugastEvents], AlignT((*Mugast).TelescopeNumber[MugastEvents], (*Mugast).DSSD_Y[MugastEvents], (*Mugast).DSSD_T[MugastEvents])))
            {
               if (AGATA_GOOD)
               {
                  //Loop over gammas
                  for (int ii = 0; ii < AddE.GetSize(); ii++)
                  {
                     if (AddE[ii] > 10)
                     {
                        Fill(pData.AGATA.mEx_DC[IdentifiedNucleus.mass][IdentifiedNucleus.name][particle], Ex[MugastEvents], CorrectDoppler(p4, AddE[ii] / 1E3, AddX[ii], AddY[ii], AddZ[ii]));
                     }
                     if (AddE[ii] > 320 && AddE[ii] < 390)
                     {
                        if(particle=="d") cout<< "Filling gammas "<<" id mass: "<<  IdentifiedNucleus.mass<<" idname: "<< IdentifiedNucleus.name<<" particle: "<< particle<<" gamma: "<< "360 keV"<<" ELab: "<<  ThetaLab[MugastEvents]<<" ELab: "<< ELab[MugastEvents]<<std::endl;
                        Fill(pData.SI.mELab_ThetaLab[IdentifiedNucleus.mass][IdentifiedNucleus.name][particle]["360 keV"], ThetaLab[MugastEvents], ELab[MugastEvents]);
                     }
                  }
               }
               Fill(pData.SI.hEx[IdentifiedNucleus.mass][IdentifiedNucleus.name][particle], Ex[MugastEvents]);
               Fill(pData.SI.mEx_TW[IdentifiedNucleus.mass][IdentifiedNucleus.name][particle], *TW, Ex[MugastEvents]);
               Fill(pData.SI.mELab_ThetaLab[IdentifiedNucleus.mass][IdentifiedNucleus.name][particle]["ANY"], ThetaLab[MugastEvents], ELab[MugastEvents]);
               Fill(pData.SI.mELab_ThetaLab["ANY"]["ANY"][particle]["ANY"], ThetaLab[MugastEvents], ELab[MugastEvents]);
               Fill(pData.SI.mECM_ThetaCM[IdentifiedNucleus.mass][IdentifiedNucleus.name][particle], ThetaCM[MugastEvents], Ecm[MugastEvents]);
            }
         }
         MugastEvents++;
      }
   }
   

   catch (std::out_of_range &e)
   {
      std::cerr << "Silicon Physics loop :" <<e.what() << std::endl;
   }
   //MUST2
   for (int ii = 0; ii < (*MUST2).Si_E.size(); ii++)
   {
      Fill(pConf.SI.mE_TOF[Form("MM%i", (*MUST2).TelescopeNumber[ii])], (*MUST2).Si_E[ii], (*MUST2).Si_T[ii]);
      if (IdentifiedNucleus.found)
      {
         Fill(pData.SI.mE_TOF[Form("MM%i", (*MUST2).TelescopeNumber[ii])][IdentifiedNucleus.mass][IdentifiedNucleus.name], (*MUST2).Si_E[ii], (*MUST2).Si_T[ii]);
      }
   }

   //AGATA
   // Fill(pConf.AGATA.hAddTS_LTS,*AddTS-*LTS);
   for (int ii = 0; ii < *nbAdd; ii++)
   {
      //    Fill(pConf.AGATA.mmAGATA3D,AddX[ii], AddY[ii], AddZ[ii]);
   }
   return kTRUE;
}

void Selector::SlaveTerminate()
{
   // The SlaveTerminate() function is called after all entries or objects
   // have been processed. When running with PROOF SlaveTerminate() is called
   // on each slave server.
   TFile *top = new TFile("my_output.root", "recreate");
   TIter iter(fOutput);
   TObject *obj;
   TCanvas *canvas = new TCanvas("canvas", "canvas");
   while ((obj = iter()))
   {
      obj->Write();
      canvas = new TCanvas("canvas","canvas");
      obj->Draw();
      std::cout << "Plot Canvas\n";
      if (debug){
         while (gROOT->FindObject("canvas") != NULL)
         {
            gSystem->Sleep(200);
            gClient->HandleInput();
            gSystem->ProcessEvents();
         }
      }
   }
   top->Close();
   std::cout << "Total 46 Ar identified :"<< count46Ar <<std::endl;
}

void Selector::Terminate()
{
   // The Terminate() function is the last function to be called during
   // a query. It always runs on the client, it can be used to present
   // the results graphically or save the results to file.

}

inline Double_t Selector::CorrectDoppler(const TLorentzVector & p4, const Double_t & Egamma, const Double_t & X, const Double_t & Y, const Double_t & Z){
   TLorentzVector pgamma(Egamma, 0, 0, Egamma);
   TVector3 PosGamma(X, Y, Z+agata_Zshift);
   pgamma.SetPhi(PosGamma.Phi());
   pgamma.SetTheta(PosGamma.Theta());
   pgamma.SetE(Egamma);
   pgamma.Boost(-p4.BoostVector());
   return pgamma.Energy();
}

inline Double_t Selector::CorrectTOF(const Double_t & tof, const TVector3 & pos, const Double_t & Ek, const std::string & coeff){
   //return ((tof)-2.99792* pos.Mag() *std::stof(coeff)*(Ek+mass["1H"])/(sqrt((Ek+mass["1H"])*(Ek+mass["1H"])-mass["1H"]*mass["1H"])));
   return (tof-std::stof(coeff))/ pos.Mag() ;
   //return tof - 3.*pos.Mag()*sqrt(mass["1H"]/(2*Ek));
}

bool Selector::GetSettings(){

   std::ifstream ifile("GraphMask.txt");
   if(ifile){  
      std::ifstream file("GraphMask.txt");
      std::string Line;
      while(std::getline(file, Line)){
            std::istringstream str(Line);
            std::string Graph;
            str >> Graph;
            std::string enabled;
            str >> enabled;
            if (!enabled.compare("false")){
               TObject* to_delete = fOutput->FindObject(Graph.c_str());
               if (to_delete) {
                  fOutput->Remove(to_delete);
                  std::cout<<"Deleted :"<< Graph <<std::endl;
                  }
            }
      }
   }else{
      std::ofstream file("GraphMask.txt");
      TIter iter(fOutput);
      TObject *obj;
   while ((obj = iter()))
   {
      file << obj->GetName() << "\t\t\t"<< "false" <<std::endl;
   }
   file.close();
   }
   return true;
}

inline bool Selector::Fill(TH1* histo, const Double_t &data1)
{
   if (!fOutput->FindObject(histo)){
      return false;
   }
   else
   {
      histo->Fill(data1);
    //  std::cout << "filling :"<<data1 <<std::endl;
   }
   return true;
}

inline bool Selector::Fill(TH2* histo, const Double_t & data1, const Double_t & data2)
{
   if (!fOutput->FindObject(histo)){
      return false;
}else
   {
      histo->Fill(data1, data2);
   }
   return true;
}

inline bool Selector::Fill(TH3* histo, const Double_t & data1, const Double_t & data2, const Double_t & data3)
{
   if (!fOutput->FindObject(histo)){
      return false;
}else
   {
      histo->Fill(data1, data2, data3);
   }
   return true;
}