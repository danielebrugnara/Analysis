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
#include <TH2.h>
#include <TStyle.h>

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

    
    cal_spec    = new TH1D("cal_spec", "cal_spec", 3000, 0, 3000);
    fOutput->Add(cal_spec);
    core_spec   = new TH1D("core_spec", "core_spec", 3000, 0, 3000);
    fOutput->Add(core_spec);
    core_spec_DC   = new TH1D("core_spec_DC", "core_spec_DC", 3000, 0, 3000);
    fOutput->Add(core_spec_DC);
    core_spec_DC_pos_1  = new TH1D("core_spec_DC_pos_1", "core_spec_DC_pos_1", 3000, 0, 3000);
    fOutput->Add(core_spec_DC_pos_1);
    core_spec_DC_pos_2  = new TH1D("core_spec_DC_pos_2", "core_spec_DC_pos_2", 3000, 0, 3000);
    fOutput->Add(core_spec_DC_pos_2);
    core_gg   = new TH2D("mgamma_gamma", "core_gg", 3000, 0, 3000, 3000, 0, 3000);
    fOutput->Add(core_gg);
    core_gg_DC   = new TH2D("core_gg_DC", "core_gg_DC", 3000, 0, 3000, 3000, 0, 3000);
    fOutput->Add(core_gg_DC);
    core_gg_DC_pos_1  = new TH2D("core_gg_DC_pos_1", "core_gg_DC_pos_1", 3000, 0, 3000, 3000, 0, 3000);
    fOutput->Add(core_gg_DC_pos_1);
    core_gg_DC_pos_2  = new TH2D("core_gg_DC_pos_2", "core_gg_DC_pos_2", 3000, 0, 3000, 3000, 0, 3000);
    fOutput->Add(core_gg_DC_pos_2);
//    addb_spec   = new TH1D("addb_spec", "addb_spec", 3000, 0, 3000);
//    fOutput->Add(addb_spec);
//    addb_spec_DC        = new TH1D("addb_spec", "addb_spec", 3000, 0, 3000);
//    fOutput->Add(addb_spec_DC);
//    addb_spec_DC_pos_1  = new TH1D("addb_spec_pos_1", "addb_spec_pos_1", 3000, 0, 3000);
//    fOutput->Add(addb_spec_DC_pos_1);
//    addb_spec_DC_pos_2  = new TH1D("addb_spec_pos_2", "addb_spec_pos_2", 3000, 0, 3000);
//    fOutput->Add(addb_spec_DC_pos_2);
    dist        = new TH1D("dist_spec", "dist_spec", 3000, 0, 300);
    fOutput->Add(dist);

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

    std::unordered_map<int,Hit> hits_analyzed;
    std::unordered_map<int,Hit>::iterator it;
    double cr;
    for (const auto & h: Hits){
        if ((it = hits_analyzed.find(cr = h.GetCrystal()))== hits_analyzed.end()){
            hits_analyzed.emplace(cr, h);
        }else{
            it->second = Hit(cr,
                             h.GetEnergy()+it->second.GetEnergy(),
                             (h.GetX()*h.GetEnergy()+it->second.GetX()*it->second.GetEnergy())/(h.GetEnergy()+it->second.GetEnergy()),
                             (h.GetY()*h.GetEnergy()+it->second.GetY()*it->second.GetEnergy())/(h.GetEnergy()+it->second.GetEnergy()),
                             (h.GetZ()*h.GetEnergy()+it->second.GetZ()*it->second.GetEnergy())/(h.GetEnergy()+it->second.GetEnergy()),
                             h.GetDetector(),
                             h.GetEnergy()>it->second.GetEnergy() ? h.GetSegment() : it->second.GetSegment());
        }
    }

    
    double tot_en = 0;
    for (const auto & hh: hits_analyzed){
        tot_en += hh.second.GetEnergy();
        TVector3 vec(hh.second.GetX(),hh.second.GetY(),hh.second.GetZ());

        core_spec->Fill(hh.second.GetEnergy());
        core_spec_DC->Fill( ComputeDoppler(vec, hh.second.GetEnergy()));
        core_spec_DC_pos_1->Fill( ComputeDoppler(vec,em_position_1, hh.second.GetEnergy()));
        core_spec_DC_pos_2->Fill( ComputeDoppler(vec,em_position_2, hh.second.GetEnergy()));

        for (const auto & hh2: hits_analyzed){
            if (hh.first == hh2.first) continue;
            TVector3 vec2(hh2.second.GetX(),hh2.second.GetY(),hh2.second.GetZ());
            double distance = (vec-vec2).Mag();
            core_gg->Fill(hh.second.GetEnergy(), hh2.second.GetEnergy());
            core_gg_DC->Fill(ComputeDoppler(vec,hh.second.GetEnergy()), 
                                ComputeDoppler(vec2, hh2.second.GetEnergy()));
            core_gg_DC_pos_1->Fill(  ComputeDoppler(vec,em_position_1,hh.second.GetEnergy()), 
                                        ComputeDoppler(vec2,em_position_1,hh2.second.GetEnergy()));
            core_gg_DC_pos_2->Fill(  ComputeDoppler(vec,em_position_2,hh.second.GetEnergy()), 
                                        ComputeDoppler(vec2,em_position_2,hh2.second.GetEnergy()));
        }
    }

    cal_spec->Fill(tot_en);
    return kTRUE;
}

void Selector::SlaveTerminate()
{
    // The SlaveTerminate() function is called after all entries or objects
    // have been processed. When running with PROOF SlaveTerminate() is called
    // on each slave server.

}

void Selector::Terminate()
{
    // The Terminate() function is the last function to be called during
    // a query. It always runs on the client, it can be used to present
    // the results graphically or save the results to file.

}

double Selector::ComputeDoppler(const TVector3 & vec, const double & en){
    return en/(sqrt(1-beta*beta)*(1. + beta*vec.CosTheta()));
}

double Selector::ComputeDoppler(const TVector3 & vec,const TVector3 & pos, const double & en){
    return en/(sqrt(1-beta*beta)*(1. + beta*(vec-pos).CosTheta()));
}
