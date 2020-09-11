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


    cal_spec    = new TH1D("cal_spec", "cal_spec", 3000, 0, 3000);
    fOutput->Add(cal_spec);

    //Core graphs
    core_spec   = new TH1D("core_spec", "core_spec", 3000, 0, 3000);
    fOutput->Add(core_spec);
    core_spec_DC   = new TH1D("core_spec_DC", "core_spec_DC", 3000, 0, 3000);
    fOutput->Add(core_spec_DC);
    core_spec_DC_pos_1  = new TH1D("core_spec_DC_pos_1", "core_spec_DC_pos_1", 3000, 0, 3000);
    fOutput->Add(core_spec_DC_pos_1);
    core_spec_DC_pos_2  = new TH1D("core_spec_DC_pos_2", "core_spec_DC_pos_2", 3000, 0, 3000);
    fOutput->Add(core_spec_DC_pos_2);
    core_gg   = new TH2D("core_gg", "core_gg", 3000, 0, 3000, 3000, 0, 3000);
    fOutput->Add(core_gg);
    core_gg_DC   = new TH2D("core_gg_DC", "core_gg_DC", 3000, 0, 3000, 3000, 0, 3000);
    fOutput->Add(core_gg_DC);
    core_gg_DC_pos_1  = new TH2D("core_gg_DC_pos_1", "core_gg_DC_pos_1", 3000, 0, 3000, 3000, 0, 3000);
    fOutput->Add(core_gg_DC_pos_1);
    core_gg_DC_pos_2  = new TH2D("core_gg_DC_pos_2", "core_gg_DC_pos_2", 3000, 0, 3000, 3000, 0, 3000);
    fOutput->Add(core_gg_DC_pos_2);

    //Addback graphs
    addb_spec   = new TH1D("addb_spec", "addb_spec", 3000, 0, 3000);
    fOutput->Add(addb_spec);
    addb_spec_DC        = new TH1D("addb_spec_DC", "addb_spec_DC", 3000, 0, 3000);
    fOutput->Add(addb_spec_DC);
    addb_spec_DC_pos_1  = new TH1D("addb_spec_pos_1", "addb_spec_pos_1", 3000, 0, 3000);
    fOutput->Add(addb_spec_DC_pos_1);
    addb_spec_DC_pos_2  = new TH1D("addb_spec_pos_2", "addb_spec_pos_2", 3000, 0, 3000);
    fOutput->Add(addb_spec_DC_pos_2);
    addb_gg   = new TH2D("addb_gg", "addb_gg", 3000, 0, 3000, 3000, 0, 3000);
    fOutput->Add(addb_gg);
    addb_gg_DC   = new TH2D("addb_gg_DC", "addb_gg_DC", 3000, 0, 3000, 3000, 0, 3000);
    fOutput->Add(addb_gg_DC);
    addb_gg_DC_pos_1  = new TH2D("addb_gg_DC_pos_1", "addb_gg_DC_pos_1", 3000, 0, 3000, 3000, 0, 3000);
    fOutput->Add(addb_gg_DC_pos_1);
    addb_gg_DC_pos_2  = new TH2D("addb_gg_DC_pos_2", "addb_gg_DC_pos_2", 3000, 0, 3000, 3000, 0, 3000);
    fOutput->Add(addb_gg_DC_pos_2);

    //Other Graphs
    dist        = new TH1D("dist_spec", "dist_spec", 3000, 0, 300);
    fOutput->Add(dist);
    coreID_coreID   = new TH2D("coreID_coreID", "coreID_coreID", 100, 0, 100, 100, 0, 100);
    fOutput->Add(coreID_coreID);
    dist_coreID        = new TH2D("dist_cireID", "dist_coreID", 3000, 0, 300, 60, 0, 60);
    fOutput->Add(dist_coreID);
    for (int ii=0; ii<60; ++ii){
        core_dist[ii] = new TH2D(Form("core_%i", ii), Form("core_%i", ii), 60, 0, 60, 3000, 0, 300);
        //fOutput->Add(core_dist.at(ii));
    }

    std::unordered_map<int, std::vector<int>> nearby_det;
    nearby_det.reserve(35);

    nearby_det[0] = {1, 2, 3, 12, 13};
    nearby_det[1] = {0, 2, 3, 5, 20, 22};
    nearby_det[2] = {0, 1, 13, 16, 17, 20};
    nearby_det[3] = {0, 1, 4, 5, 6};
    nearby_det[4] = {3, 5, 6, 8};
    nearby_det[5] = {1, 3, 4, 22, 23};
    nearby_det[6] = {3, 4, 7, 8, 9};
    nearby_det[7] = {6, 8, 9, 11, 26};
    nearby_det[8] = {4, 6, 7, 26};
    nearby_det[9] = {6, 7, 10, 11, 12};
    nearby_det[10] = {9, 11, 12, 14, 29, 31};
    nearby_det[11] = {7, 9, 10, 29};
    nearby_det[12] = {0, 9, 10, 13, 14};
    nearby_det[13] = {0, 2, 12, 14, 16, 35};
    nearby_det[14] = {10, 12, 13, 31, 32, 35};
    nearby_det[15] = {16, 17, 33};
    nearby_det[16] = {2, 13, 15, 17, 33, 35};
    nearby_det[17] = {2, 15, 16, 19, 20};
    nearby_det[18] = {19, 20, 21, 22};
    nearby_det[19] = {17, 18, 20};
    nearby_det[20] = {1, 2, 17, 18, 19, 22};
    nearby_det[21] = {18, 22, 23};
    nearby_det[22] = {1, 5, 18, 20, 21, 23};
    nearby_det[23] = {5, 21, 22};
    nearby_det[24] = {25, 26};
    nearby_det[25] = {24, 26};
    nearby_det[26] = {7, 8, 24, 25};
    nearby_det[27] = {28, 29, 30, 31};
    nearby_det[28] = {27, 29};
    nearby_det[29] = {10, 11, 27, 28, 31};
    nearby_det[30] = {27, 31, 32};
    nearby_det[31] = {10, 14, 27, 29, 30, 32};
    nearby_det[32] = {14, 30, 31, 34, 35};
    nearby_det[33] = {15, 16, 34, 35};
    nearby_det[34] = {32, 33, 35};
    nearby_det[35] = {13, 14, 16, 32, 33, 34};

    std::vector<graphEdge<bool>> edges;
    std::vector<Crystal> payloads;

    std::unordered_map<int, std::string> dictionary_crystal_name;
    dictionary_crystal_name[0] = "A";
    dictionary_crystal_name[1] = "B";
    dictionary_crystal_name[2] = "C";

    payloads.resize(nearby_det.size());
    for (const auto& it_start: nearby_det){
        for (const auto& it_end: it_start.second){
            edges.push_back({it_start.first,it_end,true});
        }
        std::stringstream str_tmp;
        str_tmp << std::setw(2) << std::setfill('0') << it_start.first/3;
        str_tmp << dictionary_crystal_name[it_start.first%3];
        payloads[it_start.first] = {it_start.first,str_tmp.str()};
    }

    addback_graph = new DiaGraph<bool, Crystal>(edges,payloads);
    //std::cout << *addback_graph;
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

    //Merge hits in the same crystal
    std::unordered_map<int,Hit> hits_analyzed;
    std::unordered_map<int,Hit>::iterator hits_analyzed_it;
    int cr;
    for (const auto & h: Hits){
        if ((hits_analyzed_it = hits_analyzed.find(cr = h.GetCrystal()))== hits_analyzed.end()){
            hits_analyzed.emplace(cr, h);
        }else{
            hits_analyzed_it->second = Hit(cr,
                                           h.GetEnergy()+hits_analyzed_it->second.GetEnergy(),
                                           (h.GetX()*h.GetEnergy()+hits_analyzed_it->second.GetX()*hits_analyzed_it->second.GetEnergy())/(h.GetEnergy()+hits_analyzed_it->second.GetEnergy()),
                                           (h.GetY()*h.GetEnergy()+hits_analyzed_it->second.GetY()*hits_analyzed_it->second.GetEnergy())/(h.GetEnergy()+hits_analyzed_it->second.GetEnergy()),
                                           (h.GetZ()*h.GetEnergy()+hits_analyzed_it->second.GetZ()*hits_analyzed_it->second.GetEnergy())/(h.GetEnergy()+hits_analyzed_it->second.GetEnergy()),
                                           h.GetEnergy()>hits_analyzed_it->second.GetEnergy() ? h.GetSegment() : hits_analyzed_it->second.GetSegment(),
                                           h.GetDetector());
        }
    }

    //Compute addback
    std::unordered_map<int, Hit> addback = hits_analyzed;
    std::unordered_map<int, Hit>::iterator addback_it1;
    std::unordered_map<int, Hit>::iterator addback_it2;
    for (addback_it1 = addback.begin(); addback_it1 != addback.end(); ++addback_it1){
        if (addback.find(addback_it1->first) == addback.end()) continue;
        auto near_crystals = addback_graph->getAdjacent(addback_it1->first);
        for (addback_it2 = addback.begin(); addback_it2 != addback.end(); ++addback_it2){
            if (addback.find(addback_it1->first) == addback.end()) break;
            bool present = false;
            for (const auto& it: near_crystals){
                if (it.second->ID == addback_it2->first)
                    present = true;
            }
            if(present) {
                int new_key = addback_it1->second.GetEnergy() > addback_it2->second.GetEnergy() ? addback_it1->first
                                                                                                : addback_it2->first;
                int delete_key = addback_it1->second.GetEnergy() < addback_it2->second.GetEnergy() ? addback_it1->first
                                                                                                   : addback_it2->first;

                addback[new_key] = Hit(new_key,
                                       addback[new_key].GetEnergy() + addback[delete_key].GetEnergy(),
                                       addback[new_key].GetX(),
                                       addback[new_key].GetY(),
                                       addback[new_key].GetZ(),
                                       addback[new_key].GetSegment(),
                                       addback[new_key].GetDetector());
                addback.erase(delete_key);
            }
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
            dist->Fill(distance);
            dist_coreID->Fill(distance, (double)hh.first);
            coreID_coreID->Fill((double)hh.first,(double)hh2.first);
            try {
                core_dist.at(hh.first)->Fill(hh2.first, distance);
            }catch(...){};

            core_gg->Fill(hh.second.GetEnergy(), hh2.second.GetEnergy());
            core_gg_DC->Fill(ComputeDoppler(vec,hh.second.GetEnergy()),
                             ComputeDoppler(vec2, hh2.second.GetEnergy()));
            core_gg_DC_pos_1->Fill(  ComputeDoppler(vec,em_position_1,hh.second.GetEnergy()),
                                     ComputeDoppler(vec2,em_position_1,hh2.second.GetEnergy()));
            core_gg_DC_pos_2->Fill(  ComputeDoppler(vec,em_position_2,hh.second.GetEnergy()),
                                     ComputeDoppler(vec2,em_position_2,hh2.second.GetEnergy()));
        }
    }

    for (const auto& hh: addback){
        tot_en += hh.second.GetEnergy();
        TVector3 vec(hh.second.GetX(),hh.second.GetY(),hh.second.GetZ());

        addb_spec->Fill(hh.second.GetEnergy());
        addb_spec_DC->Fill( ComputeDoppler(vec, hh.second.GetEnergy()));
        addb_spec_DC_pos_1->Fill( ComputeDoppler(vec,em_position_1, hh.second.GetEnergy()));
        addb_spec_DC_pos_2->Fill( ComputeDoppler(vec,em_position_2, hh.second.GetEnergy()));

        for (const auto & hh2: addback){
            if (hh.first == hh2.first) continue;
            TVector3 vec2(hh2.second.GetX(),hh2.second.GetY(),hh2.second.GetZ());
            double distance = (vec-vec2).Mag();
            dist->Fill(distance);

            addb_gg->Fill(hh.second.GetEnergy(), hh2.second.GetEnergy());
            addb_gg_DC->Fill(ComputeDoppler(vec,hh.second.GetEnergy()),
                             ComputeDoppler(vec2, hh2.second.GetEnergy()));
            addb_gg_DC_pos_1->Fill(  ComputeDoppler(vec,em_position_1,hh.second.GetEnergy()),
                                     ComputeDoppler(vec2,em_position_1,hh2.second.GetEnergy()));
            addb_gg_DC_pos_2->Fill(  ComputeDoppler(vec,em_position_2,hh.second.GetEnergy()),
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
