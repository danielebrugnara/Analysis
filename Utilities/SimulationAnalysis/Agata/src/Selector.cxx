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
    geom   = new TH3D("geom", "geom", 100, -350, 350, 100, -350, 350, 100, -350, 350);
    fOutput->Add(geom);
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
    target_spec   = new TH1D("target_spec", "target_spec", 3000, 0, 3000);
    fOutput->Add(target_spec);


    for (int i=0; i<42; ++i){
        active_crystals[i] = true;
    }

    active_crystals[4] = false;    
    active_crystals[6] = false;    
    active_crystals[7] = false;    
    active_crystals[13] = false;    
    active_crystals[20] = false;    
    active_crystals[20] = false;    
    active_crystals[26] = false;    

    std::unordered_map<int, std::vector<int>> nearby_det;
    nearby_det.reserve(35);


    nearby_det[0] = {1, 2, 3, 12, 13};
    nearby_det[1] = {0, 2, 3, 5, 20, 22};
    nearby_det[2] = {0, 1, 13, 16, 17, 20};
    nearby_det[3] = {0, 1, 4, 5, 6};
    nearby_det[4] = {3, 5, 6, 8, 24};
    nearby_det[5] = {1, 3, 4, 22, 23};
    nearby_det[6] = {3, 4, 7, 8, 9};
    nearby_det[7] = {6, 8, 9, 11, 29, 31};
    nearby_det[8] = {4, 6, 7, 24, 25, 29};
    nearby_det[9] = {6, 7, 10, 11, 12};
    nearby_det[10] = {9, 11, 12, 14, 35, 37};
    nearby_det[11] = {7, 9, 10, 31, 32, 35};
    nearby_det[12] = {0, 9, 10, 13, 14};
    nearby_det[13] = {0, 2, 12, 14, 16, 41};
    nearby_det[14] = {10, 12, 13, 37, 38, 41};
    nearby_det[15] = {16, 17, 39};
    nearby_det[16] = {2, 13, 15, 17, 39, 41};
    nearby_det[17] = {2, 15, 16, 19, 20};
    nearby_det[18] = {19, 20, 21, 22};
    nearby_det[19] = {17, 18, 20};
    nearby_det[20] = {1, 2, 17, 18, 19, 22};
    nearby_det[21] = {18, 22, 23};
    nearby_det[22] = {1, 5, 18, 20, 21, 23};
    nearby_det[23] = {5, 21, 22};
    nearby_det[24] = {4, 8, 25};
    nearby_det[25] = {8, 24, 28, 29};
    nearby_det[26] = {};
    nearby_det[27] = {28, 29, 30, 31};
    nearby_det[28] = {25, 27, 29};
    nearby_det[29] = {7, 8, 25, 27, 28, 31};
    nearby_det[30] = {27, 31, 32};
    nearby_det[31] = {7, 11, 27, 29, 30, 32};
    nearby_det[32] = {11, 30, 31, 34, 35};
    nearby_det[33] = {34, 35, 36, 37};
    nearby_det[34] = {32, 33, 35};
    nearby_det[35] = {10, 11, 32, 33, 34, 37};
    nearby_det[36] = {33, 37, 38};
    nearby_det[37] = {10, 14, 33, 35, 36, 38};
    nearby_det[38] = {14, 36, 37, 40, 41};
    nearby_det[39] = {15, 16, 40, 41};
    nearby_det[40] = {38, 39, 41};
    nearby_det[41] = {13, 14, 16, 38, 39, 40};

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

    energy_threashold = 30.;
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

    //Merge hits in the same crystal and segment
    std::unordered_map<int,std::unordered_map<int,Hit>> hits_analyzed; //key: crystal, segment
    std::unique_ptr<Hit> target_hit = nullptr;
    for (const auto & h: Hits){
        if (h.GetDetector() == 0){//Agata hit
            int crystal = h.GetCrystal();
            if (!active_crystals[crystal]) continue;
            if (hits_analyzed.find(crystal) == hits_analyzed.end()){
                //First hit in crystal
                std::unordered_map<int, Hit> tmp;
                tmp.emplace(h.GetSegment(),h);
                hits_analyzed.emplace(crystal, tmp);
            }else{
                //Crystal already hit
                int segment = h.GetSegment();
                if (hits_analyzed.at(crystal).find(segment) == hits_analyzed.at(crystal).end()){
                    //First hit in segment
                    hits_analyzed.at(crystal).emplace(segment, h);
                }else{
                    //Segment already hit
                    double esum = h.GetEnergy()+hits_analyzed.at(crystal).at(segment).GetEnergy();
                    Hit tmp_hit(    crystal,
                                    esum,
                                    (h.GetX()*h.GetEnergy()+hits_analyzed.at(crystal).at(segment).GetX()*hits_analyzed.at(crystal).at(segment).GetEnergy())/esum,
                                    (h.GetY()*h.GetEnergy()+hits_analyzed.at(crystal).at(segment).GetY()*hits_analyzed.at(crystal).at(segment).GetEnergy())/esum,
                                    (h.GetZ()*h.GetEnergy()+hits_analyzed.at(crystal).at(segment).GetZ()*hits_analyzed.at(crystal).at(segment).GetEnergy())/esum,
                                    segment,
                                    h.GetDetector());
                    hits_analyzed.at(crystal).at(segment) = tmp_hit;
                }
            }
        }else if (h.GetDetector() == 3){//Target hit
            if (target_hit == nullptr )
                target_hit.reset(new Hit(h));
            else
                target_hit.reset(new Hit(
                        h.GetCrystal(),
                        h.GetEnergy()+target_hit->GetEnergy(),
                        (h.GetX()*h.GetEnergy()+target_hit->GetX()*target_hit->GetEnergy())/(h.GetEnergy()+target_hit->GetEnergy()),
                        (h.GetY()*h.GetEnergy()+target_hit->GetY()*target_hit->GetEnergy())/(h.GetEnergy()+target_hit->GetEnergy()),
                        (h.GetZ()*h.GetEnergy()+target_hit->GetZ()*target_hit->GetEnergy())/(h.GetEnergy()+target_hit->GetEnergy()),
                        target_hit->GetSegment(),
                        h.GetDetector()));
        }else{//Nothing should be here
            throw std::runtime_error("Hit detector not recognized\n");
        }
    }

    std::unordered_map<int, Hit> cores;

    for (const auto& it_cry: hits_analyzed){
        std::pair<int,double> max_en = {-1,0}; //Segment, energy
        double total_en = 0;
        for (const auto& it_segm: hits_analyzed.at(it_cry.first)){
            total_en += it_segm.second.GetEnergy();
            if (it_segm.second.GetEnergy() > max_en.second)
                max_en = std::make_pair(it_segm.first, it_segm.second.GetEnergy());
        }
        if (max_en.first == -1) continue;
        cores.emplace(  it_cry.first,
                        Hit( it_cry.first,
                                total_en,
                                it_cry.second.at(max_en.first).GetX(),
                                it_cry.second.at(max_en.first).GetY(),
                                it_cry.second.at(max_en.first).GetZ(),
                                max_en.first,
                                it_cry.second.at(max_en.first).GetDetector()));
    }

    //Deleting under-threashold stuff
    for (auto it=cores.cbegin(); it !=cores.cend();){
        if (it->second.GetEnergy()<energy_threashold){
            cores.erase(it++);
        }else{
            ++it;
        }
    }

    //Compute addback
    std::unordered_map<int, Hit> addback;
    std::unordered_map<int, bool> visited;
    std::deque<int> visit_que;

    for (const auto& it_core: cores){
        if (visited.find(it_core.first) != visited.end()) continue;
        visit_que.push_back(it_core.first);
        std::pair<int, double> max_element = {-1, 0};
        double total_energy = 0;
        while(!visit_que.empty()){
            auto element = visit_que.front();
            visit_que.pop_front();
            visited.emplace(element,true);

            double element_energy = cores.at(element).GetEnergy();
            total_energy += element_energy;
            if (max_element.second < element_energy)
                max_element = std::make_pair(element, element_energy);

            auto near_crystals = addback_graph->getAdjacent(it_core.first);
            for (const auto& it_near: near_crystals){
                if (visited.find(it_near.second->ID) != visited.end()) continue;
                if (cores.find(it_near.second->ID) == cores.end()) continue;
                visit_que.push_back( it_near.second->ID);
            }
        }

        const Hit& maxhit = cores.at(max_element.first);
        addback.emplace(max_element.first, Hit( max_element.first,
                                                total_energy,
                                                maxhit.GetX(),
                                                maxhit.GetY(),
                                                maxhit.GetZ(),
                                                maxhit.GetSegment(),
                                                maxhit.GetDetector()));
    }


    //Filling histograms
    if (target_hit != nullptr)
        target_spec->Fill(target_hit->GetEnergy());

    double tot_en = 0;
    for (const auto & hh: cores){
        tot_en += hh.second.GetEnergy();
        TVector3 vec(hh.second.GetX(),hh.second.GetY(),hh.second.GetZ());

        geom->Fill(hh.second.GetX(), hh.second.GetY(), hh.second.GetZ());
        core_spec->Fill(hh.second.GetEnergy());
        core_spec_DC->Fill( ComputeDoppler(vec, hh.second.GetEnergy()));
        core_spec_DC_pos_1->Fill( ComputeDoppler(vec,em_position_1, hh.second.GetEnergy()));
        core_spec_DC_pos_2->Fill( ComputeDoppler(vec,em_position_2, hh.second.GetEnergy()));

        for (const auto & hh2: cores){
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
