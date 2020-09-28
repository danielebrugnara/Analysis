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


    cal_spec   = new TH1D("cal_spec", "cal_spec", 3000, 0, 3000);
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

    target_spec   = new TH1D("target_spec", "target_spec", 3000, 0, 3000);
    fOutput->Add(target_spec);

    //Individual crystal core spectra
    for(int i=0; i<=45; ++i){
        crystal_spectra.push_back(new TH1D(Form("crtstal_spectra_%i", i),
                                           Form("crtstal_spectra_%i", i),
                                           4000, 0, 4000));
        fOutput->Add(crystal_spectra.back());;
    }


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
    //std::cout << *addback_graph;

    simu_to_data[0] = 0;
    simu_to_data[1] = 1;
    simu_to_data[2] = 2;
    simu_to_data[3] = 3;
    simu_to_data[4] = 4;
    simu_to_data[5] = 5;
    simu_to_data[6] = 6;
    simu_to_data[7] = 7;
    simu_to_data[8] = 8;
    simu_to_data[9] = 9;
    simu_to_data[10] = 10;
    simu_to_data[11] = 11;
    simu_to_data[12] = 12;
    simu_to_data[13] = 13;
    simu_to_data[14] = 14;
    simu_to_data[15] = 15;
    simu_to_data[16] = 16;
    simu_to_data[17] = 17;
    simu_to_data[18] = 18;
    simu_to_data[19] = 19;
    simu_to_data[20] = 20;
    simu_to_data[21] = 21;
    simu_to_data[22] = 22;
    simu_to_data[23] = 23;
    simu_to_data[24] = 28;
    simu_to_data[25] = 29;
    //simu_to_data[26] = 26;
    simu_to_data[27] = 30;
    simu_to_data[28] = 31;
    simu_to_data[29] = 32;
    simu_to_data[30] = 33;
    simu_to_data[31] = 34;
    simu_to_data[32] = 35;
    simu_to_data[33] = 36;
    simu_to_data[34] = 37;
    simu_to_data[35] = 38;
    simu_to_data[36] = 39;
    simu_to_data[37] = 40;
    simu_to_data[38] = 41;
    simu_to_data[39] = 42;
    simu_to_data[40] = 43;
    simu_to_data[41] = 44;

    for (const auto& it: simu_to_data){
        //Making bijection of map
        data_to_simu[it.second] = it.first;
    }

    auto* threashold_file = new TFile("sampling_histos.root", "read");
    if (threashold_file->IsOpen())
        use_threasholds = true;
    else
        std::cerr << "Threasholds file not found!!\n";

    if (use_threasholds){
        for(const auto& it: data_to_simu){
            threasholds.emplace(it.second, (TH1*)threashold_file->Get(Form("probability_%i", it.second)));
        }
    }
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
    //std::unordered_map<int,std::unordered_map<int,Hit>> hits_analyzed; //key: crystal, segment
    std::map<std::pair<int, int>, Hit> hits_analyzed; //key: crystal, segment
    std::map<std::pair<int, int>, Hit>::iterator hits_analyzed_iterator;

    std::unique_ptr<Hit> target_hit = nullptr;
    double ecal = 0;
    for (const auto & h: Hits){
        if (h.GetEnergy() <= 0 ) continue;
        if (h.GetDetector() == 0){//Agata hit
            int crystal = h.GetCrystal();
            int segment = h.GetSegment();
            std::pair<int, int > key(crystal, segment);
            if (!active_crystals[crystal]) continue;
            ecal += h.GetEnergy();
            if ((hits_analyzed_iterator = hits_analyzed.find(key)) == hits_analyzed.end()){
                hits_analyzed.emplace(key, h);
            }else{
                double ensum = hits_analyzed_iterator->second.GetEnergy() + h.GetEnergy();
                auto& found = hits_analyzed_iterator->second;
                found = Hit(   crystal,
                               ensum,
                               (found.GetX()*found.GetEnergy()+h.GetX()*h.GetEnergy())/ensum,
                               (found.GetY()*found.GetEnergy()+h.GetY()*h.GetEnergy())/ensum,
                               (found.GetZ()*found.GetEnergy()+h.GetZ()*h.GetEnergy())/ensum,
                               segment,
                               h.GetDetector());
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

    std::unordered_map<int, int> max_keys_map; //key: crystal, content segment with highest energy
    std::unordered_map<int, double> etot_cores_map;
    for (const auto& it: hits_analyzed) {
        int crystal = it.first.first;
        int segment = it.first.second;
        if (etot_cores_map.find(crystal) == etot_cores_map.end()){
            max_keys_map.emplace(crystal, segment);
            etot_cores_map.emplace(crystal, it.second.GetEnergy());
        }else{
            etot_cores_map.at(crystal) += it.second.GetEnergy();
            if (it.second.GetEnergy() > hits_analyzed.find(std::make_pair(crystal, max_keys_map.at(crystal)))->second.GetEnergy()){
                max_keys_map.at(crystal) = segment;
            }
        }
    }

    std::unordered_map<int, Hit> cores;
    for (const auto& it: max_keys_map){
        auto max = hits_analyzed.find(std::make_pair(it.first,it.second));
        cores.emplace(  it.first,
                        Hit(    it.first,
                                etot_cores_map.at(it.first),
                                max->second.GetX(),
                                max->second.GetY(),
                                max->second.GetZ(),
                                max->second.GetSegment(),
                                max->second.GetDetector()
                        ));
    }

    //Deleting under-threashold stuff
    if (use_threasholds){
        for (auto it=cores.cbegin(); it !=cores.cend();){
            TH1* histo_ptr = threasholds.at(it->second.GetCrystal());
            if (it->second.GetEnergy()<150 && rand.Uniform()>histo_ptr->GetBinContent(histo_ptr->FindBin(it->second.GetEnergy()))){
                cores.erase(it++);
            }else{
                ++it;
            }
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
            //if (visited.find(element) != visited.end()) continue;
            visited.emplace(element,true);

            double element_energy = cores.at(element).GetEnergy();
            total_energy += element_energy;
            if (max_element.second < element_energy)
                max_element = std::make_pair(element, element_energy);

            auto near_crystals = addback_graph->getAdjacent(it_core.first);
            for (const auto& it_near: near_crystals){
                if (visited.find(it_near.second->ID) != visited.end()) continue;
                if (std::find(visit_que.begin(), visit_que.end(),it_near.second->ID) != visit_que.end()) continue;
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

    //Consistency check on cores and addback
    double core_tot_en = 0;
    for (const auto & it: cores) {
        core_tot_en += it.second.GetEnergy();
        if (it.second.GetEnergy() > ecal + energy_check_tolerance)
            throw std::runtime_error("Something wrong in cores energy\n");
    }
    if (!use_threasholds && abs(core_tot_en-ecal)>energy_check_tolerance)
        throw std::runtime_error("Something wrong in cores energy\n");

    double addback_tot_en = 0;
    for (const auto & it: addback) {
        addback_tot_en += it.second.GetEnergy();
        if (it.second.GetEnergy() > ecal + energy_check_tolerance)
            throw std::runtime_error("Something wrong in addback energy\n");
    }
    if (!use_threasholds && abs(addback_tot_en-ecal)>energy_check_tolerance)
        throw std::runtime_error("Something wrong in addback energy\n");

    //Filling histograms
    if (target_hit != nullptr)
        target_spec->Fill(target_hit->GetEnergy());

    for (const auto & hh: cores){
        TVector3 vec(hh.second.GetX(),hh.second.GetY(),hh.second.GetZ());

        geom->Fill(hh.second.GetX(), hh.second.GetY(), hh.second.GetZ());
        core_spec->Fill(hh.second.GetEnergy());
        core_spec_DC->Fill( ComputeDoppler(vec, hh.second.GetEnergy()));
        core_spec_DC_pos_1->Fill( ComputeDoppler(vec,em_position_1, hh.second.GetEnergy()));
        core_spec_DC_pos_2->Fill( ComputeDoppler(vec,em_position_2, hh.second.GetEnergy()));
        crystal_spectra[hh.second.GetCrystal()]->Fill(hh.second.GetEnergy());

        for (const auto & hh2: cores){
            if (hh.first == hh2.first) continue;
            TVector3 vec2(hh2.second.GetX(),hh2.second.GetY(),hh2.second.GetZ());
            //double distance = (vec-vec2).Mag();

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
        TVector3 vec(hh.second.GetX(),hh.second.GetY(),hh.second.GetZ());

        addb_spec->Fill(hh.second.GetEnergy());
        addb_spec_DC->Fill( ComputeDoppler(vec, hh.second.GetEnergy()));
        addb_spec_DC_pos_1->Fill( ComputeDoppler(vec,em_position_1, hh.second.GetEnergy()));
        addb_spec_DC_pos_2->Fill( ComputeDoppler(vec,em_position_2, hh.second.GetEnergy()));

        for (const auto & hh2: addback){
            if (hh.first == hh2.first) continue;
            TVector3 vec2(hh2.second.GetX(),hh2.second.GetY(),hh2.second.GetZ());
            //double distance = (vec-vec2).Mag();

            addb_gg->Fill(hh.second.GetEnergy(), hh2.second.GetEnergy());
            addb_gg_DC->Fill(ComputeDoppler(vec,hh.second.GetEnergy()),
                             ComputeDoppler(vec2, hh2.second.GetEnergy()));
            addb_gg_DC_pos_1->Fill(  ComputeDoppler(vec,em_position_1,hh.second.GetEnergy()),
                                     ComputeDoppler(vec2,em_position_1,hh2.second.GetEnergy()));
            addb_gg_DC_pos_2->Fill(  ComputeDoppler(vec,em_position_2,hh.second.GetEnergy()),
                                     ComputeDoppler(vec2,em_position_2,hh2.second.GetEnergy()));
        }

    }

    cal_spec->Fill(ecal);
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
    //auto* threashold_data_file = new TFile("threasholds.root", "read");
    auto* threashold_data_file = new TFile("threasholds.root", "read");
    if (threashold_data_file->IsOpen() && !use_threasholds)
        compute_threasholds = true;

    if (!compute_threasholds)
        return;

    std::unordered_map<int, TF1*> data_threasholds;
    for(const auto& it: data_to_simu){
        data_threasholds.emplace(it.second, (TF1*)threashold_data_file->Get(Form("normmodel_%i", it.first)));
    }

    //std::vector<TF1> simu_threasholds;
    bool show_canvas = true;
    std::unordered_map<int, TH1*> output_histos;
    for (unsigned int i=0; i<crystal_spectra.size(); ++i){
        if (data_threasholds.find(i) == data_threasholds.end()) continue;
        RooRealVar x("x", "x", 0, 150);
        x.setRange("fullrange", 0, 150);
        x.setRange("partial",1, 149);

        RooDataHist dh("dh", "dh", x, RooFit::Import(*(crystal_spectra[i])));

        //RooRealVar scale("scale", "scale", 1000/2., 1, 1000/2.*10);
//        RooRealVar a0("a0", "a0", 10, 0, 10);
//        RooRealVar a1("a1", "a1", 1., -100, 10.);
//        RooRealVar a2("a2", "a2", 1., -10, 10.);
        //RooChebychev model(Form("model_%i", i),
        //                    Form("model_%i", i),
        //                    x,
        //                    RooArgSet(a0, a1));
//        RooPolynomial pol(Form("pol_%i", i),
//                            Form("pol_%i", i),
//                            x,
//                            RooArgSet(a0, a1, a2));
        RooRealVar th_s("th_s", "th_s", 40, 0, 150);
        RooRealVar la_s("la_s", "la_s", 20, 1, 60);
        RooRealVar la2_s("la2_s", "la2_s", 0.05, 0.001, 0.2);
        RooGenericPdf pol(Form("pol_%i", i),
                            Form("pol_%i", i),
                            "(1+erf((th_s-x)/la_s))*(1-exp(-x*la2_s))",
                            RooArgSet(x,th_s, la_s, la2_s));

        RooRealVar threashold("threashold", "threashold", data_threasholds.at(i)->GetParameter("threashold"));
        threashold.setConstant();
        RooRealVar lambda("lambda", "lambda", data_threasholds.at(i)->GetParameter("lambda"));
        lambda.setConstant();

        RooGenericPdf erf(Form("erf_%i", i),
                            Form("erf_%i", i),
                            "(1+erf((x-threashold)/lambda))",
                            RooArgSet(x, threashold, lambda));

        RooRealVar sign_p("sign_p","sign_p", 0.5, 0.,1.);
        RooAddPdf model(Form("model_%i", i),
                        Form("model_%i", i),
                        RooArgList(pol, erf),
                            RooArgList(sign_p));
        //Plotting
        RooPlot* frame = x.frame(RooFit::Title(Form("Threashold_%i", i)));
        dh.plotOn(frame);
        auto *result = model.fitTo(dh, RooFit::Save());
        model.plotOn(frame);
        pol.plotOn(frame, RooFit::LineStyle(kDashed),RooFit::LineColor(kRed));
        erf.plotOn(frame, RooFit::LineStyle(kDashed),RooFit::LineColor(kGreen));

        RooArgList pars(*model.getParameters(RooArgSet(x)));
        RooArgSet prodset(model);
        RooProduct normPdf(Form("normmodel_%i", i),
                           Form("normmodel_%i", i),
                           prodset);

        if (show_canvas) {
            auto* cv = new TCanvas();
            frame->Draw();
            cv->WaitPrimitive();
            delete cv;
        }

        int Npts = 1000;
        data_threasholds.at(i)->SetRange(0, 150);
        data_threasholds.at(i)->SetNpx(Npts);
        auto* histo_data = data_threasholds.at(i)->GetHistogram();
        histo_data->SetName(Form("data_threasholds_%i", i));
        histo_data->Scale(1./(histo_data->GetBinContent(histo_data->GetMaximumBin())));

        auto* tmp_simu_thresh = normPdf.asTF(RooArgList(x), pars);
        tmp_simu_thresh->SetRange(0, 150);
        tmp_simu_thresh->SetNpx(Npts);
        auto* histo_simu = tmp_simu_thresh->GetHistogram();
        histo_simu->SetName(Form("simu_threasholds_%i", i));
        //histo_simu->Scale(1./(histo_simu->Integral(histo_simu->FindBin(0), histo_simu->FindBin(150))));
        histo_simu->Scale(1./(histo_simu->GetBinContent(histo_simu->GetMaximumBin())));

        auto* histo_diff = (TH1D*)histo_data->Clone(Form("probability_%i", i));
        histo_diff->Divide(histo_simu);
        histo_diff->SetLineColor(kGreen);
        output_histos.emplace(i, histo_diff);
//        histo_diff->Multiply(histo_simu);
//        histo_diff->Add(histo_simu, 1);
//        histo_diff->Scale(1./2.);


        if (show_canvas) {
            auto* cv = new TCanvas();
            histo_simu->SetMarkerColor(kRed);
            histo_simu->SetLineColor(kRed);
            histo_simu->Draw("histo");
            histo_data->SetMarkerColor(kBlue);
            histo_data->SetLineColor(kBlue);
            histo_data->Draw("same histo");
            histo_diff->SetMarkerColor(kGreen);
            histo_diff->SetLineColor(kGreen);
            histo_diff->Draw("same histo");
            cv->WaitPrimitive();
            delete cv;
        }
    }

    TFile output_file("sampling_histos.root", "recreate");
    for (const auto &it: output_histos){
        it.second->Write();
    }
    output_file.Write();
    output_file.Close();

}

double Selector::ComputeDoppler(const TVector3 & vec, const double & en){
    return en/(sqrt(1-beta*beta)*(1. + beta*vec.CosTheta()));
}

double Selector::ComputeDoppler(const TVector3 & vec,const TVector3 & pos, const double & en){
    return en/(sqrt(1-beta*beta)*(1. + beta*(vec-pos).CosTheta()));
}
