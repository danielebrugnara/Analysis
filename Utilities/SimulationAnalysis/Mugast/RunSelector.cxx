#include "RunSelector.h"

RunSelector::RunSelector(std::string simu_file_name, std::string data_file_name){

    std::string serialization_name = "serialization.txt";
    std::ifstream infile(serialization_name);
    if (!infile.is_open()) {
        std::cout << "Getting simu\n";
        std::map<std::string, Selector::Transform> trsf;
        TFile *fileSimu = new TFile(simu_file_name.c_str(), "read");
        if (fileSimu == nullptr) std::cout << "File simu not found";
        TTree *treeSimu = (TTree *) fileSimu->Get("PhysicsTree");
        if (treeSimu == nullptr) std::cout << "Tree simu not found";
        Selector *selectorSimu = new Selector();
        std::cout << "Starting selector\n";
        treeSimu->Process(selectorSimu, "option");
        trsf.emplace("simu", selectorSimu->GetTsf());

        std::cout << "Getting data\n";
        TFile *fileData = new TFile(data_file_name.c_str(), "read");
        if (fileData == nullptr) std::cout << "File data not found";
        TTree *treeData = (TTree *) fileData->Get("AnalyzedTree");
        if (treeData == nullptr) std::cout << "Tree data not found";
        SelectorData *selectorData = new SelectorData();
        std::cout << "Starting selector\n";
        treeData->Process(selectorData, "option");
        trsf.emplace("data", selectorData->GetTsf());
        std::map<Selector::Id, double> threasholdStripE = selectorData->GetThsE();
        std::map<Selector::Id, double> threasholdStripT = selectorData->GetThsT();

        std::ofstream outfile(serialization_name);
        boost::archive::binary_oarchive oarch(outfile);
        oarch << trsf;
        oarch << threasholdStripE;
        oarch << threasholdStripT;
        outfile.close();
        infile.open(serialization_name);

    }

    std::map<std::string, Selector::Transform> trsf_read;
    std::map<Selector::Id, double> threashE_read;
    std::map<Selector::Id, double> threashT_read;
    boost::archive::binary_iarchive iarch(infile);
    iarch >> trsf_read;
    iarch >> threashE_read;
    iarch >> threashT_read;

    bool runSecondSelector{true};
    if (runSecondSelector){
        TFile *fileSimu = new TFile(simu_file_name.c_str(), "read");
        if (fileSimu == nullptr) std::cout << "File simu not found";
        TTree *treeSimu = (TTree *) fileSimu->Get("PhysicsTree");
        if (treeSimu == nullptr) std::cout << "Tree simu not found";
        Selector *selectorSimu = new Selector();
        selectorSimu->SetThreasholds(threashT_read, threashE_read);
        std::cout << "Starting selector\n";
        treeSimu->Process(selectorSimu, "option");
    }

    for (auto &itTsf: trsf_read) {
        for (auto &it: itTsf.second) {
            double ncounts{0};
            Selector::Pos maxKey;
            for (auto &it2: it.second) {
                if (it2.second > ncounts) {
                    maxKey = it2.first;
                    ncounts = it2.second;
                }
            }
            for (auto it2 = it.second.begin(); it2 != it.second.end();) {
                if (it2->first != maxKey){
                    it.second.erase(it2++);
                }
                else{
                    ++it2;
                }
            }
        }
    }

    //if (trsf_read == trsf) std::cout << "happy\n";
    std::map<std::string, Selector::Transform> tsfUnsorted = trsf_read;

    for (auto& itSimu: trsf_read["simu"]){
        if (itSimu.first.mg == 9 || itSimu.first.mg == 10)
            continue;


        auto itData = trsf_read["data"].find(itSimu.first);
        if (itData != trsf_read["data"].end()){
            //if(itSimu.first.mg == 11) {
            //    std::cout << "#########found: " <<
            //              " x: " << itSimu.first.x <<
            //              " y: " << itSimu.first.y <<
            //              " mg: " << itSimu.first.mg << std::endl;
            //}
            ////if (itData->second.begin()->first ==  itSimu.second.begin()->first)
            double distance{itData->second.begin()->first.dist(itSimu.second.begin()->first, 1E-3)};
            //std::cout << "Distance "  << distance << std::endl;
            //std::cout << "x "  <<itData->second.begin()->first.x << "  " <<  itSimu.second.begin()->first.x << std::endl;
            //std::cout << "y "  <<itData->second.begin()->first.y << "  " <<  itSimu.second.begin()->first.y << std::endl;
            //std::cout << "z "  <<itData->second.begin()->first.z << "  " <<  itSimu.second.begin()->first.z << std::endl;
            if (distance> 1E-9)
                throw std::runtime_error("error in distance\n");
            else{
                tsfUnsorted["simu"].erase(itSimu.first);
                tsfUnsorted["data"].erase(itData->first);
            }

        }else{
            //if(itSimu.first.mg == 11) {
            //    std::cout << "missing: " <<
            //              " x: " << itSimu.first.x <<
            //              " y: " << itSimu.first.y <<
            //              " mg: " << itSimu.first.mg << std::endl;
            //}
        }
    }

    std::map<Selector::Id, Selector::Id> simulationToData;
    for (const auto& itSimu: tsfUnsorted["simu"]){
        for(const auto& itData: tsfUnsorted["data"]){
            double distance{itData.second.begin()->first.dist(itSimu.second.begin()->first, 1E-3)};
            if (distance < 1E-9) {
                simulationToData.emplace(itSimu.first, itData.first);
                std::cout << "found correspondence!!\n";
                tsfUnsorted["simu"].erase(itSimu.first);
                tsfUnsorted["data"].erase(itSimu.first);
                std::cout << "mg : " << itSimu.first.mg << std::endl;
            }
        }
    }

    std::map<std::pair<int, int>, int> mappingDataToSimX;
    std::map<std::pair<int, int>, int> mappingDataToSimY;


    for (const auto&it: simulationToData){
        std::pair<int, int> pairX(it.second.mg, it.second.x);
        auto itX = mappingDataToSimX.find(pairX);
        if (itX == mappingDataToSimX.end()){
            mappingDataToSimX.emplace(pairX, it.first.x);
        }else{
            if (itX->second != it.first.x){
                std::cout << "one option is: "<< itX->first.second << " to " <<  itX->second << std::endl;
                std::cout << "other option is: "<< itX->first.second << " to " <<  it.first.x << std::endl;
                //throw std::runtime_error("something wrong\n");
            }
        }


        std::pair<int, int> pairY(it.second.mg, it.second.y);
        auto itY = mappingDataToSimY.find(pairY);
        if (itY == mappingDataToSimY.end()){
            mappingDataToSimY.emplace(pairY, it.first.y);
        }else{
            if (itY->second != it.first.y)
                throw std::runtime_error("something wrong\n");
        }
    }


    std::ofstream outMapping("outMappingDataToSimulation.txt");
    for (const auto&it: mappingDataToSimX){
        outMapping  << "MG" << it.first.first
                    << " x=" << it.first.second
                    <<  " -> "
                    << " x=" << it.second
                    << std::endl;
    }

    outMapping << "\n\n";

    for (const auto&it: mappingDataToSimY){
        outMapping  << "MG" << it.first.first
                    << " y=" << it.first.second
                    <<  " -> "
                    << " y=" << it.second
                    << std::endl;
    }

    std::map<std::pair<int, int>, bool> disabledX;
    std::map<std::pair<int, int>, bool> disabledY;

    std::vector<int> mgs{1,3,4,5,7,11};

    for(int mg:mgs){
        for(int i=0; i<128; ++i){
            disabledX.emplace(std::make_pair(mg,i+1), true);
            disabledY.emplace(std::make_pair(mg,i+1), true);
        }
    }

    for(const auto& it: trsf_read["data"]){
        auto findX = disabledX.find(std::make_pair(it.first.mg, it.first.x));
        if (findX != disabledX.end())
            disabledX.erase(findX);

        auto findY = disabledY.find(std::make_pair(it.first.mg, it.first.y));
        if (findY != disabledY.end())
            disabledY.erase(findY);
    }

    std::ofstream outMissingStrips("outMissingStrips.txt");
    for (const auto&it: disabledX){
        if (it.first.first != 10 || it.first.first != 9) {
            outMissingStrips << "MG" << it.first.first
                             << " x " << it.first.second
                             << std::endl;
        }
    }
    for (const auto&it: disabledY){
        if (it.first.first != 10 || it.first.first != 9) {
            outMissingStrips << "MG" << it.first.first
                             << " y " << it.first.second
                             << std::endl;
        }
    }

}
