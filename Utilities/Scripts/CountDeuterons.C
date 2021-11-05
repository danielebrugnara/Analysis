void CountDeuterons(){
    //Create with: AnalyzedTreeDeuterons->Draw("BeamEnergy:IceThicknessFront", "(Time<313 || Time>315) && (Time<278 || Time>283) && (Time<315 || Time>318)")
    //and then save the resulting graph 

    TFile* file = new TFile("deuterons_thickness_beam.root");
    auto* graph = (TGraph*) file->Get("Graph");
    double start{15E-6};
    double end{55E-6};
    double step{5E-6};
    std::map<double, double> thickness_to_energy;
    std::map<double, double> thickness_to_nr;

    for(double val{start}; val<=end; val += step){
        thickness_to_energy[val] = graph->Eval(val);
        thickness_to_nr[val] = 0;
    }

    for(int i{0}; i<graph->GetN(); ++i){
        double pt = graph->GetPointX(i);
        std::cout <<  "pt : " << pt << std::endl;
        double val{0};
        for (const auto& it: thickness_to_nr){
            std::cout << "it " << it.first << std::endl;
            std::cout << "f+s " << it.first+step/2. << std::endl;
            if (pt>it.first + step/2.){
                std::cout << " cal\n";
                val = it.first;
            }else{
                break;
            }
        }
        thickness_to_nr[val] = thickness_to_nr[val]+1;
    }
    for (const auto& it: thickness_to_nr){
        std::cout <<"thickness: " <<it.first << "\t\tenergy: " << std::round(thickness_to_energy[it.first]) << "\t\tcounts: " << it.second << std::endl;
    }
    
}
