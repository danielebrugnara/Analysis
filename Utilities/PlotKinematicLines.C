#include <NPReaction.h>

void PlotKinematicLines(){
    std::string en = "387";
    std::vector<double> ex = {0, 0.360, 2.020, 2.287, 3.35};
    std::vector<TGraph*> graphs;
    TFile *ff = new TFile("Lines.root", "recreate");
    for(const auto & ee:ex){
        auto *reaction = new NPL::Reaction("46Ar(3He,d)47K@"+en);
        reaction->SetExcitation4(ee);
        graphs.push_back(new TGraph());
        graphs.back()->SetName(("Ar46_"+en+"_"+std::to_string(ee)).c_str());
        graphs.back()->SetTitle(("46Ar(3He,d)47K@"+en+"_"+std::to_string(ee)).c_str());
        auto *line = reaction->GetKinematicLine3();
        for (int ii=0; ii<line->GetN(); ++ii){
            double x, y;
            line->GetPoint(ii,x, y);
            graphs.back()->SetPoint(ii,x*3.14159265359/180., y);
        }
        graphs.back()->Write();
    }
    ff->Write();
    ff->Close();
}
