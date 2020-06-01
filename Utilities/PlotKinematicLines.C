#include <NPReaction.h>
#include <TFile.h>

void PlotKinematicLines(){
    std::vector<std::string> ens = {"378","334"};
    std::vector<double> ex4 = {0, 0.360, 2.020, 2.287, 3.35};
    std::vector<double> ex3 = {0, 2.224566};
    std::vector<TGraph*> graphs;
    TFile *ff = new TFile("Lines.root", "recreate");
    for (const auto & en: ens){
         for(const auto & ee4: ex4){
         for(const auto & ee3: ex3){
             auto *reaction = new NPL::Reaction("46Ar(3He,d)47K@"+en);
             reaction->SetExcitation4(ee4);
             reaction->SetExcitation3(ee3);
             graphs.push_back(new TGraph());
             graphs.back()->SetName(("Ar46_"+en+"_3_"+std::to_string(ee4)+"_"+std::to_string(ee3)).c_str());
             graphs.back()->SetTitle(("46Ar(3He,d)47K@"+en+"_3_"+std::to_string(ee4)+"_"+std::to_string(ee3)).c_str());
             auto *line = reaction->GetKinematicLine3();
             for (int ii=0; ii<line->GetN(); ++ii){
                 double x, y;
                 line->GetPoint(ii,x, y);
                 graphs.back()->SetPoint(ii,x*3.1415/180., y);
             }
             graphs.back()->Write();
             graphs.push_back(new TGraph());
             graphs.back()->SetName(("Ar46_"+en+"_4_"+std::to_string(ee4)+"_"+std::to_string(ee3)).c_str());
             graphs.back()->SetTitle(("46Ar(3He,d)47K@"+en+"_4_"+std::to_string(ee4)+"_"+std::to_string(ee3)).c_str());
             line = reaction->GetKinematicLine4();
             for (int ii=0; ii<line->GetN(); ++ii){
                 double x, y;
                 line->GetPoint(ii,x, y);
                 graphs.back()->SetPoint(ii,x*3.1415/180., y);
             }
             graphs.back()->Write();
         }
         }
    }
    ff->Write();
    ff->Close();
}
