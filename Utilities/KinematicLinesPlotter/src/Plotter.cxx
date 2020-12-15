#include "Plotter.h"


void Plotter::plotLines(){
    TStyle style;
    style.cd();
    style.SetCanvasPreferGL(kTRUE);


    std::map<std::string, std::pair<ReactionReconstruction2body<long double>*, TF1*>> graphs;
    double beam_energy = 400*UNITS::MeV;

    graphs.emplace("M48_Z19_m1_z1", std::make_pair(new ReactionReconstruction2body<long double>(
            ReactionReconstruction2body<long double>::ReactionInput2body({
                                                                                 ReactionFragment::FragmentSettings(46, 18, 0, beam_energy, 0),
                                                                                 ReactionFragment::FragmentSettings(3, 2, 0, 0, 0),
                                                                                 ReactionFragment::FragmentSettings(1, 1, 0, 0, 0),
                                                                                 ReactionFragment::FragmentSettings(48, 19, 0, 0, 0)})),
                                                   nullptr));

    graphs.emplace("M47_Z19_m2_z1", std::make_pair(new ReactionReconstruction2body<long double>(
            ReactionReconstruction2body<long double>::ReactionInput2body({
                                                                                 ReactionFragment::FragmentSettings(46, 18, 0, beam_energy, 0),
                                                                                 ReactionFragment::FragmentSettings(3, 2, 0, 0, 0),
                                                                                 ReactionFragment::FragmentSettings(2, 1, 0, 0, 0),
                                                                                 ReactionFragment::FragmentSettings(47, 19, 0,0, 0)})),
                                                   nullptr));

    graphs.emplace("M45_Z18_m4_z2", std::make_pair(new ReactionReconstruction2body<long double>(
            ReactionReconstruction2body<long double>::ReactionInput2body({
                                                                                 ReactionFragment::FragmentSettings(46, 18, 0, beam_energy, 0),
                                                                                 ReactionFragment::FragmentSettings(3, 2, 0, 0, 0),
                                                                                 ReactionFragment::FragmentSettings(4, 2, 0, 0, 0),
                                                                                 ReactionFragment::FragmentSettings(45, 18, 0,0, 0)})),
                                                   nullptr));

    graphs.emplace("M46_Z18_m3_z2", std::make_pair(new ReactionReconstruction2body<long double>(
            ReactionReconstruction2body<long double>::ReactionInput2body({
                                                                                 ReactionFragment::FragmentSettings(46, 18, 0, beam_energy, 0),
                                                                                 ReactionFragment::FragmentSettings(3, 2, 0, 0, 0),
                                                                                 ReactionFragment::FragmentSettings(3, 2, 0, 0, 0),
                                                                                 ReactionFragment::FragmentSettings(46, 18, 0,0, 0)})),
                                                   nullptr));



    int cnt=1;
    TMultiGraph multigraph;
    multigraph.SetTitle("Vamos acceptance;Angle[rad];Brho[Tm];");
    std::vector<bool> opts = {true, false};
    std::vector<int> charge_states = {13, 14, 15, 16};
//    std::vector<int> charge_states = {15};

    TLegend legend(0.7, 0.5, 1, 1);
    for(auto& it: graphs){
        cnt++;
        auto* reaction = it.second.first;
        reaction->ChooseFixed(4);
        reaction->ChooseExFixed(4);
        double brhoLost = 0.36;
        for(const auto&it_opts: opts) {
            for (const auto &it_q: charge_states) {

                it.second.second = new TF1(it.first.c_str(),
                                           [&reaction, &it_opts, &it_q, &brhoLost](double *x, double *y) {
                                               reaction->Set_Theta(*x, it_opts);
                                               //return (double) reaction->GetFreeFragment().Get_P().Vect().Mag(); //Momentum
                                               //return (double) reaction->GetFixedFragment().Get_Ek(); //Energy
                                               return (double) reaction->GetFixedFragment().GetBrhoFromEk_Q(
                                                       reaction->GetFixedFragment().Get_Ek(),
                                                       it_q
                                               )-brhoLost;
                                           },
                                           0.,
                                           reaction->Get_ThetaMax(),
                        //UNITS::CONSTANTS::pi,
                                           0);
                it.second.second->SetNpx(1000);
                it.second.second->SetLineColorAlpha(cnt, ((double)it_q)/charge_states.back()-0.65);
                auto* tmp = new TGraph(it.second.second, "");
                multigraph.Add(tmp);
                if(it_opts)
                    legend.AddEntry(tmp, it.first.c_str(), "l");
            }
        }
    }

    std::vector<double> xAcceptance = {0, 0.1, 0.1, 0, 0};
    std::vector<double> yAcceptance = {0.8, 0.8, 1.15, 1.15, 0.8};
    TPolyLine vamosAcceptance(xAcceptance.size(), &xAcceptance[0], &yAcceptance[0]);
    vamosAcceptance.SetFillColorAlpha(kRed, 0.13);
    //vamosAcceptance.SetFillColor(kRed);
    vamosAcceptance.SetLineColor(2);
    vamosAcceptance.SetLineWidth(4);

    TCanvas cv;
    multigraph.Draw("al");
    vamosAcceptance.Draw("f");
    vamosAcceptance.Draw();
    legend.Draw();
    cv.Draw();
    cv.WaitPrimitive();
    cv.SaveAs("Acceptance.root", "recreate");
}
