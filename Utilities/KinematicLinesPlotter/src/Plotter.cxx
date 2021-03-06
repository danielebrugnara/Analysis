#include "Plotter.h"

Plotter::Plotter():
        style("Modern", "Default"),
        beam_energy(400*UNITS::MeV)
{
    style.SetPalette(1);
    //style.SetPad;
    style.SetTitleW(0.9);
    style.cd();
    style.SetCanvasPreferGL(kTRUE);
    reactions.emplace("M48_Z19_m1_z1", new ReactionReconstruction2body<long double>(
            ReactionReconstruction2body<long double>::ReactionInput2body({
                                                                                 ReactionFragment::FragmentSettings(46, 18, 0, beam_energy, 0),
                                                                                 ReactionFragment::FragmentSettings(3, 2, 0, 0, 0),
                                                                                 ReactionFragment::FragmentSettings(1, 1, 0, 0, 0),
                                                                                 ReactionFragment::FragmentSettings(48, 19, 0, 0, 0)})));

    reactions.emplace("M47_Z19_m2_z1", new ReactionReconstruction2body<long double>(
            ReactionReconstruction2body<long double>::ReactionInput2body({
                                                                                 ReactionFragment::FragmentSettings(46, 18, 0, beam_energy, 0),
                                                                                 ReactionFragment::FragmentSettings(3, 2, 0, 0, 0),
                                                                                 ReactionFragment::FragmentSettings(2, 1, 0, 0, 0),
                                                                                 ReactionFragment::FragmentSettings(47, 19, 0,0, 0)})));

    reactions.emplace("M45_Z18_m4_z2", new ReactionReconstruction2body<long double>(
            ReactionReconstruction2body<long double>::ReactionInput2body({
                                                                                 ReactionFragment::FragmentSettings(46, 18, 0, beam_energy, 0),
                                                                                 ReactionFragment::FragmentSettings(3, 2, 0, 0, 0),
                                                                                 ReactionFragment::FragmentSettings(4, 2, 0, 0, 0),
                                                                                 ReactionFragment::FragmentSettings(45, 18, 0,0, 0)})));

    reactions.emplace("M46_Z18_m3_z2", new ReactionReconstruction2body<long double>(
            ReactionReconstruction2body<long double>::ReactionInput2body({
                                                                                 ReactionFragment::FragmentSettings(46, 18, 0, beam_energy, 0),
                                                                                 ReactionFragment::FragmentSettings(3, 2, 0, 0, 0),
                                                                                 ReactionFragment::FragmentSettings(3, 2, 0, 0, 0),
                                                                                 ReactionFragment::FragmentSettings(46, 18, 0,0, 0)})));

    reactions.emplace("M43_Z18_m6_z2", new ReactionReconstruction2body<long double>(
            ReactionReconstruction2body<long double>::ReactionInput2body({
                                                                                 ReactionFragment::FragmentSettings(46, 18, 0, beam_energy, 0),
                                                                                 ReactionFragment::FragmentSettings(3, 2, 0, 0, 0),
                                                                                 ReactionFragment::FragmentSettings(6, 2, 0, 0, 0),
                                                                                 ReactionFragment::FragmentSettings(43, 18, 0,0, 0)})));
    for (auto& it: reactions){
        it.second->SetPrecision(1E-4);
    }
}

void Plotter::plotVamosAcceptance(const std::string& file_name){

    std::map<std::string, std::pair<ReactionReconstruction2body<long double>*, TF1*>> graphs;
    for (auto& it: reactions){
        graphs.emplace(it.first, std::make_pair(it.second,nullptr));
    }

    int cnt=0;
    TMultiGraph multigraph;
    multigraph.SetTitle("Vamos acceptance;Angle[rad];Brho[Tm];");
    std::vector<bool> opts = {true, false};
    std::vector<int> charge_states = {13, 14, 15, 16};
//    std::vector<int> charge_states = {15};

    TLegend legend(0.7, 0.5, 1, 1);

    double brhoLost = 0.36;
    for(auto& it: graphs){
        auto* reaction = it.second.first;
        reaction->ChooseFixed(4);
        reaction->ChooseExFixed(4);
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
                it.second.second->SetLineColorAlpha(colors[cnt], ((double)it_q)/charge_states.back()-0.65);
                auto* tmp = new TGraph(it.second.second, "");
                multigraph.Add(tmp);

                if(it_opts)
                    legend.AddEntry(tmp, it.second.first->Get_Name().c_str(), "l");

                cnt++;
            }
        }
    }

    std::vector<double> xAcceptance = {0, 0.1, 0.1, 0, 0};
    std::vector<double> yAcceptance = {0.8, 0.8, 1.15, 1.15, 0.8};
    TPolyLine vamosAcceptance(xAcceptance.size(), &xAcceptance[0], &yAcceptance[0]);
    vamosAcceptance.SetFillColorAlpha(kRed, 0.13);
    vamosAcceptance.SetLineColorAlpha(kRed, 0.7);
    //vamosAcceptance.SetFillColor(kRed);
    vamosAcceptance.SetLineWidth(4);

    TCanvas cv;
    multigraph.Draw("al");
    multigraph.GetXaxis()->SetRangeUser(0, 1.2);
    vamosAcceptance.Draw("f");
    vamosAcceptance.Draw();
    legend.Draw();
    cv.Draw();
    cv.WaitPrimitive();
    cv.SaveAs("vamos_acc.pdf");
}

void Plotter::plotMugastAcceptance(const std::string& file_name) {
    std::map<std::string, std::pair<ReactionReconstruction2body<long double>*, TF1*>> graphs;
    for (auto& it: reactions){
        graphs.emplace(it.first, std::make_pair(it.second,nullptr));
    }

    int cnt=0;
    TMultiGraph multigraph;
    multigraph.SetTitle("Mugast acceptance;Angle[rad];Energy[MeV];");

    TLegend legend(0.7, 0.5, 1, 1);
    for(auto& it: graphs){
        auto* reaction = it.second.first;
        reaction->ChooseFixed(3);
        reaction->ChooseExFixed(3);
        it.second.second = new TF1(it.first.c_str(),
                                   [&reaction](double *x, double *y) {
                                       reaction->Set_Theta(*x, true);
                                       //return (double) reaction->GetFreeFragment().Get_P().Vect().Mag(); //Momentum
                                       return (double) reaction->GetFixedFragment().Get_Ek(); //Energy
                                   },
                                   0.,
                                   reaction->Get_ThetaMax(),
                                    //UNITS::CONSTANTS::pi,
                                   0);
        it.second.second->SetNpx(10000);
        it.second.second->SetLineColorAlpha(colors[cnt], 0.4);
        auto* tmp = new TGraph(it.second.second, "");
        multigraph.Add(tmp);
        legend.AddEntry(tmp, it.second.first->Get_Name().c_str(), "l");

        if (reaction->Get_ThetaMax()<UNITS::CONSTANTS::pi) {
            auto* tmp = new TF1((it.first + "_2").c_str(),
                    [&reaction](double *x, double *y) {
                        reaction->Set_Theta(*x, false);
                        //return (double) reaction->GetFreeFragment().Get_P().Vect().Mag(); //Momentum
                        return (double) reaction->GetFixedFragment().Get_Ek(); //Energy
                    },
                0.,
                    reaction->Get_ThetaMax(),
                    //UNITS::CONSTANTS::pi,
                    0);
            it.second.second->SetNpx(10000);
            it.second.second->SetLineColorAlpha(colors[cnt], 0.4);
            it.second.second->SetLineWidth(10);
            multigraph.Add(new TGraph(tmp));
        }
        cnt++;
    }

    std::vector<double> xAcceptanceMM = {0.18, 0.8, 0.8, 0.18, 0.18};
    std::vector<double> yAcceptanceMM = {1, 1, 137, 137, 1};
    TPolyLine acceptanceMM(xAcceptanceMM.size(), &xAcceptanceMM[0], &yAcceptanceMM[0]);
    acceptanceMM.SetFillColorAlpha(kRed, 0.13);
    acceptanceMM.SetLineColorAlpha(kRed, 0.7);
    //acceptanceMM.SetFillColor(kRed);
    acceptanceMM.SetLineWidth(4);

    std::vector<double> xAcceptanceMG = {2, 2.79, 2.79, 2, 2};
    std::vector<double> yAcceptanceMG = {1, 1, 29.4, 29.4, 1};
    TPolyLine acceptanceMG(xAcceptanceMG.size(), &xAcceptanceMG[0], &yAcceptanceMG[0]);
    acceptanceMG.SetFillColorAlpha(kBlue, 0.13);
    acceptanceMG.SetLineColorAlpha(kBlue, 0.7);
    //acceptanceMG.SetFillColor(kRed);
    acceptanceMG.SetLineWidth(4);

    std::vector<double> xAcceptanceMGAnular = {2.83, 2.985, 2.985, 2.83, 2.83};
    std::vector<double> yAcceptanceMGAnular = {0.8, 0.8, 29.4, 29.4, 0.8};
    TPolyLine acceptanceMGAnular(xAcceptanceMGAnular.size(), &xAcceptanceMGAnular[0], &yAcceptanceMGAnular[0]);
    acceptanceMGAnular.SetFillColorAlpha(kBlue, 0.13);
    acceptanceMGAnular.SetLineColorAlpha(kBlue, 0.7);
    //acceptanceMGAnular.SetFillColor(kRed);
    acceptanceMGAnular.SetLineWidth(4);

    TFile infile(file_name.c_str(), "read");

    TCanvas cv;
    multigraph.Draw("al");
    //multigraph.GetXaxis()->SetRangeUser(0, 1.2);
    acceptanceMM.Draw("f");
    acceptanceMM.Draw();
    acceptanceMG.Draw("f");
    acceptanceMG.Draw();
    acceptanceMGAnular.Draw("f");
    acceptanceMGAnular.Draw();
    legend.Draw();
    if (infile.IsOpen()){
        std::cout << "Plotting also histogram\n";
        infile.Get("Graph")->Draw("same");
    }

    cv.Draw();
    cv.WaitPrimitive();
    cv.SaveAs("mugast_acc.pdf");
    TFile outfile("mugast.root", "recreate");
    outfile.cd();
    for (const auto&it: graphs){
        it.second.second->Write("*");
    }
    cv.Write();
    outfile.Write();
    outfile.Close();

}


