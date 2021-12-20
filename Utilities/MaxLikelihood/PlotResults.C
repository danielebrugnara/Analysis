void PlotResults(){
    std::vector<int>bins{60, 69, 82, 90, 100, 113, 129};
    std::map<int, std::pair<TGraph*,TGraph*>>gr;
    std::map<int, int>colors;

    colors[bins[0]] = kYellow+1;
    colors[bins[1]] = kRed+1;
    colors[bins[2]] = kGreen+1;
    colors[bins[3]] = kBlue+1;
    colors[bins[4]] = kMagenta+2;
    colors[bins[5]] = kOrange;
    colors[bins[6]] = kBlack;
    
    TMultiGraph* mg = new TMultiGraph();
    auto* leg = new TLegend(0.6, 0.5, 0.9, 0.9);

    for(const auto& it:bins){
        auto* fs= new TFile(Form("cont_%i_2.root", it), "read");
        auto* fm= new TFile(Form("min_%i.root", it), "read");
        gr[it].first = (TGraph*)fm->Get(Form("min_%i", it));
        gr[it].second = (TGraph*)fs->Get(Form("cont_%i", it));
        gr[it].first->SetLineColor(colors[it]);
        gr[it].first->SetMarkerColor(colors[it]);
        gr[it].first->SetMarkerStyle(8);
        gr[it].second->SetLineColorAlpha(colors[it], 0.5);
        gr[it].second->SetLineWidth(5);
        leg->AddEntry(gr[it].second, Form("Width: %.2f", (180.-90.)/it), "lp");
        mg->Add(gr[it].second);
    }
    mg->SetTitle(";C^{2}S[L=2]/C^{2}S[L=0];C^{2}S[L=3]/C^{2}S[L=0]");
    mg->GetXaxis()->SetLabelSize(0.05);
    mg->GetXaxis()->SetTitleSize(0.05);
    mg->GetYaxis()->SetLabelSize(0.05);
    mg->GetYaxis()->SetTitleSize(0.05);
    mg->GetYaxis()->SetTitleOffset(0.9);
    
    mg->Draw("a");
    for(const auto& it:bins){
        gr[it].first->Draw("same, p");
    }
    leg->Draw();
}
