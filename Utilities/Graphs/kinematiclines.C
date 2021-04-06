void kinematiclines(){
    std::vector<double> ex = {0, 0.360, 2.02};
    std::vector<int> ls = {0, 2, 3};
    std::vector<int> colors2 = {kRed, kRed+1, kRed+2};
    std::vector<int> colors1 = {1, 12, 15};
    std::vector<int> styles = {1, 6, 10};
    std::vector<const char*> files = {"files/46Ar3Hed47K_0keV_s12.dat", "files/46Ar3Hed47K_360keV_d32.dat", "files/46Ar3Hed47K_2020keV_f72.dat"};
    std::vector<NPL::Reaction*> reactions;
    std::vector<TGraph*> graphs;
    std::vector<TGraph*> graphscmlab;
    std::vector<TGraph*> graphsangdist;

    auto* g1 = new TMultiGraph();
    auto* g2 = new TMultiGraph();
    auto* gcmlab1 = new TMultiGraph();
    auto* gangdist1 = new TMultiGraph();
    for(int i=0; i<ex.size(); ++i){
        reactions.push_back(new NPL::Reaction("46Ar(3He,d)47K@373.7"));
        reactions.back()->SetExcitation3(ex[i]);
        
        //kinematics3
        graphs.push_back(reactions.back()->GetKinematicLine3());
        graphs.back()->SetLineColor(colors1[i]);
        graphs.back()->SetLineStyle(styles[i]);
        graphs.back()->SetLineWidth(3);
        graphs.back()->SetMarkerColor(colors1[i]);
        graphs.back()->SetName(Form("d with ex = %3.2f MeV", ex[i]));
        g1->Add(graphs.back());
        
        //kinematics4
        graphs.push_back(reactions.back()->GetKinematicLine4());
        graphs.back()->SetLineColor(colors2[i]);
        graphs.back()->SetLineStyle(styles[i]);
        graphs.back()->SetLineWidth(3);
        graphs.back()->SetMarkerColor(colors2[i]);
        graphs.back()->SetName(Form("^{47}K with ex = %3.2f MeV", ex[i]));
        g2->Add(graphs.back());

        //cm lab
        graphscmlab.push_back(reactions.back()->GetThetaLabVersusThetaCM());
        graphscmlab.back()->SetLineColor(colors2[i]);
        graphscmlab.back()->SetLineStyle(styles[i]);
        graphscmlab.back()->SetLineWidth(3);
        graphscmlab.back()->SetMarkerColor(colors2[i]);
        graphscmlab.back()->SetName(Form("d with ex = %3.2f MeV", ex[i]));
        gcmlab1->Add(graphscmlab.back());

        //anggistr
        graphsangdist.push_back(new TGraph(files[i]));
        graphsangdist.back()->SetLineColor(colors2[i]);
        graphsangdist.back()->SetLineStyle(styles[i]);
        graphsangdist.back()->SetLineWidth(3);
        graphsangdist.back()->SetMarkerColor(colors2[i]);
        graphsangdist.back()->SetName(Form("d with L = %i", ls[i]));
        gangdist1->Add(graphsangdist.back());
    }

    

    auto* cv = new TCanvas();
    TPad *p1 = new TPad("p1", "", 0, 0, 1, 1);
    //p1->SetGrid();
    TPad *p2 = new TPad("p2", "", 0, 0, 1, 1);
    p2->SetFillStyle(4000); // will be transparent

    p1->Draw();
    p1->cd();
    g1->Draw("AL");
    g1->GetXaxis()->SetRangeUser(0, 180);
    g1->GetXaxis()->SetTitle("Angle [deg]");
    g1->GetYaxis()->SetTitle("Energy [MeV]");
    g1->GetXaxis()->SetLabelSize(0.05);
    g1->GetYaxis()->SetLabelSize(0.05);
    g1->GetXaxis()->SetTitleSize(0.05);
    g1->GetXaxis()->SetTitleOffset(0.9);
    g1->GetYaxis()->SetTitleOffset(0.9);
    g1->GetYaxis()->SetTitleSize(0.05);
    gPad->Update();

    p2->Draw();
    p2->cd();
    g2->Draw("LP");
    gPad->Update();

  //Double_t xmin = g2->GetHistogram()->GetXaxis()->GetXmin();
  //Double_t xmax = g2->GetHistogram()->GetXaxis()->GetXmax();
  Double_t xmin = g2->GetHistogram()->GetXaxis()->GetXmin();
  Double_t xmax = g2->GetHistogram()->GetXaxis()->GetXmax();
  Double_t dx = (xmax - xmin) / 0.8; // 10 percent margins left and right
#if 1 /* 0 or 1 */
  Double_t ymin = g2->GetHistogram()->GetMinimum();
  Double_t ymax = g2->GetHistogram()->GetMaximum();
#else /* 0 or 1 */
  Double_t ymin = g2->GetHistogram()->GetYaxis()->GetXmin();
  Double_t ymax = g2->GetHistogram()->GetYaxis()->GetXmax();
#endif /* 0 or 1 */
  Double_t dy = (ymax - ymin) / 0.8; // 10 percent margins top and bottom
  //p2->Range(xmin-0.1*dx, ymin-0.1*dy, xmax+0.1*dx, ymax+0.1*dy);
  p2->Range(xmin-0.1*dx, ymin-0.1*dy, xmax+0.1*dx, ymax+0.1*dy);
  p2->Draw();
  p2->cd();
  g2->Draw("LP");
  gPad->Update();
  
  TGaxis *xaxis = new TGaxis(xmin, ymax, xmax, ymax, xmin, xmax, 510, "-L");
  xaxis->SetLineColor(kRed);
  xaxis->SetLabelColor(kRed);
  xaxis->SetLabelSize(0.05);
  xaxis->Draw();
  gPad->Update();
  
  TGaxis *yaxis = new TGaxis(xmax, ymin, xmax, ymax, ymin, ymax, 510, "+L");
  yaxis->SetLineColor(kRed);
  yaxis->SetLabelColor(kRed);
  yaxis->SetLabelSize(0.05);
  yaxis->Draw();
  gPad->Update();


{
  TLegend *leg = new TLegend(0.5, 0.65, 0.9, 0.9);
  leg->SetFillColor(0);
  leg->SetTextSize(0.036);
  for(const auto& it: graphs){
    leg->AddEntry(it, it->GetName(), "L");
  }
  
  leg->Draw();
  gPad->Update();
}
    
  cv->SaveAs("kinematics.pdf");

{
  auto* cv = new TCanvas();
  gcmlab1->GetXaxis()->SetRangeUser(0, 180);
  gcmlab1->GetXaxis()->SetTitle("Angle CM [deg]");
  gcmlab1->GetYaxis()->SetTitle("Angle Lab [deg]");
  gcmlab1->GetXaxis()->SetLabelSize(0.05);
  gcmlab1->GetYaxis()->SetLabelSize(0.05);
  gcmlab1->GetXaxis()->SetTitleSize(0.05);
  gcmlab1->GetXaxis()->SetTitleOffset(0.9);
  gcmlab1->GetYaxis()->SetTitleOffset(0.9);
  gcmlab1->GetYaxis()->SetTitleSize(0.05);
  gcmlab1->Draw("AL");
  TLegend *leg = new TLegend(0.5, 0.75, 0.9, 0.9);
  leg->SetFillColor(0);
  leg->SetTextSize(0.036);
  for(const auto& it: graphscmlab){
    leg->AddEntry(it, it->GetName(), "L");
  }
  
  leg->Draw();
  cv->SaveAs("kinematics_cm_lab.pdf");
}

{
  auto* cv = new TCanvas();
  cv->SetLogy();
  gangdist1->GetXaxis()->SetRangeUser(0, 180);
  gangdist1->GetXaxis()->SetTitle("Angle CM [deg]");
  gangdist1->GetYaxis()->SetTitle("d#sigma/d#Omega [mb/sr]");
  gangdist1->GetXaxis()->SetLabelSize(0.05);
  gangdist1->GetYaxis()->SetLabelSize(0.05);
  gangdist1->GetXaxis()->SetTitleSize(0.05);
  gangdist1->GetXaxis()->SetTitleOffset(0.9);
  gangdist1->GetYaxis()->SetTitleOffset(0.9);
  gangdist1->GetYaxis()->SetTitleSize(0.05);
  gangdist1->Draw("AL");
  TLegend *leg = new TLegend(0.5, 0.75, 0.9, 0.9);
  leg->SetFillColor(0);
  leg->SetTextSize(0.036);
  for(const auto& it: graphsangdist){
    leg->AddEntry(it, it->GetName(), "L");
  }
  
  leg->Draw();
  cv->SaveAs("kinematics_angdistr.pdf");
}

}
