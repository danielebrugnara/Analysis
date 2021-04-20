void discr_distr_ang_likelihood_maximized()
{
//=========Macro generated from canvas: c1_n3/c1_n3
//=========  (Tue Apr 13 13:40:02 2021) by ROOT version 6.20/04

   gStyle->SetCanvasPreferGL(kTRUE);

   TCanvas *c1_n3 = new TCanvas("c1_n3", "c1_n3",838,525,700,500);
   gStyle->SetOptStat(0);
   c1_n3->Range(78.75,-22.05909,191.25,22.24457);
   c1_n3->SetFillColor(0);
   c1_n3->SetBorderMode(0);
   c1_n3->SetBorderSize(2);
   c1_n3->SetFrameBorderMode(0);
   c1_n3->SetFrameBorderMode(0);
   
   TH1D *discrepancies__3 = new TH1D("discrepancies__3","his_data",45,90,180);
   discrepancies__3->SetBinContent(15,10.20356);
   discrepancies__3->SetBinContent(16,0.8612431);
   discrepancies__3->SetBinContent(17,-6.878042);
   discrepancies__3->SetBinContent(18,1.963324);
   discrepancies__3->SetBinContent(19,-9.450924);
   discrepancies__3->SetBinContent(20,2.391376);
   discrepancies__3->SetBinContent(21,2.20421);
   discrepancies__3->SetBinContent(22,-8.824686);
   discrepancies__3->SetBinContent(23,-6.44874);
   discrepancies__3->SetBinContent(24,0.5881565);
   discrepancies__3->SetBinContent(25,-3.177985);
   discrepancies__3->SetBinContent(26,-1.038428);
   discrepancies__3->SetBinContent(27,-3.116309);
   discrepancies__3->SetBinContent(28,5.513052);
   discrepancies__3->SetBinContent(29,3.246677);
   discrepancies__3->SetBinContent(30,1.162735);
   discrepancies__3->SetBinContent(31,7.437853);
   discrepancies__3->SetBinContent(32,3.028239);
   discrepancies__3->SetBinContent(33,-2.295605);
   discrepancies__3->SetBinContent(34,2.875973);
   discrepancies__3->SetBinContent(35,4.231504);
   discrepancies__3->SetBinContent(37,-3.629008);
   discrepancies__3->SetBinContent(38,-4.239846);
   discrepancies__3->SetBinContent(39,-7.484384);
   discrepancies__3->SetBinContent(40,-0.3276769);
   discrepancies__3->SetBinContent(41,3.203727);
   discrepancies__3->SetBinError(15,5.922881);
   discrepancies__3->SetBinError(16,6.158025);
   discrepancies__3->SetBinError(17,6.140131);
   discrepancies__3->SetBinError(18,6.286133);
   discrepancies__3->SetBinError(19,6.570414);
   discrepancies__3->SetBinError(20,6.121563);
   discrepancies__3->SetBinError(21,5.775793);
   discrepancies__3->SetBinError(22,5.045463);
   discrepancies__3->SetBinError(23,4.535871);
   discrepancies__3->SetBinError(24,3.615228);
   discrepancies__3->SetBinError(25,3.159279);
   discrepancies__3->SetBinError(26,2.813465);
   discrepancies__3->SetBinError(27,2.458697);
   discrepancies__3->SetBinError(28,2.109188);
   discrepancies__3->SetBinError(29,1.930425);
   discrepancies__3->SetBinError(30,1.951735);
   discrepancies__3->SetBinError(31,1.880963);
   discrepancies__3->SetBinError(32,1.985389);
   discrepancies__3->SetBinError(33,2.064104);
   discrepancies__3->SetBinError(34,1.762235);
   discrepancies__3->SetBinError(35,1.32761);
   discrepancies__3->SetBinError(37,4.845183);
   discrepancies__3->SetBinError(38,5.808875);
   discrepancies__3->SetBinError(39,4.831633);
   discrepancies__3->SetBinError(40,5.08998);
   discrepancies__3->SetBinError(41,3.395839);
   discrepancies__3->SetEntries(575);

   Int_t ci;      // for color index setting
   TColor *color; // for color definition with alpha
   ci = TColor::GetColor("#009900");
   discrepancies__3->SetLineColor(ci);
   discrepancies__3->SetLineWidth(4);
   discrepancies__3->GetXaxis()->SetLabelFont(42);
   discrepancies__3->GetXaxis()->SetTitleOffset(1);
   discrepancies__3->GetXaxis()->SetTitleFont(42);
   discrepancies__3->GetYaxis()->SetLabelFont(42);
   discrepancies__3->GetYaxis()->SetTitleFont(42);
   discrepancies__3->GetZaxis()->SetLabelFont(42);
   discrepancies__3->GetZaxis()->SetTitleOffset(1);
   discrepancies__3->GetZaxis()->SetTitleFont(42);
   discrepancies__3->Draw("");
   
   TPaveText *pt = new TPaveText(0.4205158,0.9370168,0.5794842,0.995,"blNDC");
   pt->SetName("title");
   pt->SetBorderSize(0);
   pt->SetFillColor(0);
   pt->SetFillStyle(0);
   pt->SetTextFont(42);
   TText *pt_LaTex = pt->AddText("his_data");
   pt->Draw();
   c1_n3->Modified();
   c1_n3->cd();
   c1_n3->SetSelected(c1_n3);
}
