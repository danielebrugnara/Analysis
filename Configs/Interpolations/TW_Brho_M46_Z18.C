void TW_Brho_M46_Z18(){

   const double x[2]={242.675,327.779};
   const double y[2]={1.0765,0.865234};
   auto* graph = new TGraph(2, x, y);

   auto* spl = new TSpline3("TW_Brho_M46_Z18", graph);

   auto* file = new TFile("TW_Brho_M46_Z18.root", "recreate");
   spl->Write();
   file->Write();
   file->Close();
}
