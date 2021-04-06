void makePlots(){
    TFile f("sum.root");
    TTree *tree = (TTree*)f.Get("AnalyzedTree");
    TFile* f_out = new TFile("out.root", "recreate");
    for (int i=0;i<32;++i){
        auto* histo = new TH1D(Form("h_corrected_%i", i),Form("h%i", i),150, -15, 15);
        std::cout << "h_corrected" << i << std::endl;
        tree->Draw(Form("MugastData.Ex_corrected[0][%i]>>h_corrected%i", i, i), "MugastData.M ==2 && MugastData.Z == 1 && VamosData.id_Z == 19 && VamosData.id_M == 47");
        histo->Write();
        
    }
    for (int i=0;i<32;++i){
        auto* histo = new TH1D(Form("h_uncentered_%i", i),Form("h%i", i),150, -15, 15);
        std::cout << "h_uncentered" << i << std::endl;
        tree->Draw(Form("MugastData.Ex_uncentered[0][%i]>>h_uncentered%i", i, i), "MugastData.M ==2 && MugastData.Z == 1 && VamosData.id_Z == 19 && VamosData.id_M == 47");
        histo->Write();
        
    }
    f_out->Write();
    f_out->Close();
    

}
