void makeOtherPlots(){
    TFile f("sum.root");
    TTree *tree = (TTree*)f.Get("AnalyzedTree");
    TFile* f_out = new TFile("outOther.root", "recreate");
    for (int i=0;i<32;++i){
        auto* histo = new TH2D(Form("h_corrected_%i", i),Form("h%i", i),1000, 2, 3.1415, 1000, 1.5, 9);
        std::cout << "h_corrected" << i << std::endl;
        tree->Draw(Form("MugastData.E_corrected[0][%i]:MugastData.EmissionDirection_corrected[%i].Theta()>>h_corrected%i", i, i, i), "MugastData.M ==2 && MugastData.Z == 1 && VamosData.id_Z == 19 && VamosData.id_M == 47");
        histo->Write();
        
    }
    for (int i=0;i<32;++i){
        auto* histo = new TH2D(Form("h_uncentered_%i", i),Form("h%i", i),1000, 2, 3.1415, 1000, 1.5, 9);
        std::cout << "h_uncentered" << i << std::endl;
        tree->Draw(Form("MugastData.E_uncentered[0][%i]:MugastData.EmissionDirection_uncentered[%i].Theta()>>h_uncentered%i", i, i, i), "MugastData.M ==2 && MugastData.Z == 1 && VamosData.id_Z == 19 && VamosData.id_M == 47");
        histo->Write();
        
    }
    f_out->Write();
    f_out->Close();
    

}
