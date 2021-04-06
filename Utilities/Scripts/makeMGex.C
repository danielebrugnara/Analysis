void makeMGex(){
    auto* f = new TFile("sum.root");
    auto* t = (TTree*) f->Get("AnalyzedTree");
    std::set<int> mg = {1, 3, 4, 5, 7};

    auto* out = new TFile("MGex.root", "recreate");

    for (const auto& it: mg){
        auto* h = new TH1D(Form("h%i", it), Form("h%i", it), 100, -10, 10);
        t->Draw(Form("MugastData.Ex[0]>>h%i", it), Form("MugastData.M ==2 && MugastData.Z == 1 && VamosData.id_Z == 19 && VamosData.id_M == 47 && MugastData.MG == %i", it), "");
        h->Write();
    }
    out->Write();
    out->Close();


}
