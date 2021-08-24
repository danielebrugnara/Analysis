void AlignEx(){
    auto* file = new TFile("../../DataAnalyzed/sum.root", "read");
    auto* tree = (TTree*) file->Get("AnalyzedTreeDeuterons");
    
    std::map<int, double> shift;

    shift[1] = 0.;
    shift[3] = 0.9;
    shift[4] = 0.5;
    shift[5] = 0.6;
    shift[7] = 0.2;
    shift[11] = 0.;

    std::string command;

    for (const auto& it:shift){
        if (!command.empty()) command += " + ";
        command += Form(" (MugastData.Ex -%f)*(MugastData.MG == %i)", it.second, it.first);
    }
    command += ">>h(80,-10, 10)";

    std::string conditions;

    //conditions += " MugastData.MG == 11";

    tree->Draw(command.c_str(), conditions.c_str());

}
