void AlignEx(){
    auto* file = new TFile("../../DataAnalyzed/sum.root", "read");
    auto* tree = (TTree*) file->Get("AnalyzedTreeDeuterons");
    
    std::map<int, double> shift;

    ////Daniele
    //shift[1] = 0.;
    //shift[3] = 0.9;
    //shift[4] = 0.5;
    //shift[5] = 0.6;
    //shift[7] = 0.2;
    //shift[11] = 0.;

    ////Gottardo
    //shift[1] = 0.5;
    //shift[3] = -0.5;
    //shift[4] = -0.3;
    //shift[5] = -0.5;
    //shift[7] = 0.0;
    //shift[11] = 0.4;

    //Daniele 2
    //shift[1] =  0.8;
    //shift[3] =  1.;
    //shift[4] =  0.5;
    //shift[5] =  0.8;
    //shift[7] =  0.;
    //shift[11] = 0.;
    
    //No shift
    shift[1] =  0.;
    shift[3] =  0.;
    shift[4] =  0.;
    shift[5] =  0.;
    shift[7] =  0.;
    shift[11] = 0.;

    std::string command;

    for (const auto& it:shift){
        if (!command.empty()) command += " + ";
        command += Form(" (MugastData.Ex -%f)*(MugastData.MG == %i)", it.second, it.first);
    }
    std::string command2 = command;
    command += ">>h(65,-12, 8)";

    std::string conditions;
    //conditions += " MugastData.MG == 11";

    tree->Draw(command.c_str(), conditions.c_str());


    command2 += ":MugastData.EmissionDirection.Theta()";
    command2 += ">>h2(200, 1.5, 3.1415, 80, -10, 10)";

    (new TCanvas())->cd();
    tree->Draw(command2.c_str(), conditions.c_str());

}
