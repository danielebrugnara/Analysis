void GenerateOverlaps(){

    std::map<std::string, std::string> files;

    files["fort.20"] = "overlap_s12_-14.7785_sf659.dat";
    files["fort.22"] = "overlap_d32_-12.2673_sf0247.dat";
    files["fort.23"] = "overlap_f72_-8.6727_sf6882.dat";

    double step{0.0352};

    for (const auto& file: files){
        TGraph gr(file.second.c_str());
        std::ofstream ouf(file.first);

        int npt{0};
        std::vector<double> vals;
        for(double x{0}; x<30.; x+=step){
            npt++;
            vals.push_back(gr.Eval(x)/x);     
        }

        vals[0] = vals[1];

        ouf << "#input from " << file.second << std::endl;
        ouf << npt << "  " << step << " 0" << std::endl;
        for (const auto&it: vals){
            ouf << it << std::endl;
        }
        //two times, but this part is ignored
        for (const auto&it: vals){
            ouf << it << std::endl;
        }
        
    }
}
