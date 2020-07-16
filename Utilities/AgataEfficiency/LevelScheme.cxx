#include "LevelScheme.h"

LevelScheme::LevelScheme(const std::string & input_file_name):
        scheme(nullptr){
    ReadFile(input_file_name);
}

void LevelScheme::ReadFile(const std::string & input_file_name) {
    std::fstream fin(input_file_name, std::ios::in);

    if (!fin.is_open())
        throw std::runtime_error("csv file: "+input_file_name+" not found\n");

    std::string line;

    std::unordered_map<std::string, int> column;
    column["Elevel"] = 0;
    column["References"] = 1;
    column["SpinParity"] = 2;
    column["T12"] = 3;
    column["Egamma"] = 4;
    column["Intensity"] = 5;
    column["Multipolarity"] = 6;
    column["FinalLevel"] = 7;
    column["FinalSpinParity"] = 8;

    std::vector<std::vector<std::string>> data;
    while(std::getline(fin, line)){
        //std::cout << line << std::endl;
        auto row = ReadCSVRow(line);
        if (!row[0].compare(0, 8, "E(level)"))
            continue;
        data.push_back(row);
        //for (auto & field: row){
        //    std::cout << "field : " << field << "\t";
        //}
        //std::cout << std::endl;
    }
    //std::cout << "-------------------------\n";
    std::vector<Level> levels;
    std::vector<Gamma> gammas;
    std::unordered_map<double, int> energytoindex;

    levels.reserve(data.size());
    gammas.reserve(data.size());
    for(const auto & it: data){
        try {
            if (it[column["Elevel"]].empty())
                continue;
            levels.emplace_back(GetEnergyandError(it[column["Elevel"]]),
                                GetLifetimeandError(it[column["T12"]]),
                                GetSpinParity(it[column["SpinParity"]]));
            energytoindex[levels.back().Elevel.first] = levels.size()-1;

            std::vector<std::pair<double, double>> tmp_final_lev =
                    GetMultipleEnergiesandError(it[column["FinalLevel"]]);

            std::vector<std::pair<double, double>> tmp_gamma =
                    GetMultipleEnergiesandError(it[column["Egamma"]]);

            std::vector<std::pair<double, double>> tmp_intensities =
                    GetMultipleIntensitiesandError(it[column["Intensity"]]);
            auto branchings = ComputeBranching(tmp_intensities);

            std::vector<std::string> tmp_multipolarities =
                    GetMultipolariy(it[column["Multipolarity"]]);

            if (tmp_gamma.size()<tmp_final_lev.size())
                std::cerr << "final_lev\n";
            if (tmp_intensities.size()<tmp_final_lev.size())
                std::cerr << "intens\n";
            if (tmp_multipolarities.size()<tmp_final_lev.size())
                std::cerr << "multip\n";
            for (long unsigned ii=0; ii<tmp_final_lev.size(); ++ii){
                gammas.emplace_back(tmp_gamma[ii],
                                    std::make_pair( tmp_final_lev[ii].first,
                                                    levels.back().Elevel.first),
                                    std::make_pair( energytoindex[tmp_final_lev[ii].first],
                                                    energytoindex[levels.back().Elevel.first]),
                                    tmp_intensities[ii],
                                    branchings[ii],
                                    tmp_multipolarities[ii]
                );

            }

        }catch(const std::invalid_argument& err){
            std::cerr << "Invalid arg\n";
            std::cerr << err.what()<< "\n";
            std::cerr << "Offending level : " << it[column["Elevel"]] <<"\n";
        }catch(const std::bad_alloc& err) {
            std::cerr << "Bad alloc\n";
            std::cerr << err.what() << "\n";
            std::cerr << "Offending level : " << it[column["Elevel"]]  << "\n";
        }catch(const std::out_of_range& err) {
            std::cerr << "Out of range\n";
            std::cerr << err.what() << "\n";
            std::cerr << "Offending level : " << it[column["Elevel"]]  << "\n";
        }
    }

    std::vector<graphEdge<Gamma>> edges;

    edges.reserve(gammas.size());
    for (const auto & g:gammas)
        edges.push_back({g.Idxf_Idxi.second, g.Idxf_Idxi.first, g});

    scheme = new DiaGraph<Gamma, Level>(
            &edges[0],
            &levels[0],
            edges.size(),
            levels.size()
            );
    std::cout << *scheme;

    std::cout << "----------------------------------------------\n";
    for (int ii=0; ii<scheme->GetNumberNodes(); ++ii) {
        auto nod = scheme->getConsecutiveEdges(ii);
        for (const auto & it: nod){

            std::cout <<"~~" << *it.first <<" -> "<< *it.second<< std::endl;
        }
    }
}

std::pair<double, double> LevelScheme::GetEnergyandError(const std::string & st) {
    if (st.empty())
        return std::make_pair(0, 0);
    std::size_t found = st.find_first_of(' ');
    if (found != std::string::npos)
        return std::make_pair(std::stod(st.substr(0, found))*UNITS::keV,
                              std::stod(st.substr(found, st.size()))*UNITS::keV);
    else
        return std::make_pair(std::stod(st)*UNITS::keV,
                              0);
}

std::pair<double, double> LevelScheme::GetLifetimeandError(const std::string & st) {
    if (st == "STABLE")
        return std::make_pair(-1, 0);
    if (st.empty())
        return std::make_pair(0, 0);
    std::istringstream ss(st);
    double lifetime;
    double error = 0;
    std::string error_str;
    ss >> lifetime;
    std::string units;
    ss >> units;
    ss >> error_str;

    if (    error_str.find('+') != std::string::npos
            &&  error_str.find('-') != std::string::npos){
        std::size_t  min_found = error_str.find('-');
        error = (std::stod(error_str.substr(1, min_found))+
                 std::stod(error_str.substr(min_found+1)))/2.;
    }else{
        if (!error_str.empty())
            error = std::stod(error_str);
        else
            error = 0;
    }

    if (units =="s"){
        lifetime *=UNITS::second;
        error *=UNITS::second;
    }
    if (units =="ms") {
        lifetime *= UNITS::millisecond;
        error *= UNITS::millisecond;
    }
    if (units == "us") {
        lifetime *= UNITS::microsecond;
        error *= UNITS::microsecond;
    }
    if (units == "ns") {
        lifetime *= UNITS::nanosecond;
        error *= UNITS::nanosecond;
    }
    if (units == "ps") {
        lifetime *= UNITS::picosecond;
        error *= UNITS::picosecond;
    }
    if (units == "fs") {
        lifetime *= UNITS::femtosecond;
        error *= UNITS::femtosecond;
    }
    return std::make_pair(lifetime, error);
}

std::vector<std::pair<double, double>> LevelScheme::GetMultipleEnergiesandError(const std::string & st) {
    std::vector<std::pair<double, double>> data;
    if (st.empty())
        return data;
    std::string token;
    std::istringstream ss(st);
    while (std::getline(ss, token, ',')) {
        while (     !token.empty()
                    && (token.front() == ' '
                        ||  token.front() == '\240'
                        ||  token.front() == '\302'
                        ||  token.front() == '<')){
            token.erase(0, 1);
        }
        if (!token.empty() && token.back()=='?')
            token.erase(token.end()-1);
        if (token.empty()) {
            data.emplace_back(0, 0);
            continue;
        }
        std::size_t found = token.find_first_of(' ');
        if (found != std::string::npos)
            data.emplace_back(  std::stod(token.substr(0, found))*UNITS::keV,
                                std::stod(token.substr(found, st.size()))*UNITS::keV);
        else
            data.emplace_back(  std::stod(token)*UNITS::keV,
                                0);
    }
    return data;
}

std::vector<std::pair<double, double>> LevelScheme::GetMultipleIntensitiesandError(const std::string & st) {
    std::vector<std::pair<double, double>> data;
    if (st.empty()) {
        data.emplace_back(0, 0);
        return data;
    }
    std::string token;
    std::istringstream ss(st);
    while (std::getline(ss, token, ',')) {
        while (     !token.empty()
                    && (token.front() == ' '
                        ||  token.front() == '\240'
                        ||  token.front() == '\302'
                        ||  token.front() == '\342'
                        ||  token.front() == '\211'
                        ||  token.front() == '\244'
                        ||  token.front() == '<')){
            token.erase(0, 1);
        }
        if (!token.empty() && token.back()=='?')
            token.erase(token.end()-1);
        if (token.empty()) {
            data.emplace_back(0, 0);
            continue;
        }
        std::size_t found = token.find_first_of(' ');
        if (found != std::string::npos)
            data.emplace_back(  std::stod(token.substr(0, found)),
                                std::stod(token.substr(found, st.size())));
        else if(token=="WEAK"){
            data.emplace_back(0, 0);
        }
        else
            data.emplace_back(  std::stod(token),
                                0);
    }
    if (st.back()==',')
        data.emplace_back(0, 0);
    return data;
}

std::vector<std::string> LevelScheme::GetMultipolariy(const std::string & st) {
    std::vector<std::string> data;
    if (st.empty()){
        data.emplace_back("");
        return data;
    }
    std::string token;
    std::istringstream ss(st);
    while (std::getline(ss, token, ',')) {
        while (     !token.empty()
                    && (token.front() == ' '
                        ||  token.front() == '\240'
                        ||  token.front() == '\302'
                        ||  token.front() == '<')){
            token.erase(0, 1);
        }
        if (!token.empty() && token.back()=='?')
            token.erase(token.end()-1);
        if (token.empty()) {
            data.emplace_back(0, 0);
            continue;
        }
        if(token == "D")
            continue;
        data.emplace_back(token);
    }
    if (st.back()==',')
        data.emplace_back("");
    for (long unsigned ii=0; ii<data.size()-1; ++ii){
        if (data[ii].empty()) continue;
        if (    data[ii].front()=='[' && data[ii+1].back()==']'
                &&  data[ii].back()!=']'                            ) {
            data[ii] = data[ii] + "&" + data[ii + 1];
            data.erase(data.begin() +ii + 1);
        }
    }
    return data;
}

std::vector<std::pair<double, Level::Parity>> LevelScheme::GetSpinParity(const std::string & st) {

    std::vector<std::pair<double, Level::Parity>> found_parities;
    if (st.empty()) {
        //std::vector<std::pair<double, Level::Parity>> dummy {(0, Level::Parity::UNKNOWN)};
        //return dummy;
        return found_parities;
    }
    std::string tmp_str= st;
    //Remove unwanted characters
    std::size_t found = tmp_str.find_first_of("()");
    while (found!=std::string::npos) {
        tmp_str.erase(found, 1);
        found = tmp_str.find_first_of("()");
    }
    //Find multiple spins and parities
    //std::vector<std::pair<double, Level::Parity>> found_parities;
    std::string token;
    std::istringstream ss(tmp_str);
    while (std::getline(ss, token, ','))  {
        std::size_t par_token = token.find_first_of("+-");
        if (par_token == std::string::npos)
        {
            double spin = std::stod(token);
            if (token.find("/2")!=std::string::npos)
                spin /=2.;
            found_parities.emplace_back(spin,
                                        Level::Parity::UNKNOWN );
        }
        else{
            if (token[par_token] == '+'){
                found_parities.emplace_back(std::stod(token.erase(par_token, 1)),
                                            Level::Parity::POSITIVE_PARITY );

            }else if( token[par_token] == '-'){
                found_parities.emplace_back( std::stod(token.erase(par_token, 1)),
                                             Level::Parity::NEGATIVE_PARITY );

            } else throw std::runtime_error("Something wrong in sp-par reading\n");
        };
    };
    return found_parities;
}

//std::pair<Gamma::Multipolarity, int> LevelScheme::GetMultipolarity(const std::string &) {
//    return std::make_pair();
//}

std::vector<std::string> LevelScheme::ReadCSVRow(const std::string &row) {
    enum CSVState {
        UnquotedField,
        QuotedField,
        QuotedQuote
    };
    CSVState state = CSVState::UnquotedField;
    std::vector<std::string> fields {""};
    size_t i = 0; // index of the current field
    for (char c : row) {
        switch (state) {
            case CSVState::UnquotedField:
                switch (c) {
                    case ',': // end of field
                        fields.emplace_back(""); i++;
                        break;
                    case '"': state = CSVState::QuotedField;
                        break;
                    default:  fields[i].push_back(c);
                        break; }
                break;
            case CSVState::QuotedField:
                switch (c) {
                    case '"': state = CSVState::QuotedQuote;
                        break;
                    default:  fields[i].push_back(c);
                        break; }
                break;
            case CSVState::QuotedQuote:
                switch (c) {
                    case ',': // , after closing quote
                        fields.emplace_back(""); i++;
                        state = CSVState::UnquotedField;
                        break;
                    case '"': // "" -> "
                        fields[i].push_back('"');
                        state = CSVState::QuotedField;
                        break;
                    default:  // end of quote
                        state = CSVState::UnquotedField;
                        break; }
                break;
        }
    }
    return fields;
}

std::vector<std::pair<double, double>> LevelScheme::ComputeBranching(const std::vector<std::pair<double, double>> & intensities) {
    std::vector<std::pair<double, double>> branchings;
    double itot = 0;
    for (const auto& it: intensities){
        itot +=it.first;
    }
    branchings.reserve(intensities.size());
for (const auto& it: intensities){
        if (itot!=0)
            branchings.emplace_back(it.first/itot, it.second/itot);
        else
            branchings.emplace_back(0, 0);

    }
    return branchings;
}





