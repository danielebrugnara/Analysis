#include <fstream>
#include <iostream>
#include <regex>
#include <vector>
#include <string>

#include <TH1D.h>

enum DataType{
    A,
    AB,
    AX,
    C,
    UC,
    S,
    US,
    I,
    UI,
    L,
    UL,
    F,
    D,
    UA,
    UT
};

std::pair<std::vector<unsigned int>, std::string> fileNameParser(std::string fileName);
enum DataType getType(const std::string& type);
TH1D* getHistogram(const std::string& fileName, std::vector<unsigned int> indexes);

template<typename T>
TH1D* readData(const std::string& fileName, std::vector<unsigned int> dimensions, std::vector<unsigned int> indexes);

TH1D* tktToRoot(const std::string& fileName = "Post__5-40-16384-UI__Ener.spec",
               std::vector<unsigned int> indexes= std::vector<unsigned int>{0, 0},
               unsigned int nSpec= 1,
               double xcal=0.25){
    indexes.push_back(nSpec);

    TH1D* histo = getHistogram(fileName, indexes);
    histo->GetXaxis()->SetLimits(0, histo->GetXaxis()->GetXmax()*xcal);
    return histo;

    //histo->Draw();
}



//Functions///////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename T>
TH1D* readData(const std::string& fileName,
              std::vector<unsigned int> dimensions,
              std::vector<unsigned int> indexes){

    unsigned int size{dimensions.back()};
    dimensions.pop_back();

    unsigned int nspec{indexes.back()};
    indexes.pop_back();

    T data[size*nspec];

    std::ifstream spectra(fileName, std::ios::in|std::ios::binary);
    if (!spectra.is_open()) throw std::runtime_error("Unable to open file "+fileName+"\n");

    //std::reverse(dimensions.begin(), dimensions.end());

    unsigned int pos{0};
    for(unsigned int i{0}; i<dimensions.size(); ++i){
        int tmp = indexes[i];
        std::cout << "idx"<<i<<": " << indexes[i] <<std::endl;
        for(unsigned int j{(unsigned int)dimensions.size()-1}; j>i; --j) {
            std::cout << "dimension"<<j<<" " << dimensions[j] << std::endl;
            tmp *= dimensions[j];
        }
        pos += tmp;
    }
    pos *=sizeof(T)*size;

    std::cout << "pos: " <<  pos << std::endl;


    spectra.seekg(pos);
    spectra.read(reinterpret_cast<char*>( &data), (sizeof data));

    TH1D* spec = new TH1D("spec", "spec", size*nspec, 0, size*nspec);
    for (int i{0}; i<size*nspec; ++i){
        spec->SetBinContent(spec->FindBin(i), data[i]);
    }
    return spec;
}

std::pair<std::vector<unsigned int>, std::string> fileNameParser(std::string fileName){
    fileName = fileName.substr(0, fileName.find_last_of("_")-1);
    fileName = fileName.substr(fileName.find_first_of("_")+2);
    //std::cout << fileName << std::endl;

    auto const re = std::regex{R"(-)"};
    auto const vec = std::vector<std::string>(
            std::sregex_token_iterator{begin(fileName), end(fileName), re, -1},
            std::sregex_token_iterator{}
    );

    std::pair<std::vector<unsigned int>, std::string> output;
    for (const auto&it: vec){
        if (it != vec.back())
            output.first.push_back(std::stoi(it));
        else
            output.second = it;
    }
    return output;
}

TH1D* getHistogram(const std::string& fileName, std::vector<unsigned int> indexes){
    auto dimensions = fileNameParser(fileName);
    TH1D* histogram{nullptr};
    //std::cout << "type" << getType(dimensions.second) << std::endl;
    switch(getType(dimensions.second)){
                case A:
                    break;
                case AB:
                    break;
                case AX:
                    break;
                case C:
                    break;
                case UC:
                    break;
                case S:
                    histogram = readData<short>(fileName, dimensions.first, indexes);
                    break;
                case US:
                    histogram = readData<unsigned short>(fileName, dimensions.first, indexes);
                    break;
                case I:
                    histogram = readData<int32_t>(fileName, dimensions.first, indexes);
                    break;
                case UI:
                    histogram = readData<uint32_t>(fileName, dimensions.first, indexes);
                    break;
                case L:
                    histogram = readData<long>(fileName, dimensions.first, indexes);
                    break;
                case UL:
                    histogram = readData<unsigned long>(fileName, dimensions.first, indexes);
                    break;
                case F:
                    histogram = readData<float>(fileName, dimensions.first, indexes);
                    break;
                case D:
                    histogram = readData<double>(fileName, dimensions.first, indexes);
                    break;
                case UA:
                    break;
                case UT:
                    break;
                default:
                    break;
    }
    if (histogram != nullptr) return histogram;
    else throw std::runtime_error("shomething wrong\n");

}


enum DataType getType(const std::string& type){
    if (type == "A") return A;
    if (type == "AB") return AB;
    if (type == "AX") return AX;
    if (type == "C") return C;
    if (type == "UC") return UC;
    if (type == "S") return S;
    if (type == "US") return US;
    if (type == "I") return I;
    if (type == "UI") return UI;
    if (type == "L") return L;
    if (type == "UL") return UL;
    if (type == "F") return F;
    if (type == "D") return D;
    if (type == "UA") return UA;
    if (type == "UT") return UT;
    else throw std::runtime_error("Unknown data type");
}
