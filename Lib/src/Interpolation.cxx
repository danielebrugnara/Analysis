#ifndef INTERPOLATION_H 
#define INTERPOLATION_H

#include "Interpolation.h"

//ClassImp(Interpolation);

Interpolation::Interpolation(std::string input_file){
    ReadFile(input_file);
}

Interpolation::~Interpolation(){
    delete spline;
}

bool Interpolation::ReadFile(std::string input_file_name){
    std::ifstream in_file(input_file_name);
    if (!in_file.is_open()) return false;
    std::string line;
    std::string delimiter = ", ";
    int ii=-1;
        std::cout <<"Reading file\n";
    while (std::getline(in_file, line)){
        ii++;
        if (ii==0) continue;
        std::string x_string = line.substr(2,line.find(", ")-2);
        std::string y_string = line.substr(line.find(", ")+2, line.rfind("}")-line.find(", ")-3 );
        double x, y;
        try {
            x=std::stof(x_string);
            y=std::stof(y_string);
            if ( ii>1 && fValues_x.back()>x ) break;
            fValues_x.push_back(x);
            fValues_y.push_back(y);
            
         //   std::cout <<"Setting value : " <<fValues_x.back() << " and " <<fValues_y.back() <<std::endl;
            
        }catch(std::invalid_argument & err){
            std::cerr << "Invalid arg :" << err.what()<<"\n";
//            return false;
        }
        spline = new TSpline3(input_file_name.c_str(),&fValues_x[0], &fValues_y[0], fValues_x.size());
    }
  //  this->BuildCoeff();
    return true;
}

TSpline3* Interpolation::GetSpline(){
    return spline;
}

double Interpolation::Evaluate(double x){
    return spline->Eval(x);
}

double Interpolation::GetXPoint(int i){
    return fValues_x.at(i);
}
double Interpolation::GetYPoint(int i){
    return fValues_y.at(i);
}


#endif