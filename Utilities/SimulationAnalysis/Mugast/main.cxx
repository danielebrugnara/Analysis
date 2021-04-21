#include <iostream>
#include <string>

#include "RunSelector.h"

int main(int argc, char* argv[]){
    std::vector<std::string> files{
        "./../../../DataAnalyzed/simu/46Ar3Hed47K_0keV_s12_0_0.root",
        "./../../../DataAnalyzed/simu/46Ar3Hed47K_360keV_d32_0_0.root",
        "./../../../DataAnalyzed/simu/46Ar3Hed47K_2020keV_f72_0_0.root",
    };
    RunSelector run("../../../DataAnalyzed/simu/sumsimu.root", "../../../DataAnalyzed/sum.root", files);
    //RunSelector run("Data/anaout/anular.root", "../../../DataAnalyzed/sum.root");
    return 0;
}
