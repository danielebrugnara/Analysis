#include <iostream>

#include "RunSelector.h"

int main(int argc, char* argv[]){
    RunSelector run("../../../DataAnalyzed/sumsimu.root", "../../../DataAnalyzed/sum.root");
    //RunSelector run("Data/anaout/anular.root", "../../../DataAnalyzed/sum.root");
    return 0;
}
