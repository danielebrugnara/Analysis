#include <iostream>

#include "RunSelector.h"

int main(int argc, char* argv[]){
    std::cout << "Starting analysis now \n";
    if (argc < 2)    throw std::runtime_error("set inputs correctly\n");
    for (int ii=0; ii<argc-1;++ii){
        std::cout << "Running selector for : "<<argv[ii+1]<<" \n";
        RunSelector run(argv[ii+1]);
    }
    return 0;
}
