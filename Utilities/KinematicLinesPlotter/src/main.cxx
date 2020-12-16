#include "Plotter.h"

#include <iostream>

#include <TApplication.h>

int main(){
    TApplication theApp("app", new int(0), new char*);

    Plotter plotter;
    //plotter.plotVamosAcceptance("graph");
    plotter.plotMugastAcceptance("MM4.root");
    return 0;
}
