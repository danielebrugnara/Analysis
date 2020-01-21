#include "MugastIdentification.h"

MugastIdentification::MugastIdentification() : cuts_MG({1, 3, 4, 5, 7, 11}),
                                               cuts_M({1, 2, 4}),
                                               cuts_Z({1, 2}) {}

MugastIdentification::~MugastIdentification() {
    delete gas_thickness;
    delete havar_angle;
    delete data;
    delete fragment;
}

bool MugastIdentification::Initialize() {
    gas_thickness = new Interpolation("./Configs/Interpolations/GasThickness.txt");
    havar_angle = new Interpolation("./Configs/Interpolations/EntranceAngleHavar.txt");

    //Needed masses [M][Z]
    mass[2][1] = 2.01410177812 * AMU_TO_MEV;  //In MeV
    mass[1][1] = 1.00782503223 * AMU_TO_MEV;  //In MeV
}