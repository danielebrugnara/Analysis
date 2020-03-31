#include "AgataProcessing.h"

AgataProcessing::AgataProcessing() : ref_ts(179),
                                     z_shift(-4),
                                     gammaray(nullptr),
                                     data(nullptr)
{
}

AgataProcessing::~AgataProcessing()
{
    delete gammaray;
    delete data;
}


inline double AgataProcessing::CorrectDoppler(int ii, double additional_shift)
{
    TVector3 tmp_three_vec = gammaray->Pos[ii];
    TLorentzVector tmp_four_vec = gammaray->Pgamma[ii];
    tmp_three_vec.SetZ(gammaray->Pos[ii].Z() + additional_shift);
    tmp_four_vec.SetPhi(tmp_three_vec.Phi());
    tmp_four_vec.SetTheta(tmp_three_vec.Theta());
    tmp_four_vec.SetE(gammaray->E[ii]);
    tmp_four_vec.Boost(-data->p4->BoostVector());

    return tmp_four_vec.Energy();
}