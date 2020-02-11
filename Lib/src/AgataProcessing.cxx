#include "AgataProcessing.h"

AgataProcessing::AgataProcessing(unsigned long ref_ts,
                                    const double z_shift)
                                        :ref_ts(ref_ts),
                                        z_shift(z_shift),
                                        gammaray(nullptr),
                                        data(nullptr){
}

AgataProcessing::~AgataProcessing(){
 delete gammaray;
 delete data;    
}

inline void AgataProcessing::Process(){
    for (int ii = 0; ii < **(data->nbAdd); ++ii){                            
        ComputeDoppler(ii);
    }
}

inline void AgataProcessing::ComputeDoppler(int ii){
    gammaray->E[ii] = (*data->AddE)[ii];
    TLorentzVector *tmp_ptr = &(gammaray->Pgamma[ii]);
    gammaray->Pos[ii] = TVector3((*data->AddX)[ii],
                                    (*data->AddY)[ii],
                                    (*data->AddZ)[ii]+z_shift);
    *tmp_ptr = TLorentzVector( gammaray->E[ii],
                                            0,
                                            0,
                                            gammaray->E[ii]);

    tmp_ptr->SetPhi(gammaray->Pos[ii].Phi()); 
    tmp_ptr->SetTheta(gammaray->Pos[ii].Theta());
    tmp_ptr->SetE(gammaray->E[ii]);
    tmp_ptr->Boost(- data->p4->BoostVector());

    gammaray->EDC[ii] = tmp_ptr->Energy();
}

inline double AgataProcessing::CorrectDoppler(int  ii, double additional_shift){
    TVector3 tmp_three_vec = gammaray->Pos[ii];
    TLorentzVector tmp_four_vec = gammaray->Pgamma[ii];
    tmp_three_vec.SetZ(gammaray->Pos[ii].Z()+additional_shift);
    tmp_four_vec.SetPhi(tmp_three_vec.Phi());
    tmp_four_vec.SetTheta(tmp_three_vec.Theta());
    tmp_four_vec.SetE(gammaray->E[ii]);
    tmp_four_vec.Boost(- data->p4->BoostVector());

    return tmp_four_vec.Energy();
}


inline void AgataProcessing::SetData(Data const* data){
    if (this->data != nullptr)
        delete this->data;
    if (this->gammaray != nullptr)
        delete this->gammaray;
    this->data = data;
    gammaray = new Gamma(**(data->nbAdd));
}