#pragma once

#include <vector>

#include <TLorentzVector.h>
#include <TVector2.h>
#include <TVector3.h>

class VamosData: public TObject{
    public:
    double En;         //Energy IC
    double D_En;       //DE 1+2
    double D_En2;      //DE 1
    double Path;       //Reconstructed path
    TVector3 FocalPlanePosition;
    TVector3 PTPosition;
    TVector3 EmissionVersor;
    double T;          //Time
    double V;          //Velocity
    double Beta;       //Beta
    double Gamma;      //Gamma
    double M_Q;        //Mass over charge
    double M;          //Mass
    double Charge;     //Charge state
    double BRho;       //Brho
    bool Identified;   //Positive identification
    TLorentzVector p4; //4 momentum of recoil
    TLorentzVector p4MidTarget; //4 momentum of recoil
    int id_M;
    int id_Z;
    int id_Q;
    VamosData() : En(0),
                 D_En(0),
                 D_En2(0),
                 Path(0),
                 FocalPlanePosition(),
                 PTPosition(),
                 EmissionVersor(),
                 T(0),
                 V(0),
                 Beta(0),
                 Gamma(0),
                 M_Q(0),
                 M(0),
                 Charge(0),
                 BRho(0),
                 Identified(false),
                 p4(),
                 p4MidTarget(),
                 id_M(0),
                 id_Z(0),
                 id_Q(0){};

    ClassDef(VamosData, 3);
}; 
