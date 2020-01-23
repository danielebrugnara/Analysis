
#ifndef __MUGASTIDENTIFICATION_H__
#define __MUGASTIDENTIFICATION_H__

#include <TTreeReaderValue.h>
#include <TVector3.h>

#include <array>

#include "Identification.h"
#include "Interpolation.h"
#include "TMugastPhysics.h"

class MugastIdentification : public Identification {
   public:
    MugastIdentification();
    ~MugastIdentification();
    bool Initialize();

    struct Data {
        TTreeReaderValue<TMugastPhysics> *Mugast;
        Data(TTreeReaderValue<TMugastPhysics> *Mugast) : Mugast(Mugast){};
    };
    std::array<int, 6> cuts_MG;
    std::array<int, 3> cuts_M;
    std::array<int, 2> cuts_Z;
    std::array<std::string, 3> cuts_particles;

   private:
    struct Fragment {
        const unsigned int multiplicity;
        std::vector<TVector3> Pos;
        std::vector<double> SI_E;
        std::vector<double> SI_E2;
        std::vector<double> E;
        std::vector<double> E2;
        std::vector<double> T;
        std::vector<double> T2;
        std::vector<double> MG;
        std::vector<int> M;
        std::vector<int> Z;
        std::vector<double> Ex;
        Fragment(const unsigned int multiplicity) : multiplicity(multiplicity) {
            Pos.resize(multiplicity);
            SI_E.resize(multiplicity);
            SI_E2.resize(multiplicity);
            E.resize(multiplicity);
            E2.resize(multiplicity);
            T.resize(multiplicity);
            T2.resize(multiplicity);
            MG.resize(multiplicity);
            M.resize(multiplicity);
            Z.resize(multiplicity);
        };
    };
    Data const *data;
    Fragment *fragment;

    Interpolation *gas_thickness;
    Interpolation *havar_angle;
    const Double_t AMU_TO_MEV{931.4936148};
    std::unordered_map<int, std::unordered_map<int, double>> mass;

    bool with_cuts;

   public:
    inline void SetData(Data const *data) {
        if (this->data != nullptr)
            delete this->data;
        if (this->fragment != nullptr)
            delete this->fragment;
        this->data = data;
        fragment = new Fragment((**(data->Mugast)).DSSD_E.size());
    };

    inline bool Identify() {
#ifdef VERBOSE_DEBUG
        std::cout << "------------>MugastIdentification::Identify()\n";
#endif
        for (unsigned int ii = 0; ii < fragment->multiplicity; ++ii) {
            fragment->Pos[ii] = TVector3((**(data->Mugast)).PosX[ii],
                                         (**(data->Mugast)).PosY[ii],
                                         (**(data->Mugast)).PosZ[ii]);
            fragment->SI_E[ii] = (**(data->Mugast)).DSSD_E[ii];
            fragment->SI_E2[ii] = (**(data->Mugast)).SecondLayer_E[ii];
            fragment->T[ii] = (**(data->Mugast)).DSSD_T[ii];
            fragment->T2[ii] = (**(data->Mugast)).SecondLayer_T[ii];
            fragment->MG[ii] = (**(data->Mugast)).TelescopeNumber[ii];
        }
        try {
            for (unsigned int ii = 0; ii < fragment->multiplicity; ++ii) {
                bool already_id = false;
                for (const auto &cut_it : cuts_particles) {
                    if (with_cuts && cut_type["E_TOF"].at("E_TOF_" + cut_it + "_MG" + std::to_string(fragment->MG[ii]))->IsInside(fragment->SI_E[ii], fragment->T[ii])) {
                        if (already_id) throw std::runtime_error("Overlapping MUGAST E TOF gates");

                        fragment->M[ii] = std::stoi(cut_it.substr(cut_it.find_first_of("m") + 1,
                                                                  cut_it.find_first_of("_") - cut_it.find_first_of("m") - 1));

                        fragment->Z[ii] = std::stoi(cut_it.substr(cut_it.find_first_of("z") + 1,
                                                                  cut_it.find_first_of("_") - cut_it.find_first_of("z") - 1));

                        already_id = true;
                    }
                }
                if (!already_id) {
                    fragment->M[ii] = 0;
                    fragment->Z[ii] = 0;
                }
            };
        } catch (const std::out_of_range &err) {
            std::cerr << "Mugast cuts not found\n";
            with_cuts = false;
        }
        return true;
        //TODO: reconstruct E here

        //TODO: compute Ex here
    }

    inline int Get_Mult() { return fragment->multiplicity; };
    inline TVector3 *Get_Pos(const int &i) { return &(fragment->Pos[i]); };
    inline double Get_E(const int &i) { return fragment->E[i]; };
    inline double Get_SI_E(const int &i) { return fragment->SI_E[i]; };
    inline double Get_T(const int &i) { return fragment->T[i]; };
    inline double Get_T2(const int &i) { return fragment->T[i]; };
    inline double Get_MG(const int &i) { return fragment->MG[i]; };
    inline double Get_M(const int &i) { return fragment->M[i]; };
    inline double Get_Z(const int &i) { return fragment->Z[i]; };
};

#endif