// Tell emacs that this is a C++ source
// -*- C++ -*-.
//
#ifndef TRKRTRUTHCENTROIDBULDER_H
#define TRKRTRUTHCENTROIDBULDER_H

// Two simple helper class to build simple centroid hits in the pads of the TPC.
// Used by PHG4TpcElectronDrift.{h,cc}

// TODO:
// See if sPHENIX does 0..2*M_PI or -M_PI..M_PI boundaries
#include <array>

class TrkrTruthCentroidBuilder {
    public:
    TrkrTruthCentroidBuilder(){};
    double sum_phi   {0}, sumSq_phi   {0},
           sum_z     {0}, sumSq_z     {0},
           sumE {0};
    void add_electron(float phi, float z, float e_weight=1.);
    bool is_empty   { true }; // becomes true with any entry
    bool near_phi_boundary { false };
    std::array<float,5> build_centroid(); // {mean.phi,stddev.phi,mean.z,stddev.z,sumE}
    void reset();
    private:
    std::pair<float,float> calc_mean_var(const double& sum, const double& sumSq, const double& weight);
};

#endif
