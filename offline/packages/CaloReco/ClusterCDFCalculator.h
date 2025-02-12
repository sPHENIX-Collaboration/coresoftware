#ifndef CALORECO_CLUSTERCDFCALCULATOR_H
#define CALORECO_CLUSTERCDFCALCULATOR_H

#include <TFile.h>
#include <TMatrixD.h>

class TH1;

class ClusterCDFCalculator {
public:
    ClusterCDFCalculator() = default;
    virtual ~ClusterCDFCalculator();

    void LoadProfile(const std::string &filename);
    std::pair<double, double> GetCDF(const std::vector<double>& energies, double clusterenergy, int NMATRIXDIM);

private:
    TFile* file{nullptr};
    TH1* henbins_photon{nullptr};
    TH1* henbins_pi0{nullptr};

    std::vector<TH1*> hD2_3x3_photon;
    std::vector<TH1*> hD2_5x5_photon;
    std::vector<TH1*> hD2_7x7_photon;
    std::vector<TH1*> hD2_3x3_pi0;
    std::vector<TH1*> hD2_5x5_pi0;
    std::vector<TH1*> hD2_7x7_pi0;

    std::vector<TMatrixD> inverseCovarianceMatrix_3x3_photon;
    std::vector<TMatrixD> inverseCovarianceMatrix_5x5_photon;
    std::vector<TMatrixD> inverseCovarianceMatrix_7x7_photon;
    std::vector<TMatrixD> inverseCovarianceMatrix_3x3_pi0;
    std::vector<TMatrixD> inverseCovarianceMatrix_5x5_pi0;
    std::vector<TMatrixD> inverseCovarianceMatrix_7x7_pi0;

    std::vector<double> meanD2_3x3_photon;
    std::vector<double> meanD2_5x5_photon;
    std::vector<double> meanD2_7x7_photon;
    std::vector<double> meanD2_3x3_pi0;
    std::vector<double> meanD2_5x5_pi0;
    std::vector<double> meanD2_7x7_pi0;

    std::vector<std::vector<TH1*>> ratioHistograms_3x3_photon;
    std::vector<std::vector<TH1*>> ratioHistograms_5x5_photon;
    std::vector<std::vector<TH1*>> ratioHistograms_7x7_photon;
    std::vector<std::vector<TH1*>> ratioHistograms_3x3_pi0;
    std::vector<std::vector<TH1*>> ratioHistograms_5x5_pi0;
    std::vector<std::vector<TH1*>> ratioHistograms_7x7_pi0;

    void LoadHistogramsAndMatrices();

    double CalculateJointDeviation(const std::vector<double>& energies,
                                    const std::vector<TH1*>& ratioHistograms,
                                    const TMatrixD& inverseCovarianceMatrix,
                                    double meanD2);

    virtual std::unordered_set<int> getSubsetIndices(int gridSize, int subsetSize);
    double GetCDFValue(double jointDeviation, TH1* hD2);

    static const int NMATRIX_3x3 {9};
    static const int NMATRIX_5x5 {25};
    static const int NMATRIX_7x7 {49};
    static const int FULLGRIDSIZE {7};

    static const int NMATRIXDIM3  {3};
    static const int NMATRIXDIM5  {5};
    static const int NMATRIXDIM7  {7};

    const int gridSize {49};
    const double epsilon {1e-9};
};

#endif
