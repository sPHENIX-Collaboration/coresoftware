#ifndef CLUSTERCDFCALCULATOR_H
#define CLUSTERCDFCALCULATOR_H

#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TMatrixD.h>
#include <TVectorD.h>

#include <boost/format.hpp>

#include <cmath>
#include <iostream>
#include <vector>


class ClusterCDFCalculator
{
 private:
  TFile* file{nullptr};
  TH1D* hD2_3x3{nullptr};
  TH1D* hD2_5x5{nullptr};
  TH1D* hD2_7x7{nullptr};
  TMatrixD inverseCovarianceMatrix_3x3;
  TMatrixD inverseCovarianceMatrix_5x5;
  TMatrixD inverseCovarianceMatrix_7x7;
  double meanD2_3x3{-1};
  double meanD2_5x5{-1};
  double meanD2_7x7{-1};
  std::vector<TH1D*> ratioHistograms_3x3;
  std::vector<TH1D*> ratioHistograms_5x5;
  std::vector<TH1D*> ratioHistograms_7x7;

  void LoadHistogramsAndMatrices();

  double CalculateJointDeviation(const std::vector<double>& energies,
                                 const std::vector<TH1D*>& ratioHistograms,
                                 const TMatrixD& inverseCovarianceMatrix,
                                 double meanD2);

  virtual std::unordered_set<int> getSubsetIndices(int gridSize, int subsetSize);
  double GetCDFValue(double jointDeviation, TH1D* hD2);

  // matrix constant
  static const int NMATRIX_3x3 = 9;
  static const int NMATRIX_5x5 = 25;
  static const int NMATRIX_7x7 = 49;
  static const int FULLGRIDSIZE = 7;

  // dimensions
  static const int NMATRIXDIM3 = 3;
  static const int NMATRIXDIM5 = 5;
  static const int NMATRIXDIM7 = 7;

  // threshold and grid size
  const int gridSize = 49;
  const double epsilon = 1e-9;

 public:
  ClusterCDFCalculator() = default;
  virtual ~ClusterCDFCalculator();

  void LoadProfile(const std::string& filename);
  double GetCDF(const std::vector<double>& energies, int NMATRIX);
};

#endif
