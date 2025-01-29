#include "ClusterCDFCalculator.h"

#include <TH1.h>
#include <TH2.h>
#include <TMatrixD.h>
#include <TVectorD.h>

#include <boost/format.hpp>

#include <iostream>

ClusterCDFCalculator::~ClusterCDFCalculator()
{
  if (file)
  {
    file->Close();
    delete file;
  }
}

void ClusterCDFCalculator::LoadProfile(const std::string& filename)
{
  if (file)
  {
    std::cout << "ClusterCDFCalculator::LoadProfile Warning: a file is already loaded.. closing it first" << std::endl;
    file->Close();
    delete file;
    file = nullptr;
  }

  file = new TFile(filename.c_str(), "read");

  if (!file || file->IsZombie())
  {
    std::cout << "ClusterCDFCalculator::LoadProfile error: Could not open file!" << std::endl;
    return;
  }

  LoadHistogramsAndMatrices();
}

void ClusterCDFCalculator::LoadHistogramsAndMatrices()
{
  file->GetObject("hen",henbins);
  if (!henbins)
  {
    std::cout << "ClusterCDFCalculator::LoadHistogramsAndMatrices() FATAL ERROR : energy bins histogram does not exist.. return.." << std::endl;
    return;
  }
  const int nBins = henbins->GetNbinsX();
  hD2_3x3.resize(nBins, nullptr);
  hD2_5x5.resize(nBins, nullptr);
  hD2_7x7.resize(nBins, nullptr);

  meanD2_3x3.resize(nBins, -1);
  meanD2_5x5.resize(nBins, -1);
  meanD2_7x7.resize(nBins, -1);

  inverseCovarianceMatrix_3x3.resize(nBins, TMatrixD(NMATRIX_3x3, NMATRIX_3x3));
  inverseCovarianceMatrix_5x5.resize(nBins, TMatrixD(NMATRIX_5x5, NMATRIX_5x5));
  inverseCovarianceMatrix_7x7.resize(nBins, TMatrixD(NMATRIX_7x7, NMATRIX_7x7));

  ratioHistograms_3x3.resize(nBins, std::vector<TH1*>(NMATRIX_3x3, nullptr));
  ratioHistograms_5x5.resize(nBins, std::vector<TH1*>(NMATRIX_5x5, nullptr));
  ratioHistograms_7x7.resize(nBins, std::vector<TH1*>(NMATRIX_7x7, nullptr));

  for (int binidx = 0; binidx < nBins; ++binidx)
  {
    std::string hD2_3x3_name = (boost::format("hD2_3x3_en%d") % binidx).str();
    std::string hD2_5x5_name = (boost::format("hD2_5x5_en%d") % binidx).str();
    std::string hD2_7x7_name = (boost::format("hD2_7x7_en%d") % binidx).str();

    file->GetObject(hD2_3x3_name.c_str(),hD2_3x3[binidx]);
    file->GetObject(hD2_5x5_name.c_str(),hD2_5x5[binidx]);
    file->GetObject(hD2_7x7_name.c_str(),hD2_7x7[binidx]);

    if (!hD2_3x3[binidx] || !hD2_5x5[binidx] || !hD2_7x7[binidx])
    {
      std::cout << "ClusterCDFCalculator::LoadHistogramsAndMatrices() error: One or more hD2 histograms for bin " << binidx << " are missing. Returning..." << std::endl;
      return;
    }
    TH1 *hD2mean3 {nullptr};
    file->GetObject((boost::format("hD2mean3_en%d") % binidx).str().c_str(),hD2mean3);
    TH1 *hD2mean5{nullptr};
    file->GetObject((boost::format("hD2mean5_en%d") % binidx).str().c_str(),hD2mean5);
    TH1 *hD2mean7{nullptr};
    file->GetObject((boost::format("hD2mean7_en%d") % binidx).str().c_str(),hD2mean7);

    if (!hD2mean3 || !hD2mean5 || !hD2mean7)
    {
      std::cout << "ClusterCDFCalculator::LoadHistogramsAndMatrices() error: One or more required histograms (hD2_3x3, hD2_5x5, hD2_7x7) are missing." << std::endl;
      return;
    }

    meanD2_3x3[binidx] = hD2mean3->GetBinContent(1);
    meanD2_5x5[binidx] = hD2mean5->GetBinContent(1);
    meanD2_7x7[binidx] = hD2mean7->GetBinContent(1);

    for (int i = 0; i < NMATRIX_3x3; ++i)
    {
      std::string histName = (boost::format("heratio_3x3_en%d_%d") % binidx % i).str();
      file->GetObject(histName.c_str(),ratioHistograms_3x3[binidx][i]);
      if (!ratioHistograms_3x3[binidx][i])
      {
        std::cout << "ClusterCDFCalculator::LoadHistogramsAndMatrices() error: hist " << histName.c_str() << " is missing." << std::endl;
        return;
      }
    }

    for (int i = 0; i < NMATRIX_5x5; ++i)
    {
      std::string histName = (boost::format("heratio_5x5_en%d_%d") % binidx % i).str();
      file->GetObject(histName.c_str(),ratioHistograms_5x5[binidx][i]);
      if (!ratioHistograms_5x5[binidx][i])
      {
        std::cout << "ClusterCDFCalculator::LoadHistogramsAndMatrices() error: hist " << histName.c_str() << " is missing." << std::endl;
        return;
      }
    }

    for (int i = 0; i < NMATRIX_7x7; ++i)
    {
      std::string histName = (boost::format("heratio_7x7_en%d_%d") % binidx % i).str();
      file->GetObject(histName.c_str(),ratioHistograms_7x7[binidx][i]);
      if (!ratioHistograms_7x7[binidx][i])
      {
        std::cout << "ClusterCDFCalculator::LoadHistogramsAndMatrices() error: hist " << histName.c_str() << " is missing." << std::endl;
        return;
      }
    }
    TH2 *hCovMatrix3x3 {nullptr};
    file->GetObject((boost::format("hCovMatrix3_en%d") % binidx).str().c_str(),hCovMatrix3x3);
    TH2 *hCovMatrix5x5 {nullptr};
    file->GetObject((boost::format("hCovMatrix5_en%d") % binidx).str().c_str(),hCovMatrix5x5);
    TH2 *hCovMatrix7x7 {nullptr};
    file->GetObject((boost::format("hCovMatrix7_en%d") % binidx).str().c_str(),hCovMatrix7x7);

    if (!hCovMatrix3x3 || !hCovMatrix5x5 || !hCovMatrix7x7)
    {
      std::cout << "ClusterCDFCalculator::LoadHistogramsAndMatrices() error: hCovMatrix5x5 histograms are missing." << std::endl;
      return;
    }
    for (int i = 0; i < NMATRIX_3x3; ++i)
    {
      for (int j = 0; j < NMATRIX_3x3; ++j)
      {
        inverseCovarianceMatrix_3x3[binidx](i, j) = hCovMatrix3x3->GetBinContent(i + 1, j + 1);
      }
    }
    for (int i = 0; i < NMATRIX_5x5; ++i)
    {
      for (int j = 0; j < NMATRIX_5x5; ++j)
      {
        inverseCovarianceMatrix_5x5[binidx](i, j) = hCovMatrix5x5->GetBinContent(i + 1, j + 1);
      }
    }
    for (int i = 0; i < NMATRIX_7x7; ++i)
    {
      for (int j = 0; j < NMATRIX_7x7; ++j)
      {
        inverseCovarianceMatrix_7x7[binidx](i, j) = hCovMatrix7x7->GetBinContent(i + 1, j + 1);
      }
    }

    try
    {
      inverseCovarianceMatrix_3x3[binidx].Invert();
      inverseCovarianceMatrix_5x5[binidx].Invert();
      inverseCovarianceMatrix_7x7[binidx].Invert();
    }
    catch (const std::exception& e)
    {
      std::cout << "ClusterCDFCalculator::LoadHistogramsAndMatrices - Exception covariance invert for " << e.what() << " abort.." << std::endl;
      return;
    }
  }
}

double ClusterCDFCalculator::CalculateJointDeviation(const std::vector<double>& energies,
                                                     const std::vector<TH1*>& ratioHistograms,
                                                     const TMatrixD& inverseCovarianceMatrix,
                                                     double meanD2)
{
  int nDims = energies.size();
  double totalEnergy = 0.0;
  for (double energy : energies)
  {
    totalEnergy += energy;
  }

  TVectorD probVector(nDims);
  for (int i = 0; i < nDims; ++i)
  {
    double fraction = energies[i] / totalEnergy;
    int bin = std::min(std::max(ratioHistograms[i]->FindBin(fraction), 1), ratioHistograms[i]->GetNbinsX());
    double probability = std::max(ratioHistograms[i]->GetBinContent(bin), epsilon);
    probVector[i] = -log(probability);
  }

  return probVector * (inverseCovarianceMatrix * probVector) / meanD2;
}

double ClusterCDFCalculator::GetCDFValue(double jointDeviation, TH1* hD2)
{
  int bin = hD2->FindBin(jointDeviation);
  if (bin <= 0 || bin > hD2->GetNbinsX())
  {
    return 0.0;
  }
  return hD2->Integral(bin, hD2->GetNbinsX());
}

double ClusterCDFCalculator::GetCDF(const std::vector<double>& energies, double clusterenergy, int NMATRIXDIM)
{
  if (NMATRIXDIM != 3 && NMATRIXDIM != 5 && NMATRIXDIM != 7)
  {
    std::cout << "ClusterCDFCalculator::GetCDF error: Invalid dimension as input! only 3x3, 5x5, 7x7 are available use input 3 or 5 or 7." << std::endl;
    return -1;
  }

  int towersize = energies.size();
  if (towersize != gridSize)
  {
    std::cout << "ClusterCDFCalculator::GetCDF error: Invalid tower energy vector size! Need to be in size of 49 (7x7)." << std::endl;
    return -1;
  }

  int ibin = -1;
  if (henbins)
  {
    ibin = henbins->FindFixBin(clusterenergy) - 1;
    if (ibin < 0)
    {
      std::cout << "ClusterCDFCalculator::GetCDF fatal error: histogram bin below 0... set prob. to -1." << std::endl;
      return -1;
    }
    else if (ibin >= henbins->GetNbinsX())
    {
      ibin = henbins->GetNbinsX() - 1;
    }
  }
  else
  {
    std::cout << "ClusterCDFCalculator::GetCDF error: Energy bins histogram (henbins) is not loaded!" << std::endl;
    return -1;
  }

  std::unordered_set<int> indices = getSubsetIndices(FULLGRIDSIZE, NMATRIXDIM);
  std::vector<double> inputenergies;
  for (int is = 0; is < towersize; ++is)
  {
    if (indices.count(is))
    {
      double e = energies.at(is);
      inputenergies.push_back(e);
    }
  }

  if (NMATRIXDIM == NMATRIXDIM3)
  {
    double jointDeviation = CalculateJointDeviation(inputenergies, ratioHistograms_3x3[ibin], inverseCovarianceMatrix_3x3[ibin], meanD2_3x3[ibin]);
    return GetCDFValue(jointDeviation, hD2_3x3[ibin]);
  }
  else if (NMATRIXDIM == NMATRIXDIM5)
  {
    double jointDeviation = CalculateJointDeviation(inputenergies, ratioHistograms_5x5[ibin], inverseCovarianceMatrix_5x5[ibin], meanD2_5x5[ibin]);
    return GetCDFValue(jointDeviation, hD2_5x5[ibin]);
  }
  else if (NMATRIXDIM == NMATRIXDIM7)
  {
    double jointDeviation = CalculateJointDeviation(inputenergies, ratioHistograms_7x7[ibin], inverseCovarianceMatrix_7x7[ibin], meanD2_7x7[ibin]);
    return GetCDFValue(jointDeviation, hD2_7x7[ibin]);
  }
  else
  {
    std::cout << "ClusterCDFCalculator::GetCDF error: Invalid input size. Supported sizes are 3 (for 3x3), 5 (for 5x5), and 7 (for 7x7)." << std::endl;
    return -1;
  }
}

std::unordered_set<int> ClusterCDFCalculator::getSubsetIndices(int _gridsize, int _subsetsize)
{
  std::unordered_set<int> indices;
  if (_subsetsize > _gridsize || _subsetsize % 2 == 0)
  {
    std::cout << "Invalid subset size!" << std::endl;
    return indices;
  }

  int center = (_gridsize * _gridsize - 1) / 2;
  int halfSubset = _subsetsize / 2;

  for (int i = -halfSubset; i <= halfSubset; ++i)
  {
    for (int j = -halfSubset; j <= halfSubset; ++j)
    {
      int row = center / _gridsize + i;
      int col = center % _gridsize + j;
      if (row >= 0 && row < _gridsize && col >= 0 && col < _gridsize)
      {
        indices.insert(row * _gridsize + col);
      }
    }
  }
  return indices;
}
