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

void ClusterCDFCalculator::LoadProfile(const std::string &filename)
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
  file->GetObject("h_photon_hen", henbins_photon);
  if (!henbins_photon)
  {
    std::cout << "ClusterCDFCalculator::LoadHistogramsAndMatrices() FATAL ERROR: photon energy bins histogram does not exist.. return.." << std::endl;
    return;
  }

  file->GetObject("h_pi0_hen", henbins_pi0);
  if (!henbins_pi0)
  {
    std::cout << "ClusterCDFCalculator::LoadHistogramsAndMatrices() FATAL ERROR: pi0 energy bins histogram does not exist.. return.." << std::endl;
    return;
  }

  if (henbins_photon->GetNbinsX() != henbins_photon->GetNbinsX())
  {
    std::cout << "ClusterCDFCalculator::LoadHistogramsAndMatrices() FATAL ERROR: inconsistent energy bins for photon and pi0... return.." << std::endl;
    return;
  }

  const int nBins = henbins_photon->GetNbinsX();

  // Resize photon
  hD2_3x3_photon.resize(nBins, nullptr);
  hD2_5x5_photon.resize(nBins, nullptr);
  hD2_7x7_photon.resize(nBins, nullptr);
  meanD2_3x3_photon.resize(nBins, -1);
  meanD2_5x5_photon.resize(nBins, -1);
  meanD2_7x7_photon.resize(nBins, -1);
  inverseCovarianceMatrix_3x3_photon.resize(nBins, TMatrixD(NMATRIX_3x3, NMATRIX_3x3));
  inverseCovarianceMatrix_5x5_photon.resize(nBins, TMatrixD(NMATRIX_5x5, NMATRIX_5x5));
  inverseCovarianceMatrix_7x7_photon.resize(nBins, TMatrixD(NMATRIX_7x7, NMATRIX_7x7));
  ratioHistograms_3x3_photon.resize(nBins, std::vector<TH1 *>(NMATRIX_3x3, nullptr));
  ratioHistograms_5x5_photon.resize(nBins, std::vector<TH1 *>(NMATRIX_5x5, nullptr));
  ratioHistograms_7x7_photon.resize(nBins, std::vector<TH1 *>(NMATRIX_7x7, nullptr));

  // Resize pi0
  hD2_3x3_pi0.resize(nBins, nullptr);
  hD2_5x5_pi0.resize(nBins, nullptr);
  hD2_7x7_pi0.resize(nBins, nullptr);
  meanD2_3x3_pi0.resize(nBins, -1);
  meanD2_5x5_pi0.resize(nBins, -1);
  meanD2_7x7_pi0.resize(nBins, -1);
  inverseCovarianceMatrix_3x3_pi0.resize(nBins, TMatrixD(NMATRIX_3x3, NMATRIX_3x3));
  inverseCovarianceMatrix_5x5_pi0.resize(nBins, TMatrixD(NMATRIX_5x5, NMATRIX_5x5));
  inverseCovarianceMatrix_7x7_pi0.resize(nBins, TMatrixD(NMATRIX_7x7, NMATRIX_7x7));
  ratioHistograms_3x3_pi0.resize(nBins, std::vector<TH1 *>(NMATRIX_3x3, nullptr));
  ratioHistograms_5x5_pi0.resize(nBins, std::vector<TH1 *>(NMATRIX_5x5, nullptr));
  ratioHistograms_7x7_pi0.resize(nBins, std::vector<TH1 *>(NMATRIX_7x7, nullptr));

  for (int binidx = 0; binidx < nBins; ++binidx)
  {
    // photon hist
    std::string hD2_3x3_name_photon = (boost::format("h_photon_hD2_3x3_en%d") % binidx).str();
    std::string hD2_5x5_name_photon = (boost::format("h_photon_hD2_5x5_en%d") % binidx).str();
    std::string hD2_7x7_name_photon = (boost::format("h_photon_hD2_7x7_en%d") % binidx).str();

    file->GetObject(hD2_3x3_name_photon.c_str(), hD2_3x3_photon[binidx]);
    file->GetObject(hD2_5x5_name_photon.c_str(), hD2_5x5_photon[binidx]);
    file->GetObject(hD2_7x7_name_photon.c_str(), hD2_7x7_photon[binidx]);

    if (!hD2_3x3_photon[binidx] || !hD2_5x5_photon[binidx] || !hD2_7x7_photon[binidx])
    {
      std::cout << "ClusterCDFCalculator::LoadHistogramsAndMatrices() error: One or more photon hD2 histograms for bin " << binidx << " are missing. Returning..." << std::endl;
      return;
    }

    TH1 *hD2mean3_photon{nullptr};
    file->GetObject((boost::format("h_photon_hD2mean3_en%d") % binidx).str().c_str(), hD2mean3_photon);
    TH1 *hD2mean5_photon{nullptr};
    file->GetObject((boost::format("h_photon_hD2mean5_en%d") % binidx).str().c_str(), hD2mean5_photon);
    TH1 *hD2mean7_photon{nullptr};
    file->GetObject((boost::format("h_photon_hD2mean7_en%d") % binidx).str().c_str(), hD2mean7_photon);

    if (!hD2mean3_photon || !hD2mean5_photon || !hD2mean7_photon)
    {
      std::cout << "ClusterCDFCalculator::LoadHistogramsAndMatrices() error: One or more required photon mean histograms are missing for bin " << binidx << std::endl;
      return;
    }
    meanD2_3x3_photon[binidx] = hD2mean3_photon->GetBinContent(1);
    meanD2_5x5_photon[binidx] = hD2mean5_photon->GetBinContent(1);
    meanD2_7x7_photon[binidx] = hD2mean7_photon->GetBinContent(1);

    for (int i = 0; i < NMATRIX_3x3; ++i)
    {
      std::string histName = (boost::format("h_photon_heratio_3x3_en%d_%d") % binidx % i).str();
      file->GetObject(histName.c_str(), ratioHistograms_3x3_photon[binidx][i]);
      if (!ratioHistograms_3x3_photon[binidx][i])
      {
        std::cout << "ClusterCDFCalculator::LoadHistogramsAndMatrices() error: photon hist " << histName.c_str() << " is missing." << std::endl;
        return;
      }
    }
    for (int i = 0; i < NMATRIX_5x5; ++i)
    {
      std::string histName = (boost::format("h_photon_heratio_5x5_en%d_%d") % binidx % i).str();
      file->GetObject(histName.c_str(), ratioHistograms_5x5_photon[binidx][i]);
      if (!ratioHistograms_5x5_photon[binidx][i])
      {
        std::cout << "ClusterCDFCalculator::LoadHistogramsAndMatrices() error: photon hist " << histName.c_str() << " is missing." << std::endl;
        return;
      }
    }
    for (int i = 0; i < NMATRIX_7x7; ++i)
    {
      std::string histName = (boost::format("h_photon_heratio_7x7_en%d_%d") % binidx % i).str();
      file->GetObject(histName.c_str(), ratioHistograms_7x7_photon[binidx][i]);
      if (!ratioHistograms_7x7_photon[binidx][i])
      {
        std::cout << "ClusterCDFCalculator::LoadHistogramsAndMatrices() error: photon hist " << histName.c_str() << " is missing." << std::endl;
        return;
      }
    }
    TH2 *hCovMatrix3x3_photon{nullptr};
    file->GetObject((boost::format("h_photon_hCovMatrix3_en%d") % binidx).str().c_str(), hCovMatrix3x3_photon);
    TH2 *hCovMatrix5x5_photon{nullptr};
    file->GetObject((boost::format("h_photon_hCovMatrix5_en%d") % binidx).str().c_str(), hCovMatrix5x5_photon);
    TH2 *hCovMatrix7x7_photon{nullptr};
    file->GetObject((boost::format("h_photon_hCovMatrix7_en%d") % binidx).str().c_str(), hCovMatrix7x7_photon);

    if (!hCovMatrix3x3_photon || !hCovMatrix5x5_photon || !hCovMatrix7x7_photon)
    {
      std::cout << "ClusterCDFCalculator::LoadHistogramsAndMatrices() error: photon covariance matrix histograms are missing for bin " << binidx << std::endl;
      return;
    }
    for (int i = 0; i < NMATRIX_3x3; ++i)
    {
      for (int j = 0; j < NMATRIX_3x3; ++j)
      {
        inverseCovarianceMatrix_3x3_photon[binidx](i, j) = hCovMatrix3x3_photon->GetBinContent(i + 1, j + 1);
      }
    }
    for (int i = 0; i < NMATRIX_5x5; ++i)
    {
      for (int j = 0; j < NMATRIX_5x5; ++j)
      {
        inverseCovarianceMatrix_5x5_photon[binidx](i, j) = hCovMatrix5x5_photon->GetBinContent(i + 1, j + 1);
      }
    }
    for (int i = 0; i < NMATRIX_7x7; ++i)
    {
      for (int j = 0; j < NMATRIX_7x7; ++j)
      {
        inverseCovarianceMatrix_7x7_photon[binidx](i, j) = hCovMatrix7x7_photon->GetBinContent(i + 1, j + 1);
      }
    }
    try
    {
      inverseCovarianceMatrix_3x3_photon[binidx].Invert();
      inverseCovarianceMatrix_5x5_photon[binidx].Invert();
      inverseCovarianceMatrix_7x7_photon[binidx].Invert();
    }
    catch (const std::exception &e)
    {
      std::cout << "ClusterCDFCalculator::LoadHistogramsAndMatrices - Exception in photon covariance invert for " << e.what() << " abort.." << std::endl;
      return;
    }

    // pi0 hist
    std::string hD2_3x3_name_pi0 = (boost::format("h_pi0_hD2_3x3_en%d") % binidx).str();
    std::string hD2_5x5_name_pi0 = (boost::format("h_pi0_hD2_5x5_en%d") % binidx).str();
    std::string hD2_7x7_name_pi0 = (boost::format("h_pi0_hD2_7x7_en%d") % binidx).str();

    file->GetObject(hD2_3x3_name_pi0.c_str(), hD2_3x3_pi0[binidx]);
    file->GetObject(hD2_5x5_name_pi0.c_str(), hD2_5x5_pi0[binidx]);
    file->GetObject(hD2_7x7_name_pi0.c_str(), hD2_7x7_pi0[binidx]);

    if (!hD2_3x3_pi0[binidx] || !hD2_5x5_pi0[binidx] || !hD2_7x7_pi0[binidx])
    {
      std::cout << "ClusterCDFCalculator::LoadHistogramsAndMatrices() error: One or more pi0 hD2 histograms for bin " << binidx << " are missing. Returning..." << std::endl;
      return;
    }

    TH1 *hD2mean3_pi0{nullptr};
    file->GetObject((boost::format("h_pi0_hD2mean3_en%d") % binidx).str().c_str(), hD2mean3_pi0);
    TH1 *hD2mean5_pi0{nullptr};
    file->GetObject((boost::format("h_pi0_hD2mean5_en%d") % binidx).str().c_str(), hD2mean5_pi0);
    TH1 *hD2mean7_pi0{nullptr};
    file->GetObject((boost::format("h_pi0_hD2mean7_en%d") % binidx).str().c_str(), hD2mean7_pi0);

    if (!hD2mean3_pi0 || !hD2mean5_pi0 || !hD2mean7_pi0)
    {
      std::cout << "ClusterCDFCalculator::LoadHistogramsAndMatrices() error: One or more required pi0 mean histograms are missing for bin " << binidx << std::endl;
      return;
    }
    meanD2_3x3_pi0[binidx] = hD2mean3_pi0->GetBinContent(1);
    meanD2_5x5_pi0[binidx] = hD2mean5_pi0->GetBinContent(1);
    meanD2_7x7_pi0[binidx] = hD2mean7_pi0->GetBinContent(1);

    for (int i = 0; i < NMATRIX_3x3; ++i)
    {
      std::string histName = (boost::format("h_pi0_heratio_3x3_en%d_%d") % binidx % i).str();
      file->GetObject(histName.c_str(), ratioHistograms_3x3_pi0[binidx][i]);
      if (!ratioHistograms_3x3_pi0[binidx][i])
      {
        std::cout << "ClusterCDFCalculator::LoadHistogramsAndMatrices() error: pi0 hist " << histName.c_str() << " is missing." << std::endl;
        return;
      }
    }
    for (int i = 0; i < NMATRIX_5x5; ++i)
    {
      std::string histName = (boost::format("h_pi0_heratio_5x5_en%d_%d") % binidx % i).str();
      file->GetObject(histName.c_str(), ratioHistograms_5x5_pi0[binidx][i]);
      if (!ratioHistograms_5x5_pi0[binidx][i])
      {
        std::cout << "ClusterCDFCalculator::LoadHistogramsAndMatrices() error: pi0 hist " << histName.c_str() << " is missing." << std::endl;
        return;
      }
    }
    for (int i = 0; i < NMATRIX_7x7; ++i)
    {
      std::string histName = (boost::format("h_pi0_heratio_7x7_en%d_%d") % binidx % i).str();
      file->GetObject(histName.c_str(), ratioHistograms_7x7_pi0[binidx][i]);
      if (!ratioHistograms_7x7_pi0[binidx][i])
      {
        std::cout << "ClusterCDFCalculator::LoadHistogramsAndMatrices() error: pi0 hist " << histName.c_str() << " is missing." << std::endl;
        return;
      }
    }
    TH2 *hCovMatrix3x3_pi0{nullptr};
    file->GetObject((boost::format("h_pi0_hCovMatrix3_en%d") % binidx).str().c_str(), hCovMatrix3x3_pi0);
    TH2 *hCovMatrix5x5_pi0{nullptr};
    file->GetObject((boost::format("h_pi0_hCovMatrix5_en%d") % binidx).str().c_str(), hCovMatrix5x5_pi0);
    TH2 *hCovMatrix7x7_pi0{nullptr};
    file->GetObject((boost::format("h_pi0_hCovMatrix7_en%d") % binidx).str().c_str(), hCovMatrix7x7_pi0);

    if (!hCovMatrix3x3_pi0 || !hCovMatrix5x5_pi0 || !hCovMatrix7x7_pi0)
    {
      std::cout << "ClusterCDFCalculator::LoadHistogramsAndMatrices() error: pi0 covariance matrix histograms are missing for bin " << binidx << std::endl;
      return;
    }
    for (int i = 0; i < NMATRIX_3x3; ++i)
    {
      for (int j = 0; j < NMATRIX_3x3; ++j)
      {
        inverseCovarianceMatrix_3x3_pi0[binidx](i, j) = hCovMatrix3x3_pi0->GetBinContent(i + 1, j + 1);
      }
    }
    for (int i = 0; i < NMATRIX_5x5; ++i)
    {
      for (int j = 0; j < NMATRIX_5x5; ++j)
      {
        inverseCovarianceMatrix_5x5_pi0[binidx](i, j) = hCovMatrix5x5_pi0->GetBinContent(i + 1, j + 1);
      }
    }
    for (int i = 0; i < NMATRIX_7x7; ++i)
    {
      for (int j = 0; j < NMATRIX_7x7; ++j)
      {
        inverseCovarianceMatrix_7x7_pi0[binidx](i, j) = hCovMatrix7x7_pi0->GetBinContent(i + 1, j + 1);
      }
    }
    try
    {
      inverseCovarianceMatrix_3x3_pi0[binidx].Invert();
      inverseCovarianceMatrix_5x5_pi0[binidx].Invert();
      inverseCovarianceMatrix_7x7_pi0[binidx].Invert();
    }
    catch (const std::exception &e)
    {
      std::cout << "ClusterCDFCalculator::LoadHistogramsAndMatrices - Exception in pi0 covariance invert for " << e.what() << " abort.." << std::endl;
      return;
    }
  }
}

double ClusterCDFCalculator::CalculateJointDeviation(const std::vector<double> &energies,
                                                     const std::vector<TH1 *> &ratioHistograms,
                                                     const TMatrixD &inverseCovarianceMatrix,
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

double ClusterCDFCalculator::GetCDFValue(double jointDeviation, TH1 *hD2)
{
  int bin = hD2->FindBin(jointDeviation);
  if (bin <= 0 || bin > hD2->GetNbinsX())
  {
    return 0.0;
  }
  return hD2->Integral(bin, hD2->GetNbinsX());
}

std::pair<double, double> ClusterCDFCalculator::GetCDF(const std::vector<double> &energies, double clusterenergy, int NMATRIXDIM)
{
  if (NMATRIXDIM != 3 && NMATRIXDIM != 5 && NMATRIXDIM != 7)
  {
    std::cout << "ClusterCDFCalculator::GetCDF error: Invalid dimension as input! Only 3x3, 5x5, 7x7 are available (use input 3, 5, or 7)." << std::endl;
    return std::make_pair(-1, -1);
  }

  int towersize = energies.size();
  if (towersize != gridSize)
  {
    std::cout << "ClusterCDFCalculator::GetCDF error: Invalid tower energy vector size! It must be of size (49)." << std::endl;
    return std::make_pair(-1, -1);
  }

  int ibin = -1;
  if (henbins_photon)
  {
    ibin = henbins_photon->FindFixBin(clusterenergy) - 1;
    if (ibin < 0)
    {
      if (Verbosity() > 0) 
      {
        std::cout << "ClusterCDFCalculator::GetCDF histogram bin below 0... set prob. both to -1." << std::endl;
      }
      return std::make_pair(-1, -1);
    }
    if (ibin >= henbins_photon->GetNbinsX())
    {
      ibin = henbins_photon->GetNbinsX() - 1;
    }
  }
  else
  {
    std::cout << "ClusterCDFCalculator::GetCDF error: Photon energy bins histogram (henbins_photon) is not loaded!" << std::endl;
    return std::make_pair(-1, -1);
  }

  std::unordered_set<int> indices = getSubsetIndices(FULLGRIDSIZE, NMATRIXDIM);
  std::vector<double> inputenergies;
  for (int is = 0; is < towersize; ++is)
  {
    if (indices.count(is)) //NOLINT(readability-container-contains)
    {
      inputenergies.push_back(energies.at(is));
    }
  }

  double photonCDF = -1;
  double pi0CDF = -1;
  if (NMATRIXDIM == NMATRIXDIM3)
  {
    double jointDeviation_photon = CalculateJointDeviation(inputenergies, ratioHistograms_3x3_photon[ibin], inverseCovarianceMatrix_3x3_photon[ibin], meanD2_3x3_photon[ibin]);
    photonCDF = GetCDFValue(jointDeviation_photon, hD2_3x3_photon[ibin]);

    double jointDeviation_pi0 = CalculateJointDeviation(inputenergies, ratioHistograms_3x3_pi0[ibin], inverseCovarianceMatrix_3x3_pi0[ibin], meanD2_3x3_pi0[ibin]);
    pi0CDF = GetCDFValue(jointDeviation_pi0, hD2_3x3_pi0[ibin]);
  }
  else if (NMATRIXDIM == NMATRIXDIM5)
  {
    double jointDeviation_photon = CalculateJointDeviation(inputenergies, ratioHistograms_5x5_photon[ibin], inverseCovarianceMatrix_5x5_photon[ibin], meanD2_5x5_photon[ibin]);
    photonCDF = GetCDFValue(jointDeviation_photon, hD2_5x5_photon[ibin]);

    double jointDeviation_pi0 = CalculateJointDeviation(inputenergies, ratioHistograms_5x5_pi0[ibin], inverseCovarianceMatrix_5x5_pi0[ibin], meanD2_5x5_pi0[ibin]);
    pi0CDF = GetCDFValue(jointDeviation_pi0, hD2_5x5_pi0[ibin]);
  }
  else if (NMATRIXDIM == NMATRIXDIM7)
  {
    double jointDeviation_photon = CalculateJointDeviation(inputenergies, ratioHistograms_7x7_photon[ibin], inverseCovarianceMatrix_7x7_photon[ibin], meanD2_7x7_photon[ibin]);
    photonCDF = GetCDFValue(jointDeviation_photon, hD2_7x7_photon[ibin]);

    double jointDeviation_pi0 = CalculateJointDeviation(inputenergies, ratioHistograms_7x7_pi0[ibin], inverseCovarianceMatrix_7x7_pi0[ibin], meanD2_7x7_pi0[ibin]);
    pi0CDF = GetCDFValue(jointDeviation_pi0, hD2_7x7_pi0[ibin]);
  }
  else
  {
    std::cout << "ClusterCDFCalculator::GetCDF error: Invalid input size. Supported sizes are 3 (for 3x3), 5 (for 5x5), and 7 (for 7x7)." << std::endl;
    return std::make_pair(-1, -1);
  }

  return std::make_pair(photonCDF, pi0CDF);
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
      int row = (center / _gridsize) + i;
      int col = (center % _gridsize) + j;
      if (row >= 0 && row < _gridsize && col >= 0 && col < _gridsize)
      {
        indices.insert((row * _gridsize) + col);
      }
    }
  }
  return indices;
}
