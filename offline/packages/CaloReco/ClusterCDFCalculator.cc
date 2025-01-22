#include "ClusterCDFCalculator.h"

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
    std::cerr << "ClusterCDFCalculator::LoadProfile Warning: a file is already loaded.. closing it first" << std::endl;
    file->Close();
    delete file;
    file = nullptr;
  }

  file = new TFile(filename.c_str(), "read");

  if (!file || file->IsZombie())
  {
    std::cerr << "ClusterCDFCalculator::LoadProfile error: Could not open file!" << std::endl;
    return;
  }

  LoadHistogramsAndMatrices();
}

void ClusterCDFCalculator::LoadHistogramsAndMatrices()
{
  hD2_3x3 = (TH1D*) file->Get("hD2_3x3");
  hD2_5x5 = (TH1D*) file->Get("hD2_5x5");
  hD2_7x7 = (TH1D*) file->Get("hD2_7x7");
  if (!hD2_3x3 || !hD2_5x5 || !hD2_7x7)
  {
    std::cerr << "ClusterCDFCalculator::LoadHistogramsAndMatrices() error: One or more required histograms (hD2_3x3, hD2_5x5, hD2_7x7) are missing." << std::endl;
    return;
  }

  auto hD2mean3 = (TH1D*) file->Get("hD2mean3");
  auto hD2mean5 = (TH1D*) file->Get("hD2mean5");
  auto hD2mean7 = (TH1D*) file->Get("hD2mean7");

  if (!hD2mean3 || !hD2mean5 || !hD2mean7)
  {
    std::cerr << "ClusterCDFCalculator::LoadHistogramsAndMatrices() error: One or more required histograms (hD2_3x3, hD2_5x5, hD2_7x7) are missing." << std::endl;
    return;
  }

  meanD2_3x3 = hD2mean3->GetBinContent(1);
  meanD2_5x5 = hD2mean5->GetBinContent(1);
  meanD2_7x7 = hD2mean7->GetBinContent(1);

  for (int i = 0; i < NMATRIX_3x3; ++i)
  {
    std::string histName = (boost::format("heratio_3x3_%d") % i).str();
    auto hist = (TH1D*) file->Get(histName.c_str());
    if (!hist)
    {
      std::cerr << "ClusterCDFCalculator::LoadHistogramsAndMatrices() error: hist " << histName.c_str() << " is missing." << std::endl;
      return;
    }
    ratioHistograms_3x3.push_back(hist);
  }
  for (int i = 0; i < NMATRIX_5x5; ++i)
  {
    std::string histName = (boost::format("heratio_5x5_%d") % i).str();
    auto hist = (TH1D*) file->Get(histName.c_str());
    if (!hist)
    {
      std::cerr << "ClusterCDFCalculator::LoadHistogramsAndMatrices() error: hist " << histName.c_str() << " is missing." << std::endl;
      return;
    }
    ratioHistograms_5x5.push_back(hist);
  }
  for (int i = 0; i < NMATRIX_7x7; ++i)
  {
    std::string histName = (boost::format("heratio_7x7_%d") % i).str();
    auto hist = (TH1D*) file->Get(histName.c_str());
    if (!hist)
    {
      std::cerr << "ClusterCDFCalculator::LoadHistogramsAndMatrices() error: hist " << histName.c_str() << " is missing." << std::endl;
      return;
    }
    ratioHistograms_7x7.push_back(hist);
  }

  TH2D* hCovMatrix3x3 = (TH2D*) file->Get("hCovMatrix3");
  TH2D* hCovMatrix5x5 = (TH2D*) file->Get("hCovMatrix5");
  TH2D* hCovMatrix7x7 = (TH2D*) file->Get("hCovMatrix7");

  if (!hCovMatrix3x3 || !hCovMatrix5x5 || !hCovMatrix7x7)
  {
    std::cerr << "ClusterCDFCalculator::LoadHistogramsAndMatrices() error: hCovMatrix5x5 histograms are missing." << std::endl;
    return;
  }

  inverseCovarianceMatrix_3x3.ResizeTo(NMATRIX_3x3, NMATRIX_3x3);
  inverseCovarianceMatrix_5x5.ResizeTo(NMATRIX_5x5, NMATRIX_5x5);
  inverseCovarianceMatrix_7x7.ResizeTo(NMATRIX_7x7, NMATRIX_7x7);

  for (int i = 0; i < NMATRIX_3x3; ++i)
  {
    for (int j = 0; j < NMATRIX_3x3; ++j)
    {
      inverseCovarianceMatrix_3x3(i, j) = hCovMatrix3x3->GetBinContent(i + 1, j + 1);
    }
  }
  for (int i = 0; i < NMATRIX_5x5; ++i)
  {
    for (int j = 0; j < NMATRIX_5x5; ++j)
    {
      inverseCovarianceMatrix_5x5(i, j) = hCovMatrix5x5->GetBinContent(i + 1, j + 1);
    }
  }
  for (int i = 0; i < NMATRIX_7x7; ++i)
  {
    for (int j = 0; j < NMATRIX_7x7; ++j)
    {
      inverseCovarianceMatrix_7x7(i, j) = hCovMatrix7x7->GetBinContent(i + 1, j + 1);
    }
  }

  try
  {
    inverseCovarianceMatrix_3x3.Invert();
    inverseCovarianceMatrix_5x5.Invert();
    inverseCovarianceMatrix_7x7.Invert();
  }
  catch (const std::exception& e)
  {
    std::cerr << "ClusterCDFCalculator::LoadHistogramsAndMatrices - Exception covariance invert for " << e.what() << " abort.." << std::endl;
    return;
  }
}

double ClusterCDFCalculator::CalculateJointDeviation(const std::vector<double>& energies,
                                                     const std::vector<TH1D*>& ratioHistograms,
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

double ClusterCDFCalculator::GetCDFValue(double jointDeviation, TH1D* hD2)
{
  int bin = hD2->FindBin(jointDeviation);
  if (bin <= 0 || bin > hD2->GetNbinsX())
  {
    return 0.0;
  }
  return hD2->Integral(bin, hD2->GetNbinsX());
}

double ClusterCDFCalculator::GetCDF(const std::vector<double>& energies, int NMATRIXDIM)
{
  if (NMATRIXDIM != 3 && NMATRIXDIM != 5 && NMATRIXDIM != 7)
  {
    std::cerr << "ClusterCDFCalculator::GetCDF error: Invalid dimension as input! only 3x3, 5x5, 7x7 are available use input 3 or 5 or 7." << std::endl;
    return -1;
  }

  int towersize = energies.size();
  if (towersize != gridSize)
  {
    std::cerr << "ClusterCDFCalculator::GetCDF error: Invalid tower energy vector size! Need to be in size of 49 (7x7)." << std::endl;
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
    double jointDeviation = CalculateJointDeviation(inputenergies, ratioHistograms_3x3, inverseCovarianceMatrix_3x3, meanD2_3x3);
    return GetCDFValue(jointDeviation, hD2_3x3);
  }
  else if (NMATRIXDIM == NMATRIXDIM5)
  {
    double jointDeviation = CalculateJointDeviation(inputenergies, ratioHistograms_5x5, inverseCovarianceMatrix_5x5, meanD2_5x5);
    return GetCDFValue(jointDeviation, hD2_5x5);
  }
  else if (NMATRIXDIM == NMATRIXDIM7)
  {
    double jointDeviation = CalculateJointDeviation(inputenergies, ratioHistograms_7x7, inverseCovarianceMatrix_7x7, meanD2_7x7);
    return GetCDFValue(jointDeviation, hD2_7x7);
  }
  else
  {
    std::cerr << "ClusterCDFCalculator::GetCDF error: Invalid input size. Supported sizes are 3 (for 3x3), 5 (for 5x5), and 7 (for 7x7)." << std::endl;
    return -1;
  }
}

std::unordered_set<int> ClusterCDFCalculator::getSubsetIndices(int _gridsize, int _subsetsize)
{
  std::unordered_set<int> indices;
  if (_subsetsize > _gridsize || _subsetsize % 2 == 0)
  {
    std::cerr << "Invalid subset size!" << std::endl;
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
