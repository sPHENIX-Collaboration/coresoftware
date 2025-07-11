#include "myUtils.h"
#include "geometry_constants.h"

// root includes --
#include <TF1.h>
#include <TFitResult.h>

// c++ includes --
#include <memory>

TFitResultPtr myUtils::doGausFit(TH1 *hist, Double_t start, Double_t end, const std::string &name)
{
  // fit calib hist
  std::unique_ptr<TF1> fitFunc = std::make_unique<TF1>(name.c_str(), "gaus", start, end);

  Double_t initialAmplitude = hist->GetMaximum();
  Double_t initialMean = hist->GetMean();
  Double_t initialSigma = hist->GetRMS();

  fitFunc->SetParameter(0, initialAmplitude);
  fitFunc->SetParameter(1, initialMean);
  fitFunc->SetParameter(2, initialSigma);

  // You can also set parameter names for better readability in the stats box
  fitFunc->SetParName(0, "Amplitude");
  fitFunc->SetParName(1, "Mean");
  fitFunc->SetParName(2, "Sigma");

  // Set some visual properties for the fit line
  fitFunc->SetLineColor(kRed);
  fitFunc->SetLineWidth(2);
  fitFunc->SetLineStyle(kDashed);  // Optional: make it dashed

  TFitResultPtr fitResult = hist->Fit(fitFunc.get(), "RS");  // Fit within range, store result, quiet

  if (fitResult.Get())
  {  // Check if TFitResultPtr is valid
    std::cout << "\n----------------------------------------------------" << std::endl;
    std::cout << "Fit Results for function: " << fitFunc->GetName() << std::endl;
    std::cout << "----------------------------------------------------" << std::endl;
    std::cout << "Fit Status: " << fitResult->Status() << " (0 means successful)" << std::endl;
    if (fitResult->IsValid())
    {  // Check if the fit is valid (e.g., covariance matrix is good)
      std::cout << "Fit is Valid." << std::endl;
      for (UInt_t i = 0; i < static_cast<UInt_t>(fitFunc->GetNpar()); ++i)
      {
        std::cout << "Parameter " << fitFunc->GetParName(static_cast<Int_t>(i)) << " (" << i << "): "
                  << fitResult->Parameter(i) << " +/- " << fitResult->ParError(i) << std::endl;
      }
      std::cout << "Resolution: Sigma/Mean = " << fitResult->Parameter(2) / fitResult->Parameter(1) << std::endl;
      std::cout << "Chi^2 / NDF: " << fitResult->Chi2() << " / " << fitResult->Ndf()
                << " = " << (fitResult->Ndf() > 0 ? fitResult->Chi2() / fitResult->Ndf() : 0) << std::endl;
      std::cout << "Probability: " << TMath::Prob(fitResult->Chi2(), static_cast<Int_t>(fitResult->Ndf())) << std::endl;
    }
    else
    {
      std::cout << "Fit is NOT Valid." << std::endl;
    }
    std::cout << "----------------------------------------------------" << std::endl;
  }
  else
  {
    std::cout << "Fit did not return a valid TFitResultPtr." << std::endl;
  }

  return fitResult;
}

std::vector<std::string> myUtils::split(const std::string &s, const char delimiter)
{
  std::vector<std::string> result;

  std::stringstream ss(s);
  std::string temp;

  while (getline(ss, temp, delimiter))
  {
    if (!temp.empty())
    {
      result.push_back(temp);
    }
  }

  return result;
}

void myUtils::setEMCalDim(TH1 *hist)
{
  hist->GetXaxis()->SetLimits(0, 256);
  hist->GetXaxis()->SetNdivisions(32, false);
  hist->GetXaxis()->SetLabelSize(0.04F);
  hist->GetXaxis()->SetTickSize(0.01F);
  hist->GetYaxis()->SetTickSize(0.01F);
  hist->GetYaxis()->SetLabelSize(0.04F);
  hist->GetYaxis()->SetLimits(0, 96);
  hist->GetYaxis()->SetNdivisions(12, false);
  hist->GetYaxis()->SetTitleOffset(0.5);
  hist->GetXaxis()->SetTitleOffset(1);
}

std::pair<int, int> myUtils::getSectorIB(int iphi, int ieta)
{
  int k = iphi / CaloGeometry::CEMC_NTOW_IB_SIDE;

  int sector = (ieta < 48) ? k + 32 : k;
  int ib = (ieta < 48) ? CaloGeometry::CEMC_NIB_PER_SECTOR - (ieta / CaloGeometry::CEMC_NTOW_IB_SIDE) - 1 : (ieta - 48) / CaloGeometry::CEMC_NTOW_IB_SIDE;

  return std::make_pair(sector, ib);
}

std::pair<int, int> myUtils::getSectorIB(int towerIndex)
{
  int k = towerIndex / CaloGeometry::CEMC_NCHANNEL_PER_SECTOR;

  int sector = (k % 2) ? (k - 1) / 2 : (k / 2) + 32;
  int ib = (k % CaloGeometry::CEMC_NCHANNEL_PER_SECTOR) / CaloGeometry::CEMC_NCHANNEL_PER_IB;

  return std::make_pair(sector, ib);
}
