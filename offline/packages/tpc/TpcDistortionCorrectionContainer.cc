
/*!
 * \file TpcDistortionCorrectionContainer.cc
 * \brief stores distortion correction histograms on the node tree
 * \author Hugo Pereira Da Costa <hugo.pereira-da-costa@cea.fr>
 */

#include "TpcDistortionCorrectionContainer.h"

#include <TFile.h>
#include <TH1.h>
#include <TObject.h>

#include <iostream>
#include <memory>

//_______________________________________________________________
void TpcDistortionCorrectionContainer::load_histograms( const std::string& source )
{
  std::cout << "TpcDistortionCorrectionContainer::load_histograms - reading corrections from " << source << std::endl;
  auto *distortion_tfile = TFile::Open(source.c_str());
  if (!distortion_tfile)
  {
    std::cout << "TpcDistortionCorrectionContainer::load_histograms - cannot open " << source << std::endl;
    exit(1);
  }

  const std::array<const std::string, 2> extension = {{"_negz", "_posz"}};
  for (int j = 0; j < 2; ++j)
  {
    m_hDPint[j] = dynamic_cast<TH1*>(distortion_tfile->Get((std::string("hIntDistortionP")+extension[j]).c_str()));
    assert(m_hDPint[j]);
    m_hDRint[j] = dynamic_cast<TH1*>(distortion_tfile->Get((std::string("hIntDistortionR")+extension[j]).c_str()));
    assert(m_hDRint[j]);
    m_hDZint[j] = dynamic_cast<TH1*>(distortion_tfile->Get((std::string("hIntDistortionZ")+extension[j]).c_str()));
    assert(m_hDZint[j]);
  }
}

//_______________________________________________________________
void TpcDistortionCorrectionContainer::save_histograms( const std::string& destination ) const
{
  // save everything to root file
  std::cout << "TpcDistortionCorrectionContainer::save_histograms - writing histograms to " << destination << std::endl;
  std::unique_ptr<TFile> outputfile(TFile::Open(destination.c_str(), "RECREATE"));
  outputfile->cd();

  for (const auto& h_list : {m_hentries, m_hDRint, m_hDPint, m_hDZint})
  {
    for (const auto& h : h_list)
    {
      if (h)
      {
        h->Write(h->GetName());
      }
    }
  }

  // close TFile
  outputfile->Close();
}
