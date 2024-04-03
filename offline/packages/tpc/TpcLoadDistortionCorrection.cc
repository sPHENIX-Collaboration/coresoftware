/*!
 * \file TpcLoadDistortionCorrection.cc
 * \brief loads distortion correction histogram from file to DistortionCorrectionObject and stores on node tree
 * \author Hugo Pereira Da Costa <hugo.pereira-da-costa@cea.fr>
 */

#include "TpcLoadDistortionCorrection.h"
#include "TpcDistortionCorrectionContainer.h"

#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHDataNode.h>
#include <phool/getClass.h>

#include <TFile.h>
#include <TH1.h>

namespace
{

  // print histogram
  void print_histogram(TH1* h)
  {
    std::cout << "TpcLoadDistortionCorrection::InitRun - name: " << h->GetName() << std::endl;
    for (const auto& axis : {h->GetXaxis(), h->GetYaxis(), h->GetZaxis()})
    {
      if (axis)
      {
        std::cout
            << "  " << axis->GetName()
            << " bins: " << axis->GetNbins()
            << " min: " << axis->GetXmin()
            << " max: " << axis->GetXmax()
            << std::endl;
      }
    }
    std::cout << std::endl;
  }

}  // namespace

//_____________________________________________________________________
TpcLoadDistortionCorrection::TpcLoadDistortionCorrection(const std::string& name)
  : SubsysReco(name)
{
}

//_____________________________________________________________________
int TpcLoadDistortionCorrection::InitRun(PHCompositeNode* topNode)
{
  // look for distortion calibration object
  PHNodeIterator iter(topNode);

  std::cout << "TpcLoadDistortionCorrection::InitRun - m_phi_hist_in_radians: " << m_phi_hist_in_radians << std::endl;

  /// Get the RUN node and check
  auto runNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "RUN"));
  if (!runNode)
  {
    std::cout << "TpcLoadDistortionCorrection::InitRun - RUN Node missing, quitting" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  // create and populate the nodes for each distortion, if present:
  for (int i = 0; i < 3; i++)
  {
    if (!m_correction_in_use[i])
    {
      continue;
    }

    // get distortion correction object and create if not found
    auto distortion_correction_object = findNode::getClass<TpcDistortionCorrectionContainer>(topNode, m_node_name[i]);
    if (!distortion_correction_object)
    {
      std::cout << "TpcLoadDistortionCorrection::InitRun - creating TpcDistortionCorrectionContainer in node " << m_node_name[i] << std::endl;
      distortion_correction_object = new TpcDistortionCorrectionContainer;
      auto node = new PHDataNode<TpcDistortionCorrectionContainer>(distortion_correction_object, m_node_name[i]);
      runNode->addNode(node);
    }

    std::cout << "TpcLoadDistortionCorrection::InitRun - reading corrections from " << m_correction_filename[i] << std::endl;
    auto distortion_tfile = TFile::Open(m_correction_filename[i].c_str());
    if (!distortion_tfile && m_correction_in_use[i])
    {
      std::cout << "TpcLoadDistortionCorrection::InitRun - cannot open " << m_correction_filename[i] << std::endl;
      exit(1);
    }

    const std::array<const std::string, 2> extension = {{"_negz", "_posz"}};
    for (int j = 0; j < 2; ++j)
    {
      distortion_correction_object->m_hDPint[j] = dynamic_cast<TH1*>(distortion_tfile->Get((std::string("hIntDistortionP")+extension[j]).c_str()));
      assert(distortion_correction_object->m_hDPint[j]);
      distortion_correction_object->m_hDRint[j] = dynamic_cast<TH1*>(distortion_tfile->Get((std::string("hIntDistortionR")+extension[j]).c_str()));
      assert(distortion_correction_object->m_hDRint[j]);
      distortion_correction_object->m_hDZint[j] = dynamic_cast<TH1*>(distortion_tfile->Get((std::string("hIntDistortionZ")+extension[j]).c_str()));
      assert(distortion_correction_object->m_hDZint[j]);
    }

    // assign correction object dimension from histograms dimention, assuming all histograms have the same
    distortion_correction_object->m_dimensions = distortion_correction_object->m_hDPint[0]->GetDimension();

    // only dimensions 2 or 3 are supported
    assert(distortion_correction_object->m_dimensions == 2 || distortion_correction_object->m_dimensions == 3);

    // assign whether phi corrections (DP) should be read as radians or mm
    distortion_correction_object->m_phi_hist_in_radians = m_phi_hist_in_radians;

    if (Verbosity())
    {
      for (const auto& h : {
               distortion_correction_object->m_hDPint[0], distortion_correction_object->m_hDPint[1],
               distortion_correction_object->m_hDRint[0], distortion_correction_object->m_hDRint[1],
               distortion_correction_object->m_hDZint[0], distortion_correction_object->m_hDZint[1]})
      {
        print_histogram(h);
      }
    }
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

//_____________________________________________________________________
int TpcLoadDistortionCorrection::process_event(PHCompositeNode* /*unused*/)
{
  return Fun4AllReturnCodes::EVENT_OK;
}
