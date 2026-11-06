#include "MbdTTGainOverride.h"

#include <TF1.h>  // MbdEvent.h holds a std::unique_ptr<TF1> and forward-declares it

#include <mbd/MbdCalib.h>
#include <mbd/MbdDefs.h>
#include <mbd/MbdEvent.h>
#include <mbd/MbdReco.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/Fun4AllServer.h>

#include <phool/phool.h>

#include <iostream>

int MbdTTGainOverride::InitRun(PHCompositeNode * /*topNode*/)
{
  auto *reco = dynamic_cast<MbdReco *>(Fun4AllServer::instance()->getSubsysReco(m_mbdreco));
  MbdCalib *calib = (reco && reco->GetMbdEvent()) ? reco->GetMbdEvent()->GetCalib() : nullptr;
  if (!calib)
  {
    std::cout << PHWHERE << " " << Name() << ": no MbdCalib from " << m_mbdreco << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  const float before = calib->get_ttgain(0);
  for (int ipmt = 0; ipmt < MbdDefs::MBD_N_PMT; ++ipmt)
  {
    calib->set_ttgain(ipmt, m_gain);
  }
  std::cout << Name() << "::InitRun - MBD TTGain " << before << " -> " << m_gain << " for " << MbdDefs::MBD_N_PMT << " PMTs" << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}
