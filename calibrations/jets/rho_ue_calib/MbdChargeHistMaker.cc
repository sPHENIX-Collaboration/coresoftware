#include "MbdChargeHistMaker.h"

#include <fun4all/Fun4AllReturnCodes.h>

#include <TDirectory.h>
#include <TH1D.h>

#include <cmath>
#include <string>

int MbdChargeHistMaker::fill(PHCompositeNode * /*topNode*/)
{
  // one histogram per run, booked at its first event (the output file is the current directory)
  const std::string name = "q_" + std::to_string(run());
  auto *h = dynamic_cast<TH1D *>(gDirectory->Get(name.c_str()));
  if (!h)
  {
    const int nbins = static_cast<int>(std::lround(m_qmax));
    h = new TH1D(name.c_str(), ("MBD charge, run " + std::to_string(run()) + ";Q = q_{S} + q_{N};events").c_str(), nbins, 0, m_qmax);
  }
  h->Fill(mbdq());
  const std::string t0name = "t0_" + std::to_string(run());
  auto *ht0 = dynamic_cast<TH1D *>(gDirectory->Get(t0name.c_str()));
  if (!ht0)
  {
    ht0 = new TH1D(t0name.c_str(), ("MBD t0, run " + std::to_string(run()) + ";t0 [ns];events").c_str(), 1000, -25, 25);
  }
  ht0->Fill(mbdt0());
  return Fun4AllReturnCodes::EVENT_OK;
}
