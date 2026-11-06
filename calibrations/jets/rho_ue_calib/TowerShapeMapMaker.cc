#include "TowerShapeMapMaker.h"

#include <calobase/TowerInfo.h>
#include <calobase/TowerInfoContainer.h>

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/getClass.h>
#include <phool/phool.h>

#include <TH1D.h>
#include <TH3F.h>

#include <cmath>
#include <iostream>

void TowerShapeMapMaker::book()
{
  const int nc = n_class();
  for (int layer = 0; layer < NLAYER; ++layer)
  {
    for (int half = 0; half < 2; ++half)
    {
      const std::string tag = layer_names().at(layer) + "_h" + std::to_string(half);
      m_sum.at(layer).at(half) = new TH3F(("sum_" + tag).c_str(), ("sum of tower E_{T} [GeV], " + tag + ";i#eta;i#phi;class").c_str(),
                                          NETA, 0, NETA, NPHI, 0, NPHI, nc, 0, nc);
      m_cnt.at(layer).at(half) = new TH3F(("cnt_" + tag).c_str(), ("good tower entries, " + tag + ";i#eta;i#phi;class").c_str(),
                                          NETA, 0, NETA, NPHI, 0, NPHI, nc, 0, nc);
    }
  }
  m_nevt = new TH1D("h_nevt", "selected events per class;class = izbin * n_{Q'} + iQ';events", nc, 0, nc);
}

int TowerShapeMapMaker::fill(PHCompositeNode *topNode)
{
  m_nevt->Fill(event_class());
  const int half = evt_index() % 2;
  for (int layer = 0; layer < NLAYER; ++layer)
  {
    auto *towers = findNode::getClass<TowerInfoContainer>(topNode, tower_node(layer));
    if (!towers)
    {
      std::cout << PHWHERE << " " << Name() << ": tower node " << tower_node(layer) << " not found" << std::endl;
      return Fun4AllReturnCodes::ABORTRUN;
    }
    TH3F *hsum = m_sum.at(layer).at(half);
    TH3F *hcnt = m_cnt.at(layer).at(half);
    for (unsigned int ch = 0; ch < towers->size(); ++ch)
    {
      const unsigned int key = towers->encode_key(ch);
      const int ieta = static_cast<int>(towers->getTowerEtaBin(key));
      const int iphi = static_cast<int>(towers->getTowerPhiBin(key));
      if (ieta < 0 || ieta >= NETA || iphi < 0 || iphi >= NPHI)
      {
        continue;
      }
      auto *tower = towers->get_tower_at_channel(ch);
      const float e = tower->get_energy();
      if (!tower->get_isGood() || !std::isfinite(e))
      {
        continue;
      }
      const int bin = hsum->GetBin(ieta + 1, iphi + 1, event_class() + 1);
      hsum->AddBinContent(bin, e / std::cosh(eta_corr(layer, ieta)));
      hcnt->AddBinContent(bin, 1);
    }
  }
  return Fun4AllReturnCodes::EVENT_OK;
}
