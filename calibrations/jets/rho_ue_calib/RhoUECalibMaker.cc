#include "RhoUECalibMaker.h"

#include <jetbackground/TowerRho.h>

#include <calobase/TowerInfo.h>
#include <calobase/TowerInfoContainer.h>

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/getClass.h>
#include <phool/phool.h>

#include <TH2D.h>
#include <TH3D.h>

#include <cmath>
#include <iostream>
#include <vector>

const std::array<std::string, RhoUECalibMaker::NSUM> &RhoUECalibMaker::sum_names()
{
  static const std::array<std::string, NSUM> names = {"sw", "sx", "sy", "sxx", "sxy", "syy", "n"};
  return names;
}

void RhoUECalibMaker::book()
{
  std::vector<double> eta_edges(NETA + 1);
  for (int i = 0; i <= NETA; ++i)
  {
    eta_edges[i] = i;
  }
  for (int layer = 0; layer < NLAYER; ++layer)
  {
    for (int k = 0; k < NSUM; ++k)
    {
      const std::string name = sum_names().at(k) + "_" + layer_names().at(layer);
      m_sums.at(layer).at(k) = new TH3D(name.c_str(), (name + ";i#eta;z_{vtx} [cm];Q'").c_str(), NETA, eta_edges.data(),
                                        n_zvtx(), zvtx_edges().data(), n_mbdq(), mbdq_edges().data());
    }
  }
  m_zsum = new TH2D("zsum", "sum of z_{vtx} per class;z_{vtx} [cm];Q'", n_zvtx(), zvtx_edges().data(), n_mbdq(), mbdq_edges().data());
  m_qsum = new TH2D("qsum", "sum of Q' per class;z_{vtx} [cm];Q'", n_zvtx(), zvtx_edges().data(), n_mbdq(), mbdq_edges().data());
  m_nsum = new TH2D("nsum", "events per class;z_{vtx} [cm];Q'", n_zvtx(), zvtx_edges().data(), n_mbdq(), mbdq_edges().data());
}

int RhoUECalibMaker::fill(PHCompositeNode *topNode)
{
  // Q' clamped into the class range, so the centres stay inside their bins
  const double qclamped = std::min(std::max(qprime(), mbdq_edges().front()), mbdq_edges().back() - 1e-9);
  m_zsum->Fill(zvtx(), qclamped, zvtx());
  m_qsum->Fill(zvtx(), qclamped, qclamped);
  m_nsum->Fill(zvtx(), qclamped);

  for (int layer = 0; layer < NLAYER; ++layer)
  {
    auto *towerrho = findNode::getClass<TowerRho>(topNode, m_rho_nodes.at(layer));
    auto *towers = findNode::getClass<TowerInfoContainer>(topNode, tower_node(layer));
    if (!towerrho || !towers)
    {
      std::cout << PHWHERE << " " << Name() << ": " << m_rho_nodes.at(layer) << " or " << tower_node(layer) << " not found" << std::endl;
      return Fun4AllReturnCodes::ABORTRUN;
    }
    if (towerrho->get_method() != TowerRho::Method::MULT)
    {
      std::cout << PHWHERE << " " << Name() << ": " << m_rho_nodes.at(layer) << " is not MULT rho" << std::endl;
      return Fun4AllReturnCodes::ABORTRUN;
    }

    // strip energy and good towers
    std::array<double, NETA> esum{};
    std::array<int, NETA> ngood{};
    for (unsigned int ch = 0; ch < towers->size(); ++ch)
    {
      const unsigned int key = towers->encode_key(ch);
      const int ieta = static_cast<int>(towers->getTowerEtaBin(key));
      if (ieta < 0 || ieta >= NETA)
      {
        continue;
      }
      auto *tower = towers->get_tower_at_channel(ch);
      const float e = tower->get_energy();
      if (tower->get_isGood() && std::isfinite(e))
      {
        esum.at(ieta) += e;
        ++ngood.at(ieta);
      }
    }

    const double rho = towerrho->get_rho();
    auto &sums = m_sums.at(layer);
    for (int ieta = 0; ieta < NETA; ++ieta)
    {
      if (ngood.at(ieta) == 0)
      {
        continue;
      }
      const double n = ngood.at(ieta);
      const double x = rho * std::cosh(eta_corr(layer, ieta));
      const double y = esum.at(ieta) / n;
      const int bin = sums[0]->GetBin(ieta + 1, izvtx() + 1, imbdq() + 1);
      const std::array<double, NSUM> add = {n, n * x, n * y, n * x * x, n * x * y, n * y * y, 1};
      for (int k = 0; k < NSUM; ++k)
      {
        sums.at(k)->AddBinContent(bin, add.at(k));
      }
    }
  }
  return Fun4AllReturnCodes::EVENT_OK;
}
