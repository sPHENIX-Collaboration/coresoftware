#include "EventPlaneReco.h"

#include "Eventplaneinfo.h"
#include "EventplaneinfoMap.h"
#include "EventplaneinfoMapv1.h"
#include "Eventplaneinfov1.h"

#include <calobase/TowerInfo.h>
#include <calobase/TowerInfoContainer.h>
#include <calobase/TowerInfoDefs.h>

#include <epd/EpdGeom.h>

#include <mbd/MbdGeom.h>
#include <mbd/MbdOut.h>
#include <mbd/MbdPmtContainer.h>
#include <mbd/MbdPmtHit.h>

#include <globalvertex/GlobalVertex.h>
#include <globalvertex/GlobalVertexMap.h>
#include <globalvertex/MbdVertex.h>
#include <globalvertex/MbdVertexMap.h>

#include <cdbobjects/CDBHistos.h>
#include <cdbobjects/CDBTTree.h>
#include <ffamodules/CDBInterface.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h> // for SubsysReco

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNode.h> // for PHNode
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h> // for PHObject
#include <phool/getClass.h>
#include <phool/phool.h> // for PHWHERE
#include <phool/recoConsts.h>

#include <TProfile2D.h>

#include <array> // for array
#include <cfloat>
#include <cmath>
#include <cstdlib> // for exit
#include <format>
#include <iostream>
#include <set>     // for _Rb_tree_const_iterator
#include <utility> // for pair
#include <vector>  // for vector

EventPlaneReco::EventPlaneReco(const std::string &name) : SubsysReco(name) {

  south_q.resize(m_MaxOrder);
  north_q.resize(m_MaxOrder);
  northsouth_q.resize(m_MaxOrder);

  south_q_subtract.resize(m_MaxOrder);
  north_q_subtract.resize(m_MaxOrder);
  northsouth_q_subtract.resize(m_MaxOrder);

  shift_north.resize(m_MaxOrder);
  shift_south.resize(m_MaxOrder);
  shift_northsouth.resize(m_MaxOrder);
  tmp_south_psi.resize(m_MaxOrder);
  tmp_north_psi.resize(m_MaxOrder);
  tmp_northsouth_psi.resize(m_MaxOrder);

  for (auto &vec : south_q) {
    vec.resize(2);
  }

  for (auto &vec : north_q) {
    vec.resize(2);
  }

  for (auto &vec : northsouth_q) {
    vec.resize(2);
  }

  for (auto &vec : south_q_subtract) {
    vec.resize(2);
  }

  for (auto &vec : north_q_subtract) {
    vec.resize(2);
  }

  for (auto &vec : northsouth_q_subtract) {
    vec.resize(2);
  }

  ring_q_north.resize(nRings);
  ring_q_south.resize(nRings);

  for (auto &rq : ring_q_north) {
    rq.resize(m_MaxOrder, std::vector<double>(2, 0.0));
  }
  for (auto &rq : ring_q_south) {
    rq.resize(m_MaxOrder, std::vector<double>(2, 0.0));
  }

  all_ring_Qvecs_north.assign(
      nRings, std::vector<std::pair<double, double>>(m_MaxOrder, {0.0, 0.0}));

  all_ring_Qvecs_south.assign(
      nRings, std::vector<std::pair<double, double>>(m_MaxOrder, {0.0, 0.0}));
}

int EventPlaneReco::InitRun(PHCompositeNode *topNode) {

  FileName = "EVENTPLANE_CORRECTION";
  if (_isSim) {
    FileName = "EVENTPLANE_CORRECTION_SIM";
  }

  std::string calibdir = CDBInterface::instance()->getUrl(FileName);

  if (calibdir.empty()) {
    std::cout << PHWHERE << "No Eventplane calibration file for domain "
              << FileName << " found" << std::endl;
    std::cout << PHWHERE
              << "Will only produce raw Q vectors and event plane angles "
              << std::endl;
  }

  CDBHistos *cdbhistosIn = new CDBHistos(calibdir);
  cdbhistosIn->LoadCalibrations();

  // Get phiweights
  h_phi_weight_south_input =
      dynamic_cast<TH1F *>(cdbhistosIn->getHisto("h_phi_weight_south", false));
  h_phi_weight_north_input =
      dynamic_cast<TH1F *>(cdbhistosIn->getHisto("h_phi_weight_north", false));

  // Get recentering histograms
  for (unsigned int order = 0; order < m_MaxOrder; order++) {
    tprof_mean_cos_south_epd_input[order] =
        dynamic_cast<TProfile2D *>(cdbhistosIn->getHisto(
            std::format("tprof_mean_cos_south_epd_order_{}", order), false));
    tprof_mean_sin_south_epd_input[order] =
        dynamic_cast<TProfile2D *>(cdbhistosIn->getHisto(
            std::format("tprof_mean_sin_south_epd_order_{}", order), false));
    tprof_mean_cos_north_epd_input[order] =
        dynamic_cast<TProfile2D *>(cdbhistosIn->getHisto(
            std::format("tprof_mean_cos_north_epd_order_{}", order), false));
    tprof_mean_sin_north_epd_input[order] =
        dynamic_cast<TProfile2D *>(cdbhistosIn->getHisto(
            std::format("tprof_mean_sin_north_epd_order_{}", order), false));
    tprof_mean_cos_northsouth_epd_input[order] =
        dynamic_cast<TProfile2D *>(cdbhistosIn->getHisto(
            std::format("tprof_mean_cos_northsouth_epd_order_{}", order),
            false));
    tprof_mean_sin_northsouth_epd_input[order] =
        dynamic_cast<TProfile2D *>(cdbhistosIn->getHisto(
            std::format("tprof_mean_sin_northsouth_epd_order_{}", order),
            false));
  }

  // Get shifting histograms
  for (unsigned int order = 0; order < m_MaxOrder; order++) {
    for (int p = 0; p < _imax; p++) {
      tprof_cos_north_epd_shift_input[order][p] =
          dynamic_cast<TProfile2D *>(cdbhistosIn->getHisto(
              std::format("tprof_cos_north_epd_shift_order_{}_{}", order, p),
              false));
      tprof_sin_north_epd_shift_input[order][p] =
          dynamic_cast<TProfile2D *>(cdbhistosIn->getHisto(
              std::format("tprof_sin_north_epd_shift_order_{}_{}", order, p),
              false));
      tprof_cos_south_epd_shift_input[order][p] =
          dynamic_cast<TProfile2D *>(cdbhistosIn->getHisto(
              std::format("tprof_cos_south_epd_shift_order_{}_{}", order, p),
              false));
      tprof_sin_south_epd_shift_input[order][p] =
          dynamic_cast<TProfile2D *>(cdbhistosIn->getHisto(
              std::format("tprof_sin_south_epd_shift_order_{}_{}", order, p),
              false));
      tprof_cos_northsouth_epd_shift_input[order][p] =
          dynamic_cast<TProfile2D *>(cdbhistosIn->getHisto(
              std::format("tprof_cos_northsouth_epd_shift_order_{}_{}", order,
                          p),
              false));
      tprof_sin_northsouth_epd_shift_input[order][p] =
          dynamic_cast<TProfile2D *>(cdbhistosIn->getHisto(
              std::format("tprof_sin_northsouth_epd_shift_order_{}_{}", order,
                          p),
              false));
    }
  }

  if (Verbosity() > 1) {
    cdbhistosIn->Print();
  }

  return CreateNodes(topNode);
}

int EventPlaneReco::process_event(PHCompositeNode *topNode) {
  if (Verbosity() > 1) {
    std::cout << "EventPlaneReco::process_event -- entered" << std::endl;
  }

  //---------------------------------
  // Get Objects off of the Node Tree
  //---------------------------------

  if (_isSim) {
    // Use GlobalVertexMap for simulation
    GlobalVertexMap *vertexmap =
        findNode::getClass<GlobalVertexMap>(topNode, "GlobalVertexMap");
    if (!vertexmap) {
      std::cout << PHWHERE << "::ERROR - cannot find GlobalVertexMap"
                << std::endl;
      exit(-1);
    }

    if (!vertexmap->empty()) {
      GlobalVertex *vtx = vertexmap->begin()->second;
      if (vtx) {
        _mbdvtx = vtx->get_z();
      }
    }
  } else {
    // Use MbdVertexMap for data
    MbdVertexMap *mbdvtxmap =
        findNode::getClass<MbdVertexMap>(topNode, "MbdVertexMap");
    if (!mbdvtxmap) {
      std::cout << PHWHERE << "::ERROR - cannot find MbdVertexMap" << std::endl;
      exit(-1);
    }

    MbdVertex *mvertex = nullptr;
    if (mbdvtxmap) {
      for (MbdVertexMap::ConstIter mbditer = mbdvtxmap->begin();
           mbditer != mbdvtxmap->end(); ++mbditer) {
        mvertex = mbditer->second;
      }
      if (mvertex) {
        _mbdvtx = mvertex->get_z();
      }
    }
  }

  EventplaneinfoMap *epmap =
      findNode::getClass<EventplaneinfoMap>(topNode, "EventplaneinfoMap");
  if (!epmap) {
    std::cout << PHWHERE << "::ERROR - cannot find EventplaneinfoMap"
              << std::endl;
    exit(-1);
  }

  if (_sepdEpReco) {

    TowerInfoContainer *epd_towerinfo =
        findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_SEPD");
    if (!epd_towerinfo) {
      epd_towerinfo = findNode::getClass<TowerInfoContainer>(
          topNode, "TOWERINFO_CALIB_EPD");
      if (!epd_towerinfo) {
        std::cout << PHWHERE
                  << "::ERROR - cannot find sEPD Calibrated TowerInfoContainer"
                  << std::endl;
        exit(-1);
      }
    }

    EpdGeom *_epdgeom = findNode::getClass<EpdGeom>(topNode, "TOWERGEOM_EPD");
    if (!_epdgeom) {
      std::cout << PHWHERE << "::ERROR - cannot find TOWERGEOM_EPD"
                << std::endl;
      exit(-1);
    }

    ResetMe();

    if ((std::fabs(_mbdvtx) < _mbd_vertex_cut)) {

      unsigned int ntowers = epd_towerinfo->size();
      for (unsigned int ch = 0; ch < ntowers; ch++) {
        TowerInfo *_tower = epd_towerinfo->get_tower_at_channel(ch);
        float epd_e = _tower->get_energy();
        bool isZS = _tower->get_isZS();
        if (!isZS) // exclude ZS
        {
          unsigned int key = TowerInfoDefs::encode_epd(ch);
          int arm = TowerInfoDefs::get_epd_arm(key);
          if (arm == 0) {
            _ssum += epd_e;
          } else if (arm == 1) {
            _nsum += epd_e;
          }
        }
      }

      if (_ssum > _epd_charge_min && _nsum > _epd_charge_min &&
          _ssum < _epd_charge_max && _nsum < _epd_charge_max) {
        _do_ep = true;
      }

      if (_do_ep) {

        // Apply phi weights in builiding ring Q-vectors
        for (unsigned int ch = 0; ch < ntowers; ch++) {
          TowerInfo *_tower = epd_towerinfo->get_tower_at_channel(ch);
          float epd_e = _tower->get_energy();
          bool isZS = _tower->get_isZS();
          if (!isZS) // exclude ZS
          {
            if (epd_e < 0.2) // expecting Nmips
            {
              continue;
            }
            unsigned int key = TowerInfoDefs::encode_epd(ch);
            float tile_phi = _epdgeom->get_phi(key);
            int arm = TowerInfoDefs::get_epd_arm(key);
            int rbin = TowerInfoDefs::get_epd_rbin(key);
            int phibin = TowerInfoDefs::get_epd_phibin(key);

            float truncated_e =
                (epd_e < _epd_e) ? epd_e : _epd_e; // set cutoff at _epd_e

            float TileWeight = truncated_e; // default

            if (h_phi_weight_south_input && h_phi_weight_north_input) {
              if (arm == 0) {
                TileWeight =
                    truncated_e * h_phi_weight_south_input->GetBinContent(
                                      phibin + 1); // scale by 1/<N(φ_{i})>
              } else if (arm == 1) {
                TileWeight =
                    truncated_e * h_phi_weight_north_input->GetBinContent(
                                      phibin + 1); // scale by 1/<N(φ_{i})>
              }
            }

            for (unsigned int order = 0; order < m_MaxOrder; ++order) {
              double Cosine = cos(tile_phi * (double)(order + 1));
              double Sine = sin(tile_phi * (double)(order + 1));

              // Arm-specific Q-vectors
              if (arm == 0) {
                south_q[order][0] += truncated_e * Cosine;
                south_q[order][1] += truncated_e * Sine;
                ring_q_south[rbin][order][0] += TileWeight * Cosine;
                ring_q_south[rbin][order][1] += TileWeight * Sine;

              } else if (arm == 1) {
                north_q[order][0] += truncated_e * Cosine;
                north_q[order][1] += truncated_e * Sine;
                ring_q_north[rbin][order][0] += TileWeight * Cosine;
                ring_q_north[rbin][order][1] += TileWeight * Sine;
              }

              // Combined Q-vectors
              northsouth_q[order][0] += truncated_e * Cosine;
              northsouth_q[order][1] += truncated_e * Sine;
            }
          }
        }

        _totalcharge = _nsum + _ssum;

        // Get recentering histograms and do recentering
        // Recentering: subtract Qn,x and Qn,y values averaged over all events
        for (unsigned int order = 0; order < m_MaxOrder; order++) {
          if (tprof_mean_cos_south_epd_input[order]) // check if recentering
                                                     // histograms exist
          {
            // south
            TAxis *south_xaxis =
                tprof_mean_cos_south_epd_input[order]->GetXaxis();
            TAxis *south_yaxis =
                tprof_mean_cos_south_epd_input[order]->GetYaxis();
            int xbin_south = south_xaxis->FindBin(_ssum);
            int ybin_south = south_yaxis->FindBin(_mbdvtx);

            double event_ave_cos_south =
                tprof_mean_cos_south_epd_input[order]->GetBinContent(
                    xbin_south, ybin_south);
            double event_ave_sin_south =
                tprof_mean_sin_south_epd_input[order]->GetBinContent(
                    xbin_south, ybin_south);
            south_q_subtract[order][0] = _ssum * event_ave_cos_south;
            south_q_subtract[order][1] = _ssum * event_ave_sin_south;
            south_q[order][0] -= south_q_subtract[order][0];
            south_q[order][1] -= south_q_subtract[order][1];

            // north
            TAxis *north_xaxis =
                tprof_mean_cos_north_epd_input[order]->GetXaxis();
            TAxis *north_yaxis =
                tprof_mean_cos_north_epd_input[order]->GetYaxis();
            int xbin_north = north_xaxis->FindBin(_nsum);
            int ybin_north = north_yaxis->FindBin(_mbdvtx);

            double event_ave_cos_north =
                tprof_mean_cos_north_epd_input[order]->GetBinContent(
                    xbin_north, ybin_north);
            double event_ave_sin_north =
                tprof_mean_sin_north_epd_input[order]->GetBinContent(
                    xbin_north, ybin_north);
            north_q_subtract[order][0] = _nsum * event_ave_cos_north;
            north_q_subtract[order][1] = _nsum * event_ave_sin_north;
            north_q[order][0] -= north_q_subtract[order][0];
            north_q[order][1] -= north_q_subtract[order][1];

            // northsouth
            TAxis *northsouth_xaxis =
                tprof_mean_cos_northsouth_epd_input[order]->GetXaxis();
            TAxis *northsouth_yaxis =
                tprof_mean_cos_northsouth_epd_input[order]->GetYaxis();
            int xbin_northsouth = northsouth_xaxis->FindBin(_totalcharge);
            int ybin_northsouth = northsouth_yaxis->FindBin(_mbdvtx);

            double event_ave_cos_northsouth =
                tprof_mean_cos_northsouth_epd_input[order]->GetBinContent(
                    xbin_northsouth, ybin_northsouth);
            double event_ave_sin_northsouth =
                tprof_mean_sin_northsouth_epd_input[order]->GetBinContent(
                    xbin_northsouth, ybin_northsouth);
            northsouth_q_subtract[order][0] =
                _totalcharge * event_ave_cos_northsouth;
            northsouth_q_subtract[order][1] =
                _totalcharge * event_ave_sin_northsouth;
            northsouth_q[order][0] -= northsouth_q_subtract[order][0];
            northsouth_q[order][1] -= northsouth_q_subtract[order][1];
          }
        }

        // Get recentered psi_n
        Eventplaneinfo *epinfo = new Eventplaneinfov1();
        for (unsigned int order = 0; order < m_MaxOrder; order++) {
          double n = order + 1.0;
          if (tprof_mean_cos_south_epd_input[order]) // if present, Qs are
                                                     // recentered
          {
            tmp_south_psi[order] =
                epinfo->GetPsi(south_q[order][0], south_q[order][1], n);
            tmp_north_psi[order] =
                epinfo->GetPsi(north_q[order][0], north_q[order][1], n);
            tmp_northsouth_psi[order] = epinfo->GetPsi(
                northsouth_q[order][0], northsouth_q[order][1], n);
          } else {
            tmp_south_psi[order] = NAN;
            tmp_north_psi[order] = NAN;
            tmp_northsouth_psi[order] = NAN;
          }
        }

        // Get shifting histograms and calculate shift
        for (unsigned int order = 0; order < m_MaxOrder; order++) {
          for (int p = 0; p < _imax; p++) {
            if (tprof_cos_south_epd_shift_input[order][p]) // check if shifting
                                                           // histograms exist
            {
              double terms = p + 1.0;
              double n = order + 1.0;
              double tmp = n * terms;
              double prefactor = 2.0 / terms;

              // south
              TAxis *south_xaxis =
                  tprof_cos_south_epd_shift_input[order][p]->GetXaxis();
              TAxis *south_yaxis =
                  tprof_cos_south_epd_shift_input[order][p]->GetYaxis();
              int xbin_south = south_xaxis->FindBin(_ssum);
              int ybin_south = south_yaxis->FindBin(_mbdvtx);

              // north
              TAxis *north_xaxis =
                  tprof_cos_north_epd_shift_input[order][p]->GetXaxis();
              TAxis *north_yaxis =
                  tprof_cos_north_epd_shift_input[order][p]->GetYaxis();
              int xbin_north = north_xaxis->FindBin(_nsum);
              int ybin_north = north_yaxis->FindBin(_mbdvtx);

              //                   // northsouth
              TAxis *northsouth_xaxis =
                  tprof_cos_northsouth_epd_shift_input[order][p]->GetXaxis();
              TAxis *northsouth_yaxis =
                  tprof_cos_northsouth_epd_shift_input[order][p]->GetYaxis();
              int xbin_northsouth = northsouth_xaxis->FindBin(_totalcharge);
              int ybin_northsouth = northsouth_yaxis->FindBin(_mbdvtx);

              //                    Equation (6) of arxiv:nucl-ex/9805001
              //                    i = terms; n = order; i*n = tmp
              //                    (2 / i ) * <cos(i*n*psi_n)> * sin(i*n*psi_n)
              //                    - <sin(i*n*psi_n)> * cos(i*n*psi_n)

              // north
              shift_north[order] +=
                  prefactor *
                  (tprof_cos_north_epd_shift_input[order][p]->GetBinContent(
                       xbin_north, ybin_north) *
                       sin(tmp * tmp_north_psi[order]) -
                   tprof_sin_north_epd_shift_input[order][p]->GetBinContent(
                       xbin_north, ybin_north) *
                       cos(tmp * tmp_north_psi[order]));

              // south
              shift_south[order] +=
                  prefactor *
                  (tprof_cos_south_epd_shift_input[order][p]->GetBinContent(
                       xbin_south, ybin_south) *
                       sin(tmp * tmp_south_psi[order]) -
                   tprof_sin_south_epd_shift_input[order][p]->GetBinContent(
                       xbin_south, ybin_south) *
                       cos(tmp * tmp_south_psi[order]));

              //                     // northsouth
              shift_northsouth[order] +=
                  prefactor *
                  (tprof_cos_northsouth_epd_shift_input[order][p]
                           ->GetBinContent(xbin_northsouth, ybin_northsouth) *
                       sin(tmp * tmp_northsouth_psi[order]) -
                   tprof_sin_northsouth_epd_shift_input[order][p]
                           ->GetBinContent(xbin_northsouth, ybin_northsouth) *
                       cos(tmp * tmp_northsouth_psi[order]));
            }
          }
        }

        // n * deltapsi_n = (2 / i ) * <cos(i*n*psi_n)> * sin(i*n*psi_n) -
        // <sin(i*n*psi_n)> * cos(i*n*psi_n) Divide out n
        for (unsigned int order = 0; order < m_MaxOrder; order++) {
          double n = order + 1.0;
          shift_north[order] /= n;
          shift_south[order] /= n;
          shift_northsouth[order] /= n;
        }

        // Now add shift to psi_n to flatten it
        for (unsigned int order = 0; order < m_MaxOrder; order++) {
          if (tprof_cos_north_epd_shift_input[0][0]) {
            tmp_south_psi[order] += shift_south[order];
            tmp_north_psi[order] += shift_north[order];
            tmp_northsouth_psi[order] += shift_northsouth[order];
          }
        }

        // Now enforce the range
        for (unsigned int order = 0; order < m_MaxOrder; order++) {
          if (tprof_cos_north_epd_shift_input[0][0]) {
            double range = M_PI / (double)(order + 1);
            if (tmp_south_psi[order] < -1.0 * range) {
              tmp_south_psi[order] += 2.0 * range;
            }
            if (tmp_south_psi[order] > range) {
              tmp_south_psi[order] -= 2.0 * range;
            }
            if (tmp_north_psi[order] < -1.0 * range) {
              tmp_north_psi[order] += 2.0 * range;
            }
            if (tmp_north_psi[order] > range) {
              tmp_north_psi[order] -= 2.0 * range;
            }
            if (tmp_northsouth_psi[order] < -1.0 * range) {
              tmp_northsouth_psi[order] += 2.0 * range;
            }
            if (tmp_northsouth_psi[order] > range) {
              tmp_northsouth_psi[order] -= 2.0 * range;
            }
          }
        }

        for (unsigned int order = 0; order < m_MaxOrder; order++) {
          south_Qvec.emplace_back(south_q[order][0], south_q[order][1]);
          north_Qvec.emplace_back(north_q[order][0], north_q[order][1]);
          northsouth_Qvec.emplace_back(northsouth_q[order][0],
                                       northsouth_q[order][1]);
        }

        for (int rbin = 0; rbin < nRings; ++rbin) {
          for (unsigned int order = 0; order < m_MaxOrder; ++order) {
            all_ring_Qvecs_north[rbin][order] = std::make_pair(
                ring_q_north[rbin][order][0], ring_q_north[rbin][order][1]);
            all_ring_Qvecs_south[rbin][order] = std::make_pair(
                ring_q_south[rbin][order][0], ring_q_south[rbin][order][1]);
          }
        }

        if (epd_towerinfo) {
          Eventplaneinfo *sepds = new Eventplaneinfov1();
          sepds->set_qvector(south_Qvec);
          sepds->set_shifted_psi(tmp_south_psi);
          epmap->insert(sepds, EventplaneinfoMap::sEPDS);

          Eventplaneinfo *sepdn = new Eventplaneinfov1();
          sepdn->set_qvector(north_Qvec);
          sepdn->set_shifted_psi(tmp_north_psi);
          epmap->insert(sepdn, EventplaneinfoMap::sEPDN);

          Eventplaneinfo *sepdns = new Eventplaneinfov1();
          sepdns->set_qvector(northsouth_Qvec);
          sepdns->set_shifted_psi(tmp_northsouth_psi);
          epmap->insert(sepdns, EventplaneinfoMap::sEPDNS);

          Eventplaneinfo *epring_south = new Eventplaneinfov1();
          epring_south->set_ring_qvector(all_ring_Qvecs_south);
          epmap->insert(epring_south, EventplaneinfoMap::sEPDRING_SOUTH);

          Eventplaneinfo *epring_north = new Eventplaneinfov1();
          epring_north->set_ring_qvector(all_ring_Qvecs_north);
          epmap->insert(epring_north, EventplaneinfoMap::sEPDRING_NORTH);

          if (Verbosity() > 1) {
            sepds->identify();
            sepdn->identify();
            sepdns->identify();
            epring_south->identify();
            epring_north->identify();
          }
        }
      }
    }
  }

  if (_mbdEpReco) {
    ResetMe();

    MbdPmtContainer *mbdpmts =
        findNode::getClass<MbdPmtContainer>(topNode, "MbdPmtContainer");
    if (!mbdpmts) {
      std::cout << PHWHERE << "::ERROR - cannot find MbdPmtContainer"
                << std::endl;
      exit(-1);
    }

    MbdGeom *mbdgeom = findNode::getClass<MbdGeom>(topNode, "MbdGeom");
    if (!mbdgeom) {
      std::cout << PHWHERE << "::ERROR - cannot find MbdGeom" << std::endl;
      exit(-1);
    }

    if (mbdpmts) {
      if (Verbosity()) {
        std::cout << "EventPlaneReco::process_event -  mbdpmts" << std::endl;
      }

      for (int ipmt = 0; ipmt < mbdpmts->get_npmt(); ipmt++) {
        float mbd_q = mbdpmts->get_pmt(ipmt)->get_q();
        _mbdQ += mbd_q;
      }

      for (int ipmt = 0; ipmt < mbdpmts->get_npmt(); ipmt++) {
        float mbd_q = mbdpmts->get_pmt(ipmt)->get_q();
        float phi = mbdgeom->get_phi(ipmt);
        int arm = mbdgeom->get_arm(ipmt);

        if (_mbdQ < _mbd_e) {
          continue;
        }

        if (arm == 0) {
          for (unsigned int order = 0; order < m_MaxOrder; order++) {
            double Cosine = cos(phi * (double)(order + 1));
            double Sine = sin(phi * (double)(order + 1));
            south_q[order][0] += mbd_q * Cosine; // south Qn,x
            south_q[order][1] += mbd_q * Sine;   // south Qn,y
          }
        } else if (arm == 1) {
          for (unsigned int order = 0; order < m_MaxOrder; order++) {
            double Cosine = cos(phi * (double)(order + 1));
            double Sine = sin(phi * (double)(order + 1));
            north_q[order][0] += mbd_q * Cosine; // north Qn,x
            north_q[order][1] += mbd_q * Sine;   // north Qn,y
          }
        }
      }
    }

    for (unsigned int order = 0; order < m_MaxOrder; order++) {
      south_Qvec.emplace_back(south_q[order][0], south_q[order][1]);
      north_Qvec.emplace_back(north_q[order][0], north_q[order][1]);
    }

    if (mbdpmts) {
      Eventplaneinfo *mbds = new Eventplaneinfov1();
      mbds->set_qvector(south_Qvec);
      epmap->insert(mbds, EventplaneinfoMap::MBDS);

      Eventplaneinfo *mbdn = new Eventplaneinfov1();
      mbdn->set_qvector(north_Qvec);
      epmap->insert(mbdn, EventplaneinfoMap::MBDN);

      if (Verbosity() > 1) {
        mbds->identify();
        mbdn->identify();
      }
    }

    ResetMe();
  }

  if (Verbosity()) {
    epmap->identify();
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int EventPlaneReco::CreateNodes(PHCompositeNode *topNode) {
  PHNodeIterator iter(topNode);

  PHCompositeNode *dstNode =
      dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode) {
    std::cout << PHWHERE << "DST Node missing, doing nothing." << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  PHCompositeNode *globalNode = dynamic_cast<PHCompositeNode *>(
      iter.findFirst("PHCompositeNode", "GLOBAL"));
  if (!globalNode) {
    globalNode = new PHCompositeNode("GLOBAL");
    dstNode->addNode(globalNode);
  }

  EventplaneinfoMap *eps =
      findNode::getClass<EventplaneinfoMap>(topNode, "EventplaneinfoMap");
  if (!eps) {
    eps = new EventplaneinfoMapv1();
    PHIODataNode<PHObject> *EpMapNode =
        new PHIODataNode<PHObject>(eps, "EventplaneinfoMap", "PHObject");
    globalNode->addNode(EpMapNode);
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

void EventPlaneReco::ResetMe() {
  for (auto &vec : south_q) {
    std::fill(vec.begin(), vec.end(), 0.);
  }

  for (auto &vec : north_q) {
    std::fill(vec.begin(), vec.end(), 0.);
  }

  for (auto &vec : northsouth_q) {
    std::fill(vec.begin(), vec.end(), 0.);
  }

  for (auto &order_vec : ring_q_north) {
    for (auto &xy_vec : order_vec) {
      std::fill(xy_vec.begin(), xy_vec.end(), 0.0);
    }
  }

  for (auto &order_vec : ring_q_south) {
    for (auto &xy_vec : order_vec) {
      std::fill(xy_vec.begin(), xy_vec.end(), 0.0);
    }
  }

  south_Qvec.clear();
  north_Qvec.clear();
  northsouth_Qvec.clear();

  for (auto &ring : all_ring_Qvecs_north) {
    for (auto &q : ring) {
      q = {0.0, 0.0};
    }
  }

  for (auto &ring : all_ring_Qvecs_south) {
    for (auto &q : ring) {
      q = {0.0, 0.0};
    }
  }

  for (auto &vec : south_q_subtract) {
    std::fill(vec.begin(), vec.end(), 0.);
  }

  for (auto &vec : north_q_subtract) {
    std::fill(vec.begin(), vec.end(), 0.);
  }

  for (auto &vec : northsouth_q_subtract) {
    std::fill(vec.begin(), vec.end(), 0.);
  }

  std::fill(shift_north.begin(), shift_north.end(), 0.);
  std::fill(shift_south.begin(), shift_south.end(), 0.);
  std::fill(shift_northsouth.begin(), shift_northsouth.end(), 0.);

  std::fill(tmp_south_psi.begin(), tmp_south_psi.end(), NAN);
  std::fill(tmp_north_psi.begin(), tmp_north_psi.end(), NAN);
  std::fill(tmp_northsouth_psi.begin(), tmp_northsouth_psi.end(), NAN);

  _nsum = 0.;
  _ssum = 0.;
  _do_ep = false;
  _mbdQ = 0.;
  _totalcharge = 0.;
}

int EventPlaneReco::End(PHCompositeNode * /*topNode*/) {

  std::cout << " EventPlaneReco::End() " << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}
