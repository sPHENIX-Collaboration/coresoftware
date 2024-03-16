#include "EventPlaneReco.h"

#include "Eventplaneinfo.h"
#include "EventplaneinfoMap.h"
#include "EventplaneinfoMapv1.h"
#include "Eventplaneinfov1.h"

#include <calobase/TowerInfo.h>
#include <calobase/TowerInfoContainer.h>
#include <calobase/TowerInfoDefs.h>

#include <globalvertex/GlobalVertex.h>
#include <globalvertex/GlobalVertexMap.h>

#include <epd/EpdGeom.h>

#include <mbd/MbdGeom.h>
#include <mbd/MbdOut.h>
#include <mbd/MbdPmtContainer.h>
#include <mbd/MbdPmtHit.h>

#include <cdbobjects/CDBHistos.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>  // for SubsysReco

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNode.h>  // for PHNode
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>  // for PHObject
#include <phool/getClass.h>
#include <phool/phool.h>  // for PHWHERE

#include <array>  // for array
#include <cfloat>
#include <cmath>
#include <cstdlib>  // for exit
#include <iostream>
#include <set>      // for _Rb_tree_const_iterator
#include <utility>  // for pair
#include <vector>   // for vector
#include <TProfile.h>
#include <TH1D.h>
#include <boost/format.hpp>

EventPlaneReco::EventPlaneReco(const std::string &name)
  : SubsysReco(name)
  , OutFileName("eventplane_correction_histograms")
{
  mbd_e_south = 0.;
  mbd_e_north = 0.;
  mbdQ = 0.;
  tmp_south_psi.resize(m_MaxOrder);
  tmp_north_psi.resize(m_MaxOrder);
  shift_north.resize(m_MaxOrder);
  shift_south.resize(m_MaxOrder);
  south_q.resize(m_MaxOrder);
  north_q.resize(m_MaxOrder);
  south_q_subtract.resize(m_MaxOrder);
  north_q_subtract.resize(m_MaxOrder);

  for (auto &vec : south_q)
  {
    vec.resize(2);
  }

  for (auto &vec : north_q)
  {
    vec.resize(2);
  }
    
 for (auto &vec : south_q_subtract)
 {
    vec.resize(2);
 }

 for (auto &vec : north_q_subtract)
 {
    vec.resize(2);
 }
    
}

int EventPlaneReco::Init(PHCompositeNode *topNode)
{

    OutFileName =  boost::str(boost::format("eventplane_correction_histograms_run_%d.root") %m_runNo).c_str();

    //-----------------------------------load calibration histograms-----------------------------------------//
    CDBHistos *cdbhistosIn = new CDBHistos(OutFileName);
    cdbhistosIn->LoadCalibrations();
    
    // Get recentering histograms
    for (unsigned int order = 1; order < m_MaxOrder + 1; order++)
    {
        if(!(cdbhistosIn->getHisto(boost::str(boost::format("tprof_mean_cos_south_mbd_order_%d") %order).c_str())))
        {
            tprof_mean_cos_south_mbd_input[order - 1] = NULL;
            tprof_mean_sin_south_mbd_input[order - 1] =  NULL;
            tprof_mean_cos_north_mbd_input[order - 1] =  NULL;
            tprof_mean_sin_north_mbd_input[order - 1] =  NULL;
        }
        else
        {
            tprof_mean_cos_south_mbd_input[order - 1] = dynamic_cast<TProfile *> (cdbhistosIn->getHisto(boost::str(boost::format("tprof_mean_cos_south_mbd_order_%d") %order).c_str()));
            tprof_mean_sin_south_mbd_input[order - 1] = dynamic_cast<TProfile *> (cdbhistosIn->getHisto(boost::str(boost::format("tprof_mean_sin_south_mbd_order_%d") %order).c_str()));
            tprof_mean_cos_north_mbd_input[order - 1] = dynamic_cast<TProfile *> (cdbhistosIn->getHisto(boost::str(boost::format("tprof_mean_cos_north_mbd_order_%d") %order).c_str()));
            tprof_mean_sin_north_mbd_input[order - 1] = dynamic_cast<TProfile *> (cdbhistosIn->getHisto(boost::str(boost::format("tprof_mean_sin_north_mbd_order_%d") %order).c_str()));
        }
    }

    // Get shifting histograms
    for (unsigned int order = 1; order < m_MaxOrder + 1; order++)
    {
        for ( int p = 0; p < _imax; p++ )
        {
            if(!(cdbhistosIn->getHisto(boost::str(boost::format("tprof_cos_north_mbd_shift_order_%d_%d") %order %p).c_str())))
            {
                tprof_cos_north_mbd_shift_input[order - 1][p] = NULL;
                tprof_sin_north_mbd_shift_input[order - 1][p] =  NULL;
                tprof_cos_south_mbd_shift_input[order - 1][p] =  NULL;
                tprof_sin_south_mbd_shift_input[order - 1][p] =  NULL;
            }
            else
            {
                tprof_cos_north_mbd_shift_input[order - 1][p] = dynamic_cast<TProfile *> (cdbhistosIn->getHisto(boost::str(boost::format("tprof_cos_north_mbd_shift_order_%d_%d") %order %p).c_str()));
                tprof_sin_north_mbd_shift_input[order - 1][p] = dynamic_cast<TProfile *> (cdbhistosIn->getHisto(boost::str(boost::format("tprof_sin_north_mbd_shift_order_%d_%d") %order %p).c_str()));
                tprof_cos_south_mbd_shift_input[order - 1][p] = dynamic_cast<TProfile *> (cdbhistosIn->getHisto(boost::str(boost::format("tprof_cos_south_mbd_shift_order_%d_%d") %order %p).c_str()));
                tprof_sin_south_mbd_shift_input[order - 1][p] = dynamic_cast<TProfile *> (cdbhistosIn->getHisto(boost::str(boost::format("tprof_sin_south_mbd_shift_order_%d_%d") %order %p).c_str()));
            }
        }
      }
    
    cdbhistosIn->Print();
    //-----------------------------------create correction histograms-----------------------------------------//

    cdbhistosOut = new CDBHistos(OutFileName);

    // Create and register recentering histograms
    for (unsigned int order = 1; order < m_MaxOrder + 1; order++)
    {
        //total charge in MBD in data max out ~2500, sim ~ factor 100
        tprof_mean_cos_south_mbd[order - 1] =  new TProfile(boost::str(boost::format("tprof_mean_cos_south_mbd_order_%d") %order).c_str(), "", 125*400, 0, 250000, -1e10, 1e10);
        tprof_mean_sin_south_mbd[order - 1] =  new TProfile(boost::str(boost::format("tprof_mean_sin_south_mbd_order_%d") %order).c_str(), "", 125*400, 0, 250000, -1e10, 1e10);
        tprof_mean_cos_north_mbd[order - 1] =  new TProfile(boost::str(boost::format("tprof_mean_cos_north_mbd_order_%d") %order).c_str(), "", 125*400, 0, 250000, -1e10, 1e10);
        tprof_mean_sin_north_mbd[order - 1] =  new TProfile(boost::str(boost::format("tprof_mean_sin_north_mbd_order_%d") %order).c_str(), "", 125*400, 0, 250000, -1e10, 1e10);
       
        cdbhistosOut->registerHisto(tprof_mean_cos_south_mbd[order - 1]);
        cdbhistosOut->registerHisto(tprof_mean_sin_south_mbd[order - 1]);
        cdbhistosOut->registerHisto(tprof_mean_cos_north_mbd[order - 1]);
        cdbhistosOut->registerHisto(tprof_mean_sin_north_mbd[order - 1]);
        
    }

    // Create and register shifting histograms
    for (unsigned int order = 1; order < m_MaxOrder + 1; order++)
    {
        // Request that recentering histograms already exists before creating
        if(cdbhistosIn->getHisto(boost::str(boost::format("tprof_mean_cos_south_mbd_order_%d") %order).c_str()))
        {
            for ( int p = 0; p < _imax; p++ )
            {
                tprof_cos_north_mbd_shift[order - 1][p] =  new TProfile(boost::str(boost::format("tprof_cos_north_mbd_shift_order_%d_%d") %order %p).c_str(), "", 125*400, 0, 250000, -1e10, 1e10);
                tprof_sin_north_mbd_shift[order - 1][p] =  new TProfile(boost::str(boost::format("tprof_sin_north_mbd_shift_order_%d_%d") %order %p).c_str(), "", 125*400, 0, 250000, -1e10, 1e10);
                tprof_cos_south_mbd_shift[order - 1][p] =  new TProfile(boost::str(boost::format("tprof_cos_south_mbd_shift_order_%d_%d") %order %p).c_str(), "", 125*400, 0, 250000, -1e10, 1e10);
                tprof_sin_south_mbd_shift[order - 1][p] =  new TProfile(boost::str(boost::format("tprof_sin_south_mbd_shift_order_%d_%d") %order %p).c_str(), "", 125*400, 0, 250000, -1e10, 1e10);
                
                cdbhistosOut->registerHisto(tprof_cos_north_mbd_shift[order - 1][p]);
                cdbhistosOut->registerHisto(tprof_sin_north_mbd_shift[order - 1][p]);
                cdbhistosOut->registerHisto(tprof_cos_south_mbd_shift[order - 1][p]);
                cdbhistosOut->registerHisto(tprof_sin_south_mbd_shift[order - 1][p]);
            }
        }
        
    }
 
    hvertex = new TH1D("hvertex","z-vertex [cm]",2100,-1050,1050);
    cdbhistosOut->registerHisto(hvertex);
 
    return CreateNodes(topNode);
    
}

int EventPlaneReco::InitRun(PHCompositeNode *topNode)
{
  if (Verbosity() > 0)
  {
    std::cout << "======================= EventPlaneReco::InitRun() =======================" << std::endl;
    std::cout << "===========================================================================" << std::endl;
  }

  return CreateNodes(topNode);
}

int EventPlaneReco::process_event(PHCompositeNode *topNode)
{
 
  if (Verbosity() > 1)
  {
    std::cout << "EventPlaneReco::process_event -- entered" << std::endl;
  }

  //---------------------------------
  // Get Objects off of the Node Tree
  //---------------------------------
    
  GlobalVertexMap *vertexmap = findNode::getClass<GlobalVertexMap>(topNode, "GlobalVertexMap");
  if (!vertexmap)
  {
    std::cout << PHWHERE << " Fatal Error - GlobalVertexMap node is missing" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
    
  if(vertexmap->empty())
  {
    _mbd_vertex = 999.0; //track events with empty vertices
  }
  else
  {
    GlobalVertex *vtx = vertexmap->begin()->second;
    if(vtx)
    {   
      _mbd_vertex = vtx->get_z();
    }
  }

  hvertex->Fill(_mbd_vertex);
  EventplaneinfoMap *epmap = findNode::getClass<EventplaneinfoMap>(topNode, "EventplaneinfoMap");
  if (!epmap)
  {
    std::cout << PHWHERE << "::ERROR - cannot find EventplaneinfoMap" << std::endl;
    exit(-1);
  }

  if (_sepdEpReco)
  {
    ResetMe();
    TowerInfoContainer *epd_towerinfo = findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_EPD");
    if (!epd_towerinfo)
    {
      std::cout << PHWHERE << "::ERROR - cannot find TOWERINFO_CALIB_EPD" << std::endl;
      exit(-1);
    }
    EpdGeom *_epdgeom = findNode::getClass<EpdGeom>(topNode, "TOWERGEOM_EPD");
    if (!_epdgeom)
    {
      std::cout << PHWHERE << "::ERROR - cannot find TOWERGEOM_EPD" << std::endl;
      exit(-1);
    }

    if (epd_towerinfo)
    {
      if (Verbosity())
      {
        std::cout << "EventPlaneReco::process_event -  epd_towerinfo" << std::endl;
      }

      unsigned int ntowers = epd_towerinfo->size();
      for (unsigned int ch = 0; ch < ntowers; ch++)
      {
        TowerInfo *_tower = epd_towerinfo->get_tower_at_channel(ch);
        unsigned int key = TowerInfoDefs::encode_epd(ch);
        float epd_e = _tower->get_energy();
        if (epd_e < 0.5)
        {
          continue;
        }
        float tile_phi = _epdgeom->get_phi(key);
        int arm = TowerInfoDefs::get_epd_arm(key);
        float truncated_e = (epd_e < _epd_e) ? epd_e : _epd_e;
        if (arm == 0)
        {
          for (unsigned int order = 1; order < m_MaxOrder + 1; order++)
          {
            double Cosine = cos(tile_phi * (double) order);
            double Sine = sin(tile_phi * (double) order);
            south_q[order - 1][0] += truncated_e * Cosine;  // south Qn,x
            south_q[order - 1][1] += truncated_e * Sine;    // south Qn,y
          }
        }
        else if (arm == 1)
        {
          for (unsigned int order = 1; order < m_MaxOrder + 1; order++)
          {
            double Cosine = cos(tile_phi * (double) order);
            double Sine = sin(tile_phi * (double) order);
            north_q[order - 1][0] += truncated_e * Cosine;  // north Qn,x
            north_q[order - 1][1] += truncated_e * Sine;    // north Qn,y
          }
        }
      }
    }
    for (unsigned int order = 1; order < m_MaxOrder + 1; order++)
    {
      south_Qvec.emplace_back(south_q[order - 1][0], south_q[order - 1][1]);
      north_Qvec.emplace_back(north_q[order - 1][0], north_q[order - 1][1]);
    }

    if (epd_towerinfo)
    {
      Eventplaneinfo *sepds = new Eventplaneinfov1();
      sepds->set_qvector(south_Qvec);
      epmap->insert(sepds, EventplaneinfoMap::sEPDS);

      Eventplaneinfo *sepdn = new Eventplaneinfov1();
      sepdn->set_qvector(north_Qvec);
      epmap->insert(sepdn, EventplaneinfoMap::sEPDN);

      if (Verbosity() > 1)
      {
        sepds->identify();
        sepdn->identify();
      }
    }

    ResetMe();
  }

    //with requiring that the map only store ep info for events that meets this z-vertex cut
    //there will be events with no stored ep info
  if (_mbdEpReco && (std::fabs(_mbd_vertex) < _vertex_cut))
  {
      ResetMe();
      MbdPmtContainer *mbdpmts = findNode::getClass<MbdPmtContainer>(topNode, "MbdPmtContainer");
      if (!mbdpmts)
      {
          std::cout << PHWHERE << "::ERROR - cannot find MbdPmtContainer" << std::endl;
          exit(-1);
      }
      
      MbdGeom *mbdgeom = findNode::getClass<MbdGeom>(topNode, "MbdGeom");
      if (!mbdgeom)
      {
          std::cout << PHWHERE << "::ERROR - cannot find MbdGeom" << std::endl;
          exit(-1);
      }
      
      
      if (mbdpmts)
      {
          if (Verbosity())
          {
              std::cout << "EventPlaneReco::process_event -  mbdpmts" << std::endl;
          }
          
          mbd_e_south = 0.;
          mbd_e_north = 0.;
          mbdQ = 0.;
          
          for (int ipmt = 0; ipmt < mbdpmts->get_npmt(); ipmt++)
          {
              float mbd_q = mbdpmts->get_pmt(ipmt)->get_q();
              mbdQ += mbd_q;
          }
          
          for (int ipmt = 0; ipmt < mbdpmts->get_npmt(); ipmt++)
          {
              float mbd_q = mbdpmts->get_pmt(ipmt)->get_q();
              float phi = mbdgeom->get_phi(ipmt);
              int arm = mbdgeom->get_arm(ipmt);
              
              if(mbdQ < _mbd_e) continue;
              
              if (arm == 0)
              {
                  mbd_e_south += mbd_q;
                  for (unsigned int order = 1; order < m_MaxOrder + 1; order++)
                  {
                      double Cosine = cos(phi * (double) order);
                      double Sine = sin(phi * (double) order);
                      south_q[order - 1][0] += mbd_q * Cosine;  // south Qn,x
                      south_q[order - 1][1] += mbd_q * Sine;    // south Qn,y
                  }
              }
              else if (arm == 1)
              {
                  mbd_e_north += mbd_q;
                  for (unsigned int order = 1; order < m_MaxOrder + 1; order++)
                  {
                      double Cosine = cos(phi * (double) order);
                      double Sine = sin(phi * (double) order);
                      north_q[order - 1][0] += mbd_q * Cosine;  // north Qn,x
                      north_q[order - 1][1] += mbd_q * Sine;    // north Qn,y
                  }
              }
            }
          }
      
          // Filled during first run
          for (unsigned int order = 1; order < m_MaxOrder + 1; order++)
          {
              // Fill recentering histograms by order
              tprof_mean_cos_south_mbd[order - 1]->Fill(mbd_e_south,south_q[order - 1][0]/mbd_e_south);
              tprof_mean_sin_south_mbd[order - 1]->Fill(mbd_e_south,south_q[order - 1][1]/mbd_e_south);
              tprof_mean_cos_north_mbd[order - 1]->Fill(mbd_e_north,north_q[order - 1][0]/mbd_e_north);
              tprof_mean_sin_north_mbd[order - 1]->Fill(mbd_e_north,north_q[order - 1][1]/mbd_e_north);
          }
          
          
          // Get recentering histograms and do recentering
          // Recentering: subtract Qn,x and Qn,y values averaged over all events
          // Recentering histogram should be available on second run
          for (unsigned int order = 1; order < m_MaxOrder + 1; order++)
          {
              if (tprof_mean_cos_south_mbd_input[order-1]) //check if recentering histograms exist
              {
                  // south
                  int bin_south = tprof_mean_cos_south_mbd_input[order - 1]->FindBin(mbd_e_south);
                  double event_ave_cos_south = tprof_mean_cos_south_mbd_input[order - 1]->GetBinContent(bin_south);
                  double event_ave_sin_south = tprof_mean_sin_south_mbd_input[order - 1]->GetBinContent(bin_south);
                  south_q_subtract[order - 1][0] = mbd_e_south*event_ave_cos_south;
                  south_q_subtract[order - 1][1] = mbd_e_south*event_ave_sin_south;
                  south_q[order - 1][0] -= south_q_subtract[order - 1][0];
                  south_q[order - 1][1] -= south_q_subtract[order - 1][1];
                  
                  // north
                  int bin_north = tprof_mean_cos_north_mbd_input[order - 1]->FindBin(mbd_e_north);
                  double event_ave_cos_north = tprof_mean_cos_north_mbd_input[order - 1]->GetBinContent(bin_north);
                  double event_ave_sin_north = tprof_mean_sin_north_mbd_input[order - 1]->GetBinContent(bin_north);
                  north_q_subtract[order - 1][0] = mbd_e_north*event_ave_cos_north;
                  north_q_subtract[order - 1][1] = mbd_e_north*event_ave_sin_north;
                  north_q[order - 1][0] -= north_q_subtract[order - 1][0];
                  north_q[order - 1][1] -= north_q_subtract[order - 1][1];
              }
          }
          
          

      // Higher order harmonics might still be present after recentering, so do event-by-event shifting of the planes
      // First get temp psi_n. Important to start with recentered planes, so ask for recentered histograms
      Eventplaneinfo *epinfo = new Eventplaneinfov1();
      for (unsigned int order = 1; order < m_MaxOrder + 1; order++)
      {
          if(tprof_mean_cos_south_mbd_input[order-1]) // if present, Qs are recentered
          {
              tmp_south_psi[order - 1] = epinfo->GetPsi(south_q[order - 1][0], south_q[order - 1][1], order);
              tmp_north_psi[order - 1] = epinfo->GetPsi(north_q[order - 1][0], north_q[order - 1][1], order);
          }
          else
          {
              tmp_south_psi[order - 1] = NAN;
              tmp_north_psi[order - 1] = NAN;
          }

      }
      

    // Filled during second run
    for (unsigned int order = 1; order < m_MaxOrder + 1; order++)
    {
        if( tprof_mean_cos_south_mbd_input[order-1]) // if present, Qs are recentered
        {
            // Fill shifting histograms by order and terms
            for ( int p = 0; p < _imax; p++ )
            {
                double terms = p+1.0;
                double tmp = (double)(order*terms);
                
                tprof_cos_south_mbd_shift[order - 1][p]->Fill(mbd_e_south, cos(tmp*tmp_south_psi[order - 1])); // south <cos(i*n*psi_n)>
                tprof_sin_south_mbd_shift[order - 1][p]->Fill(mbd_e_south, sin(tmp*tmp_south_psi[order - 1])); // south <sin(i*n*psi_n)>
                tprof_cos_north_mbd_shift[order - 1][p]->Fill(mbd_e_north, cos(tmp*tmp_north_psi[order - 1])); // north <cos(i*n*psi_n)>
                tprof_sin_north_mbd_shift[order - 1][p]->Fill(mbd_e_north, sin(tmp*tmp_north_psi[order - 1])); // north <sin(i*n*psi_n)>
                   
            }
        }
    }

     // Get shifting histograms and calculate shift
     // This was filled at the end of the second run
     // In third run, shifting histograms should be available
      for (unsigned int order = 1; order < m_MaxOrder + 1; order++)
      {
        for ( int p = 0; p < _imax; p++ )
        {
            if (tprof_cos_south_mbd_shift_input[order-1][p]) // check if shifting histograms exist
            {
                double terms = p+1.0;
                double tmp = (double)(order*terms);
                double prefactor = 2.0/terms;

                int bin_north = tprof_cos_north_mbd_shift_input[order - 1][p]->FindBin(mbd_e_north);
                int bin_south = tprof_cos_south_mbd_shift_input[order - 1][p]->FindBin(mbd_e_south);
 
                // Equation (6) of arxiv:nucl-ex/9805001
                // i = terms; n = order; i*n = tmp
                // (2 / i ) * <cos(i*n*psi_n)> * sin(i*n*psi_n) - <sin(i*n*psi_n)> * cos(i*n*psi_n)
                
                // north
                shift_north[order - 1] += prefactor*(tprof_cos_north_mbd_shift_input[order - 1][p]->GetBinContent(bin_north) * sin(tmp*tmp_north_psi[order-1]) - tprof_sin_north_mbd_shift_input[order - 1][p]->GetBinContent(bin_north) * cos(tmp*tmp_north_psi[order-1]));
                
                // south
                shift_south[order - 1] += prefactor*(tprof_cos_south_mbd_shift_input[order - 1][p]->GetBinContent(bin_south) * sin(tmp*tmp_south_psi[order-1]) - tprof_sin_south_mbd_shift_input[order - 1][p]->GetBinContent(bin_south) * cos(tmp*tmp_south_psi[order-1]));
                
            }
        }
      }
      
      // n * deltapsi_n = (2 / i ) * <cos(i*n*psi_n)> * sin(i*n*psi_n) - <sin(i*n*psi_n)> * cos(i*n*psi_n)
      // Divide out n
      for (unsigned int order = 1; order < m_MaxOrder + 1; order++)
      {
          shift_north[order - 1] /= order;
          shift_south[order - 1] /= order;
      }
      
      
      // Now add shift to psi_n to flatten it
      // Verify that this is done only in third run by asking for shifting hists
      for (unsigned int order = 1; order < m_MaxOrder + 1; order++)
      {
          if(tprof_cos_north_mbd_shift_input[0][0])
          {
            tmp_south_psi[order - 1] += shift_south[order - 1];
            tmp_north_psi[order - 1] += shift_north[order - 1];
          }
      }
      
      // Now enforce the range
      for (unsigned int order = 1; order < m_MaxOrder + 1; order++)
      {
          if(tprof_cos_north_mbd_shift_input[0][0])
          {
            double range = M_PI/(double)order;
            if(tmp_south_psi[order - 1] < -1.0*range) tmp_south_psi[order - 1] += 2.0*range;
            if(tmp_south_psi[order - 1] > range) tmp_south_psi[order - 1] -= 2.0*range;
            if(tmp_north_psi[order - 1] < -1.0*range) tmp_north_psi[order - 1] += 2.0*range;
            if(tmp_north_psi[order - 1] > range) tmp_north_psi[order - 1] -= 2.0*range;
          }
      }

    //-------------------------------- store all MBD eventplane information ---------------------------------//
    for (unsigned int order = 1; order < m_MaxOrder + 1; order++)
    {
         south_Qvec.emplace_back(south_q[order - 1][0], south_q[order - 1][1]);
         north_Qvec.emplace_back(north_q[order - 1][0], north_q[order - 1][1]); 
    }
                                                                                                    

    if (mbdpmts)
    { 
      Eventplaneinfo *mbds = new Eventplaneinfov1();
      mbds->set_qvector(south_Qvec);
      mbds->set_shifted_psi(tmp_south_psi);
      epmap->insert(mbds, EventplaneinfoMap::MBDS);

      Eventplaneinfo *mbdn = new Eventplaneinfov1();
      mbdn->set_qvector(north_Qvec);
      mbdn->set_shifted_psi(tmp_north_psi);
      epmap->insert(mbdn, EventplaneinfoMap::MBDN);

      if (Verbosity() > 1)
      {
        mbds->identify();
        mbdn->identify();
      }
    }
      
    ResetMe();
  }

  if (Verbosity())
  {
    epmap->identify();
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int EventPlaneReco::CreateNodes(PHCompositeNode *topNode)
{
  PHNodeIterator iter(topNode);

  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    std::cout << PHWHERE << "DST Node missing, doing nothing." << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  PHCompositeNode *globalNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "GLOBAL"));
  if (!globalNode)
  {
    globalNode = new PHCompositeNode("GLOBAL");
    dstNode->addNode(globalNode);
  }

  EventplaneinfoMap *eps = findNode::getClass<EventplaneinfoMap>(topNode, "EventplaneinfoMap");
  if (!eps)
  {
    eps = new EventplaneinfoMapv1();
    PHIODataNode<PHObject> *EpMapNode = new PHIODataNode<PHObject>(eps, "EventplaneinfoMap", "PHObject");
    globalNode->addNode(EpMapNode);
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

void EventPlaneReco::ResetMe()
{
  for (auto &vec : south_q)
  {
    std::fill(vec.begin(), vec.end(), 0.);
  }

  for (auto &vec : north_q)
  {
    std::fill(vec.begin(), vec.end(), 0.);
  }
    
  for (auto &vec : south_q_subtract)
  {
    std::fill(vec.begin(), vec.end(), 0.);
  }

  for (auto &vec : north_q_subtract)
  {
    std::fill(vec.begin(), vec.end(), 0.);
  }

  std::fill(shift_north.begin(), shift_north.end(), 0.);
  std::fill(shift_south.begin(), shift_south.end(), 0.);

  std::fill(tmp_south_psi.begin(), tmp_south_psi.end(), NAN);
  std::fill(tmp_north_psi.begin(), tmp_north_psi.end(), NAN);

  south_Qvec.clear();
  north_Qvec.clear();
    
}

int EventPlaneReco::End(PHCompositeNode * /*topNode*/)
{
  cdbhistosOut->WriteCDBHistos();
  delete cdbhistosOut;
  std::cout << " EventPlaneReco::End() " << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}
