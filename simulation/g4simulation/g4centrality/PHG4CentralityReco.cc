#include "PHG4CentralityReco.h"

#include <centrality/CentralityInfo.h>    // for CentralityInfo, CentralityI...
#include <centrality/CentralityInfov1.h>  // for CentralityInfov1

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>  // for PHObject
#include <phool/getClass.h>
#include <phool/phool.h>  // for PHWHERE

#include <ffaobjects/EventHeaderv1.h>
#include <g4main/PHG4Hit.h>
#include <g4main/PHG4HitContainer.h>

#include <iostream>  // for operator<<, basic_ostream
#include <map>       // for _Rb_tree_const_iterator
#include <sstream>
#include <stdexcept>  // for runtime_error
#include <utility>    // for pair

#include "EpFinder.h"

PHG4CentralityReco::PHG4CentralityReco(const std::string &name)
  : SubsysReco(name)
  , _centrality_calibration_params(name)
{
    EPD_EpFinderN = NULL; //north Epd
    EPD_EpFinderS = NULL;
    EPD_EpFinderN_Trunc = NULL; //north Epd with truncation applied
    EPD_EpFinderS_Trunc = NULL;
}

int PHG4CentralityReco::InitRun(PHCompositeNode *topNode)
{
  if (Verbosity() >= 1)
    std::cout << "PHG4CentralityReco::InitRun : enter " << std::endl;

  try
  {
    CreateNode(topNode);
  }
  catch (std::exception &e)
  {
    std::cout << PHWHERE << ": " << e.what() << std::endl;
    throw;
  }

  auto bhits = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_BBC");
  if (!bhits)
    std::cout << "PHG4CentralityReco::InitRun : cannot find G4HIT_BBC, will not use MBD centrality" << std::endl;

  auto ehits = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_EPD");
  if (!ehits)
    std::cout << "PHG4CentralityReco::InitRun : cannot find G4HIT_EPD, will not use sEPD centrality" << std::endl;

  if (_centrality_calibration_params.exist_string_param("description"))
  {
    if (Verbosity() >= 1)
    {
      std::cout << "PHG4CentralityReco::InitRun : Centrality calibration description : " << std::endl
                << "    ";
      std::cout << _centrality_calibration_params.get_string_param("description") << std::endl;
    }
    // search for possible centile definitions
    for (int n = 0; n < 101; n++)
    {
      std::ostringstream s1;
      s1 << "epd_centile_" << n;
      if (_centrality_calibration_params.exist_double_param(s1.str().c_str()))
      {
        _cent_cal_epd[_centrality_calibration_params.get_double_param(s1.str().c_str())] = n;
        if (Verbosity() >= 2)
          std::cout << "PHG4CentralityReco::InitRun : sEPD centrality calibration, centile " << n << "% is " << _centrality_calibration_params.get_double_param(s1.str().c_str()) << std::endl;
      }
    }
    for (int n = 0; n < 101; n++)
    {
      std::ostringstream s2;
      s2 << "mbd_centile_" << n;
      if (_centrality_calibration_params.exist_double_param(s2.str().c_str()))
      {
        _cent_cal_mbd[_centrality_calibration_params.get_double_param(s2.str().c_str())] = n;
        if (Verbosity() >= 2)
          std::cout << "PHG4CentralityReco::InitRun : MBD centrality calibration, centile " << n << "% is " << _centrality_calibration_params.get_double_param(s2.str().c_str()) << std::endl;
      }
    }
    for (int n = 0; n < 101; n++)
    {
      std::ostringstream s3;
      s3 << "bimp_centile_" << n;
      if (_centrality_calibration_params.exist_double_param(s3.str().c_str()))
      {
        _cent_cal_bimp[_centrality_calibration_params.get_double_param(s3.str().c_str())] = n;
        if (Verbosity() >= 2)
          std::cout << "PHG4CentralityReco::InitRun : b (impact parameter) centrality calibration, centile " << n << "% is " << _centrality_calibration_params.get_double_param(s3.str().c_str()) << std::endl;
      }
    }
  }
  else
  {
    std::cout << "PHG4CentralityReco::InitRun : no centrality calibration found!" << std::endl;
  }
    
    //binning (x,y); (eta, phi) : (16, 24)
    EPD_EpFinderN = new EpFinder(1,"EPD_OUTPUT.root", "EPD_INPUT.root", 16, 24);
    EPD_EpFinderS = new EpFinder(1, "EPDS_OUTPUT.root", "EPDS_INPUT.root", 16, 24);
    EPD_EpFinderN_Trunc = new EpFinder(1, "EPD_Trunc_OUTPUT.root", "EPD_Trunc_INPUT.root", 16, 24);
    EPD_EpFinderS_Trunc = new EpFinder(1, "EPDS_Trunc_OUTPUT.root", "EPDS_Trunc_INPUT.root", 16, 24);

  return Fun4AllReturnCodes::EVENT_OK;
}

int PHG4CentralityReco::GetPhiBin(float tphi, float numPhiDivisions)
{

  float sphi = ((2.0*M_PI)/numPhiDivisions);

  if(tphi>=(2.0*M_PI))
  {
    tphi -= (2.0*M_PI);
  }
  else if(tphi<0.0)
  {
    tphi += (2.0*M_PI);
  }
       
  return (int)(tphi/sphi);

}

float PHG4CentralityReco::GetMeanPhi(int iphi, float numPhiDivisions)
{

  float sphi = ((2.0*M_PI)/numPhiDivisions);
  float tphi = ((float)iphi + 0.5)*sphi;

  if(tphi>=(2.0*M_PI))
  {
    tphi -= (2.0*M_PI);
  }

  else if(tphi<0.0)
  {
    tphi += (2.0*M_PI);
  }
       
  return tphi;

}
 
int PHG4CentralityReco::GetEtaBin(float teta, float eta_low, float eta_high, float numEtaDivisions)
{

   float seta = fabs((eta_high-eta_low)/numEtaDivisions);
   int ieta = ((teta-eta_low)/seta);
   return fabs(ieta);
}


int PHG4CentralityReco::process_event(PHCompositeNode *topNode)
{
  
  _bimp = 101;
  auto event_header = findNode::getClass<EventHeaderv1>(topNode, "EventHeader");
  if (event_header)
  {
    _bimp = event_header->get_floatval("bimp");

    if (Verbosity() >= 5)
      std::cout << "PHG4CentralityReco::process_event : Hijing impact parameter b = " << _bimp << std::endl;
  }
  else
  {
    if (Verbosity() >= 5)
      std::cout << "PHG4CentralityReco::process_event : No Hijing impact parameter info, setting b = 101" << std::endl;
  }

  _mbd_N = 0;
  _mbd_S = 0;
  _mbd_NS = 0;

  auto bhits = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_BBC");

  if (bhits)
  {
    auto brange = bhits->getHits();
    for (auto it = brange.first; it != brange.second; ++it)
    {
      if ((it->second->get_t(0) > -50) && (it->second->get_t(1) < 50))
      {
        _mbd_NS += it->second->get_edep();
        int id = it->second->get_layer();
        if ((id & 0x40) == 0)
          _mbd_N += it->second->get_edep();
        else
          _mbd_S += it->second->get_edep();
      }
    }

    if (Verbosity() >= 5)
      std::cout << "PHG4CentralityReco::process_event : MBD Sum Charge N / S / N+S = " << _mbd_N << " / " << _mbd_S << " / " << _mbd_NS << std::endl;
  }
  else
  {
    if (Verbosity() >= 5)
      std::cout << "PHG4CentralityReco::process_event : No MBD info, setting all Sum Charges = 0" << std::endl;
  }

  _epd_N = 0;
  _epd_S = 0;
  _epd_NS = 0;

  auto ehits = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_EPD");

  if (ehits)
  {
    auto erange = ehits->getHits();
    for (auto it = erange.first; it != erange.second; ++it)
      if ((it->second->get_t(0) > -50) && (it->second->get_t(1) < 50))
      {
        _epd_NS += it->second->get_edep();
        int id = it->second->get_scint_id();
        if ((id & 0x200) == 0)
          _epd_N += it->second->get_edep();
        else
          _epd_S += it->second->get_edep();
      }

    if (Verbosity() >= 5)
      std::cout << "PHG4CentralityReco::process_event : sEPD Sum Energy N / S / N+S = " << _epd_N << " / " << _epd_S << " / " << _epd_NS << std::endl;
  }
  else
  {
    if (Verbosity() >= 5)
      std::cout << "PHG4CentralityReco::process_event : No sEPD info, setting all Sum Energies = 0" << std::endl;
  }

  if (Verbosity() >= 1)
    std::cout << "PHG4CentralityReco::process_event : summary MBD (N, S, N+S) = (" << _mbd_N << ", " << _mbd_S << ", " << _mbd_NS << "), sEPD (N, S, N+S) = (" << _epd_N << ", " << _epd_S << ", " << _epd_NS << ")" << std::endl;

  if (_do_centrality_calibration)
  {
    // sEPD centrality
    float low_epd_val = -10000;
    float high_epd_val = 10000;
    int low_epd_centile = -1;
    int high_epd_centile = -1;

    for (std::map<float, int>::iterator it = _cent_cal_epd.begin(); it != _cent_cal_epd.end(); ++it)
    {
      float signal = it->first;
      int cent = it->second;

      if (signal < _epd_NS && signal > low_epd_val)
      {
        low_epd_val = signal;
        low_epd_centile = cent;
      }
      if (signal > _epd_NS && signal < high_epd_val)
      {
        high_epd_val = signal;
        high_epd_centile = cent;
      }

    }  // close iterate through sEPD cuts

    if (low_epd_centile >= 0 && high_epd_centile >= 0)
    {
      _epd_cent = (low_epd_centile + high_epd_centile) / 2.0;
      if (Verbosity() >= 10)
        std::cout << "PHG4CentralityReco::process_event : lower EPD value is " << low_epd_val << " (" << low_epd_centile << "%), higher is " << high_epd_val << " (" << high_epd_centile << "%), assigning " << _epd_cent << "%" << std::endl;
    }
    else
    {
      _epd_cent = 101;
      if (Verbosity() >= 5)
        std::cout << "PHG4CentralityReco::process_event : not able to map EPD value to a centrality. debug info = " << low_epd_val << "/" << low_epd_centile << "/" << high_epd_val << "/" << high_epd_centile << std::endl;
    }

    // MBD centrality
    float low_mbd_val = -10000;
    float high_mbd_val = 10000;
    int low_mbd_centile = -1;
    int high_mbd_centile = -1;

    for (std::map<float, int>::iterator it = _cent_cal_mbd.begin(); it != _cent_cal_mbd.end(); ++it)
    {
      float signal = it->first;
      int cent = it->second;

      if (signal < _mbd_NS && signal > low_mbd_val)
      {
        low_mbd_val = signal;
        low_mbd_centile = cent;
      }
      if (signal > _mbd_NS && signal < high_mbd_val)
      {
        high_mbd_val = signal;
        high_mbd_centile = cent;
      }

    }  // close iterate through MBD cuts

    if (low_mbd_centile >= 0 && high_mbd_centile >= 0)
    {
      _mbd_cent = (low_mbd_centile + high_mbd_centile) / 2.0;
      if (Verbosity() >= 10)
        std::cout << "PHG4CentralityReco::process_event : lower MBD value is " << low_mbd_val << " (" << low_mbd_centile << "%), higher is " << high_mbd_val << " (" << high_mbd_centile << "%), assigning " << _mbd_cent << "%" << std::endl;
    }
    else
    {
      _mbd_cent = 101;
      if (Verbosity() >= 5)
        std::cout << "PHG4CentralityReco::process_event : not able to map MBD value to a centrality. debug info = " << low_mbd_val << "/" << low_mbd_centile << "/" << high_mbd_val << "/" << high_mbd_centile << std::endl;
    }

    // b (impact parameter) centrality
    float low_bimp_val = -1;
    float high_bimp_val = 10000;
    int low_bimp_centile = -1;
    int high_bimp_centile = -1;

    for (std::map<float, int>::iterator it = _cent_cal_bimp.begin(); it != _cent_cal_bimp.end(); ++it)
    {
      float signal = it->first;
      int cent = it->second;

      if (signal < _bimp && signal > low_bimp_val)
      {
        low_bimp_val = signal;
        low_bimp_centile = cent;
      }
      if (signal > _bimp && signal < high_bimp_val)
      {
        high_bimp_val = signal;
        high_bimp_centile = cent;
      }

    }  // close iterate through bimp cuts

    if (low_bimp_centile >= 0 && high_bimp_centile >= 0)
    {
      _bimp_cent = (low_bimp_centile + high_bimp_centile) / 2.0;
      if (Verbosity() >= 10)
        std::cout << "PHG4CentralityReco::process_event : lower b value is " << low_bimp_val << " (" << low_bimp_centile << "%), higher is " << high_bimp_val << " (" << high_bimp_centile << "%), assigning " << _bimp_cent << "%" << std::endl;
    }
    else
    {
      _bimp_cent = 101;
      if (Verbosity() >= 5)
        std::cout << "PHG4CentralityReco::process_event : not able to map b value to a centrality. debug info = " << low_bimp_val << "/" << low_bimp_centile << "/" << high_bimp_val << "/" << high_bimp_centile << std::endl;
    }

  }  // close centrality calibration

GetEventPlanes(topNode);
FillNode(topNode);

  return Fun4AllReturnCodes::EVENT_OK;
}


void PHG4CentralityReco::GetEventPlanes(PHCompositeNode *topNode)
{
    
    //Read the detector geometry and set up arrays for phi weighting
          
      static bool first = true;

      if(first){


       for(int i=0; i<24; i++)
       {
          Nepd_phi_list[i].clear();
          Sepd_phi_list[i].clear();
          TNepd_phi_list[i].clear();
          TSepd_phi_list[i].clear();
        }


    //generate a list of all towers in the same phi range
    //all eta bins are at the same phi

    for(int i=0; i<24; i++)
     {
        for(int j=0; j<16; j++)
        {
          std::pair<int,int> newPair(j,i);
          Nepd_phi_list[i].push_back(newPair);
          Sepd_phi_list[i].push_back(newPair);
          TNepd_phi_list[i].push_back(newPair);
          TSepd_phi_list[i].push_back(newPair);
        }
     }

        first = false;
    }

    int thisphibin = -1;
    int thisetabin = -1;
    int iarm = -1;
    float EPDEnergy[2][16][24] = {{{0.0}}};
    float EPDEnergyTrunc[2][16][24] = {{{0.0}}};
    float thisMip = 0.;
    
    _epd_N_ep = 0.;
    _epd_S_ep = 0.;
    _epd_t_N_ep = 0.;
    _epd_t_S_ep = 0.;
    
    
    //array for Nmip cut off (per centrality, per ring)
    float Nmipcutoff[10][16] =
    {
           { 8., 8., 8., 8., 8., 8., 8., 8., 8., 10., 10., 10., 10., 10., 10., 10., },
           { 7., 7., 7., 7., 7., 7., 7., 7., 7., 9., 9., 9., 9., 10., 10.,10., },
           { 5.,  5.,  5.,  5.,  5.,  5.,  5.,  5., 5., 7., 7., 7., 7., 7., 7., 6., },
           { 4.,  4.,  4.,  4.,  4.,  4.,  4.,  4., 4., 5., 5., 5., 5., 5., 5., 8., },
           { 3.,  3.,  3.,  3.,  3.,  3.,  3.,  3., 3., 4., 4., 4., 4., 4., 4., 4., },
           { 2.,  2.,  2.,  2.,  2.,  2.,  2.,  2., 2., 3., 3., 3., 3., 2., 2., 2., },
           { 2.,  2.,  2.,  2.,  2.,  2.,  2.,  2., 2., 2., 2., 2., 2., 2., 2., 2., },
           { 2.,  2.,  2.,  2.,  2.,  2.,  2.,  2., 2., 2., 2., 2., 2., 2., 2., 2., },
           { 2.,  2.,  2.,  2.,  2.,  2.,  2.,  2., 2., 2., 2., 2., 2., 2., 2., 2., },
           { 2.,  2.,  2.,  2.,  2.,  2.,  2.,  2., 2., 2., 2., 2., 2., 2., 2., 2., },
    };

  
  //get centrality info
  int _b = -1;
  
  if( (_epd_cent>=0) && (_epd_cent <10)) _b = 0;
  if( (_epd_cent>=10) && (_epd_cent <20)) _b = 1;
  if( (_epd_cent>=20) && (_epd_cent <30)) _b = 2;
  if( (_epd_cent>=30) && (_epd_cent <40)) _b = 3;
  if( (_epd_cent>=40) && (_epd_cent <50)) _b = 4;
  if( (_epd_cent>=50) && (_epd_cent <60)) _b = 5;
  if( (_epd_cent>=60) && (_epd_cent <70)) _b = 6;
  if( (_epd_cent>=70) && (_epd_cent <80)) _b = 7;
  if( (_epd_cent>=80) && (_epd_cent <90)) _b = 8;
  if( (_epd_cent>=90) && (_epd_cent <100)) _b = 9;

  if(_b < 0) _b = 9;
    
    auto _epd_hit_container = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_EPD");
    if (!_epd_hit_container)
      std::cout << "PHG4CentralityReco::InitRun : cannot find G4HIT_EPD" << std::endl;

   PHG4HitContainer::ConstRange e_range = _epd_hit_container->getHits();
   for(PHG4HitContainer::ConstIterator e_itr = e_range.first; e_itr != e_range.second; e_itr++)
   {

    PHG4Hit *this_epdhit = e_itr->second;
    if(!this_epdhit){ continue;}

    if(this_epdhit->get_light_yield()>0.0)
    {

       TVector3 hitPos(this_epdhit->get_avg_x(), this_epdhit->get_avg_y(), this_epdhit->get_avg_z());

       if(this_epdhit->get_z(0) > 0) {iarm = 1; }
       else{iarm = 0;}

       if(TMath::Abs(hitPos.Eta()) < 2.01) continue;
       if(TMath::Abs(hitPos.Eta()) > 4.96) continue;

      //Tile segmentation in eta/phi

       if(iarm==1)
       {

       if((hitPos.Eta() >= 2.01) && (hitPos.Eta() < 4.29))
       {

       thisphibin = GetPhiBin(hitPos.Phi(), 24.0);

       if((hitPos.Eta() >= 2.01) && (hitPos.Eta() < 2.07))
       {thisetabin = 0;}

       if((hitPos.Eta() >= 2.07) && (hitPos.Eta() < 2.14))
       {thisetabin = 1;}

       if((hitPos.Eta() >= 2.14) && (hitPos.Eta() < 2.21))
       {thisetabin = 2;}

       if((hitPos.Eta() >= 2.21) && (hitPos.Eta() < 2.28))
       {thisetabin = 3;}

       if((hitPos.Eta() >= 2.28) && (hitPos.Eta() < 2.37))
       {thisetabin = 4;}

       if((hitPos.Eta() >= 2.37) && (hitPos.Eta() < 2.46))
       {thisetabin = 5;}

       if((hitPos.Eta() >= 2.46) && (hitPos.Eta() < 2.56))
       {thisetabin = 6;}

       if((hitPos.Eta() >= 2.56) && (hitPos.Eta() < 2.68))
       {thisetabin = 7;}

       if((hitPos.Eta() >= 2.68) && (hitPos.Eta() < 2.81))
       {thisetabin = 8;}

       if((hitPos.Eta() >= 2.81) && (hitPos.Eta() < 2.95))
       {thisetabin = 9;}

       if((hitPos.Eta() >= 2.95) && (hitPos.Eta() < 3.13))
       {thisetabin = 10;}

       if((hitPos.Eta() >= 3.13) && (hitPos.Eta() < 3.34))
       {thisetabin = 11;}

       if((hitPos.Eta() >= 3.34) && (hitPos.Eta() < 3.61))
       {thisetabin = 12;}

       if((hitPos.Eta() >= 3.61) && (hitPos.Eta() < 3.9))
       {thisetabin = 13;}

       if((hitPos.Eta() >= 3.9) && (hitPos.Eta() < 4.29))
       {thisetabin = 14;}

       }
       //innermost ring
       if ((hitPos.Eta() >= 4.29) && (hitPos.Eta() < 4.96))
       {
        thisetabin = 15;
        thisphibin = GetPhiBin(hitPos.Phi(), 12.0);
       }

       }//north


       if(iarm==0)
       {

       if((hitPos.Eta() <= -2.01) && (hitPos.Eta() > -4.29))
       {

       thisphibin = GetPhiBin(hitPos.Phi(), 24.0);

       if((hitPos.Eta() <= -2.01) && (hitPos.Eta() > -2.07))
       {thisetabin = 0;}

       if((hitPos.Eta() <= -2.07) && (hitPos.Eta() > -2.14))
       {thisetabin = 1;}

       if((hitPos.Eta() <= -2.14) && (hitPos.Eta() > -2.21))
       {thisetabin = 2;}

       if((hitPos.Eta() <= - 2.21) && (hitPos.Eta() > -2.28))
       {thisetabin = 3;}

       if((hitPos.Eta() <= -2.28) && (hitPos.Eta() > -2.37))
       {thisetabin = 4;}

       if((hitPos.Eta() <= -2.37) && (hitPos.Eta() > -2.46))
       {thisetabin = 5;}

       if((hitPos.Eta() <= -2.46) && (hitPos.Eta() > -2.56))
       {thisetabin = 6;}

       if((hitPos.Eta() <= -2.56) && (hitPos.Eta() > -2.68))
       {thisetabin = 7;}

       if((hitPos.Eta() <= -2.68) && (hitPos.Eta() > -2.81))
       {thisetabin = 8;}

       if((hitPos.Eta() <= -2.81) && (hitPos.Eta() > -2.95))
       {thisetabin = 9;}

       if((hitPos.Eta() <= -2.95) && (hitPos.Eta() > -3.13))
       {thisetabin = 10;}

       if((hitPos.Eta() <= -3.13) && (hitPos.Eta() > -3.34))
       {thisetabin = 11;}

       if((hitPos.Eta() <= -3.34) && (hitPos.Eta() > -3.61))
       {thisetabin = 12;}

       if((hitPos.Eta() <= -3.61) && (hitPos.Eta() > -3.9))
       {thisetabin = 13;}

       if((hitPos.Eta() <= -3.9) && (hitPos.Eta() > -4.29))
       {thisetabin = 14;}
       }

       //innermost ring
       if ((hitPos.Eta() <= -4.29) && (hitPos.Eta() > -4.96))
       {
       thisphibin = GetPhiBin(hitPos.Phi(), 12.0);
       thisetabin = 15;
       }
      }//south

        //normalized by MPV
       thisMip = this_epdhit->get_light_yield()/2.05924e-6;


       if((thisphibin >= 0) && (thisetabin >= 0)) EPDEnergy[iarm][thisetabin][thisphibin] += this_epdhit->get_light_yield();
       if((thisphibin >= 0) && (thisetabin >= 0)) EPDEnergyTrunc[iarm][thisetabin][thisphibin] += thisMip;

    }//energy
   }//end loop over hits
    
    
    // -------------------------------------
    // Run the EPD Event Plane Finder
    // -------------------------------------


    std::vector<EpHit> Nepdsimhits;
    Nepdsimhits.clear();
  
    
    //NORTH
    /*******************************************************************/
    
         for(int i=0; i<16; i++){
            for(int j=0; j<24; j++){
    
              if(EPDEnergy[1][i][j]>0.0){
              float meanPhi = GetMeanPhi(j, 24.0);
              if(i==15) meanPhi = GetMeanPhi(j, 12.0);
    
              EpHit newepdHit;
              newepdHit.phi = meanPhi;
              newepdHit.nMip = EPDEnergy[1][i][j];
              newepdHit.iy = j;
              newepdHit.ix = i;
              newepdHit.samePhi = &Nepd_phi_list[newepdHit.iy];
              Nepdsimhits.push_back(newepdHit);}
        }
          }

        
     EpInfo EPDN_EpResult = EPD_EpFinderN->Results(&Nepdsimhits,0);
     _epd_N_ep = EPDN_EpResult.RawPsi(2);
    
    
    std::vector<EpHit> Sepdsimhits;
    Sepdsimhits.clear();
    
    //SOUTH
    /*******************************************************************/

         for(int i=0; i<16; i++){
            for(int j=0; j<24; j++){
             if(EPDEnergy[0][i][j]>0.0){
             
              float meanPhi = GetMeanPhi(j, 24.0);
               if(i==15) meanPhi = GetMeanPhi(j, 12.0);
              
              EpHit newepdHit;
              newepdHit.phi = meanPhi;
              newepdHit.nMip = EPDEnergy[0][i][j];
              newepdHit.iy = j;
              newepdHit.ix = i;
              newepdHit.samePhi = &Sepd_phi_list[newepdHit.iy];
              Sepdsimhits.push_back(newepdHit);}
        }
      }
           
       EpInfo SEPD_EpResult = EPD_EpFinderS->Results(&Sepdsimhits,0);
       _epd_S_ep = SEPD_EpResult.RawPsi(2);

    
    
    std::vector<EpHit> TNepdsimhits;
    TNepdsimhits.clear();
    
    //NORTH-MIP CALIB
    /*******************************************************************/
    
        float thisCalE = 0.0;
         for(int i=0; i<16; i++){
            for(int j=0; j<24; j++){
          
             if(EPDEnergyTrunc[1][i][j]>0.0){
     
              float meanPhi = GetMeanPhi(j, 24.0);
              if(i==15) meanPhi = GetMeanPhi(j, 12.0);
              
              if (EPDEnergyTrunc[1][i][j]<0.2) continue;
              thisCalE = (EPDEnergyTrunc[1][i][j]<Nmipcutoff[_b][i])?EPDEnergyTrunc[1][i][j]:Nmipcutoff[_b][i];
              
              EpHit newepdHit;
              newepdHit.phi = meanPhi;
              newepdHit.nMip = thisCalE;
              newepdHit.iy = j;
              newepdHit.ix = i;
              newepdHit.samePhi = &TNepd_phi_list[newepdHit.iy];
              TNepdsimhits.push_back(newepdHit);}
        }
      }
     
      EpInfo EPDN_Trunc_EpResult = EPD_EpFinderN_Trunc->Results(&TNepdsimhits,0);
      _epd_t_N_ep = EPDN_Trunc_EpResult.RawPsi(2);


    std::vector<EpHit> TSepdsimhits;
    TSepdsimhits.clear();
    
    //SOUTH-MIP CALIB
    /*******************************************************************/

    for(int i=0; i<16; i++){
       for(int j=0; j<24; j++){
        if(EPDEnergyTrunc[0][i][j]>0.0){
        
         float meanPhi = GetMeanPhi(j, 24.0);
          if(i==15) meanPhi = GetMeanPhi(j, 12.0);


         if (EPDEnergyTrunc[0][i][j]<0.2) continue;
         thisCalE = (EPDEnergyTrunc[0][i][j]<Nmipcutoff[_b][i])?EPDEnergyTrunc[0][i][j]:Nmipcutoff[_b][i];
        
         EpHit newepdHit;
         newepdHit.phi = meanPhi;
         newepdHit.nMip = thisCalE;
         newepdHit.iy = j;
         newepdHit.ix = i;
         newepdHit.samePhi = &TSepd_phi_list[newepdHit.iy];
         TSepdsimhits.push_back(newepdHit);}
   }
}

  EpInfo SEPD_Trunc_EpResult = EPD_EpFinderS_Trunc->Results(&TSepdsimhits,0);
 _epd_t_S_ep = SEPD_Trunc_EpResult.RawPsi(2);

  return;

}


int PHG4CentralityReco::End(PHCompositeNode * /*topNode*/)
{
    EPD_EpFinderN->Finish();
    EPD_EpFinderS->Finish();
    EPD_EpFinderN_Trunc->Finish();
    EPD_EpFinderS_Trunc->Finish();

    delete EPD_EpFinderN;
    delete EPD_EpFinderS;
    delete EPD_EpFinderN_Trunc;
    delete EPD_EpFinderS_Trunc;
    
  return Fun4AllReturnCodes::EVENT_OK;
}

void PHG4CentralityReco::FillNode(PHCompositeNode *topNode)
{
 
  CentralityInfo *cent = findNode::getClass<CentralityInfo>(topNode, "CentralityInfo");
  if (!cent)
  {
    std::cout << " ERROR -- can't find CentralityInfo node after it should have been created" << std::endl;
    return;
  }
  else
  {
    cent->set_quantity(CentralityInfo::PROP::mbd_N, _mbd_N);
    cent->set_quantity(CentralityInfo::PROP::mbd_S, _mbd_S);
    cent->set_quantity(CentralityInfo::PROP::mbd_NS, _mbd_NS);
    cent->set_quantity(CentralityInfo::PROP::epd_N, _epd_N);
    cent->set_quantity(CentralityInfo::PROP::epd_S, _epd_S);
    cent->set_quantity(CentralityInfo::PROP::epd_NS, _epd_NS);
    cent->set_quantity(CentralityInfo::PROP::bimp, _bimp);
      
    cent->set_quantity(CentralityInfo::PROP::epd_N_EP, _epd_N_ep);
    cent->set_quantity(CentralityInfo::PROP::epd_S_EP, _epd_S_ep);
    cent->set_quantity(CentralityInfo::PROP::epd_N_trunc_EP, _epd_t_N_ep);
    cent->set_quantity(CentralityInfo::PROP::epd_S_trunc_EP, _epd_t_S_ep);

    cent->set_centile(CentralityInfo::PROP::epd_NS, _epd_cent);
    cent->set_centile(CentralityInfo::PROP::mbd_NS, _mbd_cent);
    cent->set_centile(CentralityInfo::PROP::bimp, _bimp_cent);

    cent->set_centile(CentralityInfo::PROP::epd_NS, _epd_cent);
    cent->set_centile(CentralityInfo::PROP::mbd_NS, _mbd_cent);
    cent->set_centile(CentralityInfo::PROP::bimp, _bimp_cent);
    cent->set_centile(CentralityInfo::PROP::bimp, _bimp_cent);

  }
}


void PHG4CentralityReco::CreateNode(PHCompositeNode *topNode)
{
  PHNodeIterator iter(topNode);

  PHCompositeNode *dstNode = static_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    std::cerr << PHWHERE << "DST Node missing, doing nothing." << std::endl;
    throw std::runtime_error("Failed to find DST node in PHG4CentralityReco::CreateNode");
  }

  PHNodeIterator dstiter(dstNode);
  PHCompositeNode *DetNode = dynamic_cast<PHCompositeNode *>(dstiter.findFirst("PHCompositeNode", "GLOBAL"));
  if (!DetNode)
  {
    DetNode = new PHCompositeNode("GLOBAL");
    dstNode->addNode(DetNode);
  }

  CentralityInfo *cent = new CentralityInfov1();
  PHIODataNode<PHObject> *centNode = new PHIODataNode<PHObject>(cent, "CentralityInfo", "PHObject");
  DetNode->addNode(centNode);
 
 
}
