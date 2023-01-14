#include "fillSpaceChargeMaps.h"

#include "Shifter.h"

#include <g4main/PHG4Hit.h>
#include <g4main/PHG4HitContainer.h>

#include <fun4all/Fun4AllHistoManager.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>  // for SubsysReco

#include <phool/getClass.h>

#include <TAxis.h>  // for TAxis
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TTree.h>
#include <TVector3.h>

#include <algorithm>  // for max
//#include <cmath>      // for sin, asin, cos, floor, M_PI
#include <cstdio>     // for sprintf, printf
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <utility>  // for pair
#include <vector>

//____________________________________________________________________________..
fillSpaceChargeMaps::fillSpaceChargeMaps(const std::string &name, const std::string &filename)
  : SubsysReco(name)
  , _filename(filename)
{
  std::cout << "fillSpaceChargeMaps::fillSpaceChargeMaps(const std::string &name) Calling ctor" << std::endl;
}

//____________________________________________________________________________..
fillSpaceChargeMaps::~fillSpaceChargeMaps()
{
  //std::cout << "fillSpaceChargeMaps::~fillSpaceChargeMaps() Calling dtor" << std::endl;
  // delete whatever is created
  delete hm;
}

//____________________________________________________________________________..
int fillSpaceChargeMaps::Init(PHCompositeNode * /*topNode*/)
{
  int nz = 72;
  double z_rdo = 108 * cm;
  outfile = new TFile(_filename.c_str(), "RECREATE");

  hm = new Fun4AllHistoManager("HITHIST");
  const int r_bins_N = 66;  //51;
  double r_bins[r_bins_N + 1] = {217.83,
                                 223.83, 229.83, 235.83, 241.83, 247.83, 253.83, 259.83, 265.83, 271.83, 277.83, 283.83, 289.83, 295.83, 301.83, 306.83,
                                 311.05, 317.92, 323.31, 329.27, 334.63, 340.59, 345.95, 351.91, 357.27, 363.23, 368.59, 374.55, 379.91, 385.87, 391.23, 397.19, 402.49,
                                 411.53, 421.70, 431.90, 442.11, 452.32, 462.52, 472.73, 482.94, 493.14, 503.35, 513.56, 523.76, 533.97, 544.18, 554.39, 564.59, 574.76,
                                 583.67, 594.59, 605.57, 616.54, 627.51, 638.48, 649.45, 660.42, 671.39, 682.36, 693.33, 704.30, 715.27, 726.24, 737.21, 748.18, 759.11};
  double r_bins_edges[r_bins_N + 3] = {217.83-2,217.83,
                                 223.83, 229.83, 235.83, 241.83, 247.83, 253.83, 259.83, 265.83, 271.83, 277.83, 283.83, 289.83, 295.83, 301.83, 306.83,
                                 311.05, 317.92, 323.31, 329.27, 334.63, 340.59, 345.95, 351.91, 357.27, 363.23, 368.59, 374.55, 379.91, 385.87, 391.23, 397.19, 402.49,
                                 411.53, 421.70, 431.90, 442.11, 452.32, 462.52, 472.73, 482.94, 493.14, 503.35, 513.56, 523.76, 533.97, 544.18, 554.39, 564.59, 574.76,
                                 583.67, 594.59, 605.57, 616.54, 627.51, 638.48, 649.45, 660.42, 671.39, 682.36, 693.33, 704.30, 715.27, 726.24, 737.21, 748.18, 759.11, 759.11+2};

  const int nphi = 205;

  double phi_bins[nphi + 1] = { 0., 6.3083-2 * M_PI, 6.3401-2 * M_PI, 6.372-2 * M_PI, 6.4039-2 * M_PI, 6.4358-2 * M_PI, 6.4676-2 * M_PI, 6.4995-2 * M_PI, 6.5314-2 * M_PI, 
                                0.2618, 0.2937, 0.3256, 0.3574, 0.3893, 0.4212, 0.453, 0.4849, 0.5168, 0.5487, 0.5805, 0.6124, 0.6443, 0.6762, 0.7081, 
                                0.7399, 0.7718, 0.7854, 0.8173, 0.8491, 0.881, 0.9129, 0.9448, 0.9767, 1.0085, 1.0404, 1.0723, 1.1041, 1.136, 1.1679, 
                                1.1998, 1.2317, 1.2635, 1.2954, 1.309, 1.3409, 1.3727, 1.4046, 1.4365, 1.4684, 1.5002, 1.5321, 1.564, 1.5959, 1.6277, 
                                1.6596, 1.6915, 1.7234, 1.7552, 1.7871, 1.819, 1.8326, 1.8645, 1.8963, 1.9282, 1.9601, 1.992, 2.0238, 2.0557, 2.0876, 
                                2.1195, 2.1513, 2.1832, 2.2151, 2.247, 2.2788, 2.3107, 2.3426, 2.3562, 2.3881, 2.42, 2.4518, 2.4837, 2.5156, 2.5474, 
                                2.5793, 2.6112, 2.6431, 2.6749, 2.7068, 2.7387, 2.7706, 2.8024, 2.8343, 2.8662, 2.8798, 2.9117, 2.9436, 2.9754, 3.0073, 
                                3.0392, 3.0711, 3.1029, 3.1348, 3.1667, 3.1986, 3.2304, 3.2623, 3.2942, 3.326, 3.3579, 3.3898, 3.4034, 3.4353, 3.4671, 
                                3.499, 3.5309, 3.5628, 3.5946, 3.6265, 3.6584, 3.6903, 3.7221, 3.754, 3.7859, 3.8178, 3.8496, 3.8815, 3.9134, 3.927, 
                                3.9589, 3.9907, 4.0226, 4.0545, 4.0864, 4.1182, 4.1501, 4.182, 4.2139, 4.2457, 4.2776, 4.3095, 4.3414, 4.3732, 4.4051, 
                                4.437, 4.4506, 4.4825, 4.5143, 4.5462, 4.5781, 4.61, 4.6418, 4.6737, 4.7056, 4.7375, 4.7693, 4.8012, 4.8331, 4.865, 
                                4.8968, 4.9287, 4.9606, 4.9742, 5.0061, 5.0379, 5.0698, 5.1017, 5.1336, 5.1654, 5.1973, 5.2292, 5.2611, 5.2929, 5.3248, 
                                5.3567, 5.3886, 5.4204, 5.4523, 5.4842, 5.4978, 5.5297, 5.5615, 5.5934, 5.6253, 5.6572, 5.689, 5.7209, 5.7528, 5.7847, 
                                5.8165, 5.8484, 5.8803, 5.9122, 5.944, 5.9759, 6.0078, 6.0214, 6.0533, 6.0851, 6.117, 6.1489, 6.1808, 6.2127, 6.2445, 
                                6.2764, 2 * M_PI};


  double z_bins[2 * nz + 1];
  for (int z = 0; z <= 2 * nz; z++)
  {
    z_bins[z] = -z_rdo + z_rdo / nz * z;
  }

  _h_R = new TH1F("_h_R", "_h_R;R, [m]", r_bins_N+2, r_bins_edges);
  _h_hits = new TH1F("_h_hits", "_h_hits;N, [hit]", 1000, 0, 1e6);
  _h_DC_E = new TH2F("_h_DC_E", "_h_DC_E;SC;E#times10^{6}", 2000, -100, 2e5 - 100, 1000, -50, 1e3 - 50);
  char name[100];
  char name_ax[100];
  for (int iz = 0; iz < 30; iz++)
  {
    sprintf(name, "_h_SC_ibf_%d", iz);
    sprintf(name_ax, "_h_SC_ibf_%d;#phi, [rad];R, [m];Z, [m]", iz);
    _h_SC_ibf[iz] = new TH3F(name, name_ax, nphi, phi_bins, r_bins_N, r_bins, 2 * nz, z_bins);
    sprintf(name, "_h_SC_prim_%d", iz);
    sprintf(name_ax, "_h_SC_prim_%d;#phi, [rad];R, [m];Z, [m]", iz);
    _h_SC_prim[iz] = new TH3F(name, name_ax, nphi, phi_bins, r_bins_N, r_bins, 2 * nz, z_bins);

    hm->registerHisto(_h_SC_prim[iz]);
    hm->registerHisto(_h_SC_ibf[iz]);
  }
  hm->registerHisto(_h_hits);
  hm->registerHisto(_h_R);
  hm->registerHisto(_h_DC_E);

  //_event_timestamp = 0;
  _hit_eion = 0;
  _hit_r = 0;
  _hit_phi = 0;
  _hit_z = 0;
  _ibf_vol = 0;
  _amp_ele_vol = 0;
  if (_fSliming == 1)
  {
    _rawHits = new TTree("hTree", "tpc hit tree for ionization");
    _rawHits->Branch("isOnPlane", &_isOnPlane);
    _rawHits->Branch("hit_z", &_hit_z);
    _rawHits->Branch("hit_r", &_hit_r);
    _rawHits->Branch("hit_phi", &_hit_phi);
    _rawHits->Branch("hit_eion", &_hit_eion);
    _rawHits->Branch("ibf_vol", &_ibf_vol);
    _rawHits->Branch("amp_ele_vol", &_amp_ele_vol);

    _rawHits->Branch("event_timestamp", &_event_timestamp);
    _rawHits->Branch("event_bunchXing", &_event_bunchXing);
  }
  return 0;
}

//____________________________________________________________________________..
int fillSpaceChargeMaps::InitRun(PHCompositeNode * /*topNode*/)
{
  //std::cout << "fillSpaceChargeMaps::InitRun(PHCompositeNode *topNode) Initializing for Run XXX" << std::endl;
  std::string line;
  //AA collisions timestamps
  std::string txt_file = "./data/timestamps_50kHz_1M.txt";
  int start_line = 3;
  if (_collSyst == 1)
  {
    //pp collisions timestamps
    txt_file = "/phenix/u/hpereira/sphenix/work/g4simulations/timestamps_3MHz.txt";
    start_line = 2;
  }
  std::ifstream InputFile(txt_file);
  if (InputFile.is_open())
  {
    int n_line = 0;
    while (getline(InputFile, line))
    {
      n_line++;
      if (n_line > start_line)
      {
        std::istringstream is(line);
        double n[2] = {0, 0};
        int i = 0;
        while (is >> n[i])
        {
          i++;
        }
        _timestamps[n[0]] = n[1];
        if (n_line < 10)
        {
          std::cout << n[1] << std::endl;
        }
        _keys.push_back(int(n[0]));
      }
    }
    InputFile.close();
  }

  else
    std::cout << "Unable to open file:" << txt_file << std::endl;

  if (_fUseIBFMap)
  {
    TFile *MapsFile = new TFile("./data/IBF_Map.root", "READ");
    if (MapsFile->IsOpen()) printf("File opened successfully\n");
    _h_modules_anode = (TH2F *) MapsFile->Get("h_modules_anode")->Clone("_h_modules_anode");
    _h_modules_measuredibf = (TH2F *) MapsFile->Get("h_modules_measuredibf")->Clone("_h_modules_measuredibf");
  }

  _mbRate = _freqKhz * kHz;
  _xingRate = 9.383 * MHz;
  _mean = mbRate / xingRate;

  return 0;
}

//____________________________________________________________________________..
int fillSpaceChargeMaps::process_event(PHCompositeNode *topNode)
{
  //double tpc_frame_r3_outer = 759.11*mm;  //758.4;//mm inner edge of larger-r frame of r3
  //double tpc_frame_r3_inner = 583.67*mm;  //583.5;//mm outer edge of smaller-r frame of r3

  //double tpc_frame_r2_outer = 574.76*mm;  //574.9;//mm inner edge of larger-r frame of r3
  //double tpc_frame_r2_inner = 411.53*mm;  //411.4;//mm outer edge of smaller-r frame of r3

  //double tpc_frame_r1_outer = 402.49*mm;  //402.6;//mm inner edge of larger-r frame of r3
  //double tpc_frame_r1_inner = 217.83*mm;  //221.0;//mm outer edge of smaller-r frame of r3


  double z_bias_avg = 0;
  int bemxingsInFile = _keys.size();
  int key = -1;
  if (_fAvg == 1)
  {
    z_bias_avg = 1.055 * m * (float) rand() / RAND_MAX;
  }
  if (_evtstart >= bemxingsInFile) _evtstart = _evtstart - bemxingsInFile;
  key = _keys.at(_evtstart);
  _event_timestamp = (float) _timestamps[key] * ns;  //units in seconds
  _event_bunchXing = key;

  if (_evtstart % 10 == 0) std::cout << "_evtstart = " << _evtstart << std::endl;
  _evtstart++;

  Shifter shifter("/sphenix/user/rcorliss/distortion_maps/2021.04/apr07.average.real_B1.4_E-400.0.ross_phi1_sphenix_phislice_lookup_r26xp40xz40.distortion_map.hist.root");

  PHG4HitContainer *hits = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_TPC");
  int n_hits = 0;
  if (hits)
  {
    PHG4HitContainer::ConstRange hit_range = hits->getHits();
    for (PHG4HitContainer::ConstIterator hit_iter = hit_range.first; hit_iter != hit_range.second; hit_iter++)
    {
      int f_fill_prim[30] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
      int f_fill_ibf[30] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

      float hit_x0 = hit_iter->second->get_x(0);
      float hit_y0 = hit_iter->second->get_y(0);
      float hit_z0 = hit_iter->second->get_z(0);
      float hit_x1 = hit_iter->second->get_x(1);
      float hit_y1 = hit_iter->second->get_y(1);
      float hit_z1 = hit_iter->second->get_z(1);

      float hit_eion = hit_iter->second->get_eion();
      float N_electrons = hit_eion * Tpc_ElectronsPerGeV;
      float x = (hit_x0 + f * (hit_x1 - hit_x0)) * cm;
      float y = (hit_y0 + f * (hit_y1 - hit_y0)) * cm;
      float z = (hit_z0 + f * (hit_z1 - hit_z0)) * cm;

      float r = sqrt(x * x + y * y);
      float phi = atan2(x, y);
      if (phi < 0) phi += 2 * M_PI;

      // Shift electrons according to the field maps
      TVector3 oldPos(x / cm, y / cm, z / cm);
      TVector3 newPos;
      if (_shiftElectrons == 1)
      {
        if (oldPos.z() < 0)
        {
          oldPos.SetZ(abs(oldPos.z()));
          newPos = shifter.ShiftForward(oldPos);
          newPos.SetZ(newPos.z() * -1);
        }
        else
        {
          newPos = shifter.ShiftForward(oldPos);
        }
      }
      else
      {
        newPos.SetZ(oldPos.Z());
        newPos.SetX(oldPos.X());
        newPos.SetY(oldPos.Y());
      }

      //Reading IBF and Gain weights according to X-Y position
      float w_ibf = 1.;
      float w_gain = 1.;

      if (_fUseIBFMap)
      {
        int bin_x = _h_modules_anode->GetXaxis()->FindBin(x / mm);
        int bin_y = _h_modules_anode->GetYaxis()->FindBin(y / mm);
        w_ibf = _h_modules_measuredibf->GetBinContent(bin_x, bin_y);
        w_gain = _h_modules_anode->GetBinContent(bin_x, bin_y);
      }
      float ionsPerEle = w_gain * _ampGain * w_ibf * _ampIBFfrac;

      //Check that it is on the frame
      _isOnPlane = 0;
      double dr_bin = -1;
      double dphi_bin = -1;
      double new_phi = newPos.Phi();
      double new_r = newPos.Perp() * cm;
      if (new_phi < 0) new_phi += 2 * M_PI;
      if (!IsOverFrame(new_r, new_phi))
      {
        _isOnPlane = 1;
      }
      else
      {
        std::vector<double> r_phi_bin = putOnPlane(new_r / mm, new_phi);
        dr_bin = r_phi_bin[0];
        dphi_bin = r_phi_bin[1];
      }
      _hit_z = z;
      _hit_r = r;
      _hit_phi = phi;
      _hit_eion = hit_eion;
      _ibf_vol = N_electrons * ionsPerEle;
      _amp_ele_vol = w_gain * _ampGain;
      if (new_r < 210 || new_r > 770) continue;
      if (_fSliming == 1) _rawHits->Fill();
      double z_prim[30] = {-1 * 1e10, -1 * 1e10, -1 * 1e10, -1 * 1e10, -1 * 1e10, -1 * 1e10, -1 * 1e10, -1 * 1e10, -1 * 1e10, -1 * 1e10,
                           -1 * 1e10, -1 * 1e10, -1 * 1e10, -1 * 1e10, -1 * 1e10, -1 * 1e10, -1 * 1e10, -1 * 1e10, -1 * 1e10, -1 * 1e10,
                           -1 * 1e10, -1 * 1e10, -1 * 1e10, -1 * 1e10, -1 * 1e10, -1 * 1e10, -1 * 1e10, -1 * 1e10, -1 * 1e10, -1 * 1e10};
      double z_ibf[30] = {-1 * 1e10, -1 * 1e10, -1 * 1e10, -1 * 1e10, -1 * 1e10, -1 * 1e10, -1 * 1e10, -1 * 1e10, -1 * 1e10, -1 * 1e10,
                          -1 * 1e10, -1 * 1e10, -1 * 1e10, -1 * 1e10, -1 * 1e10, -1 * 1e10, -1 * 1e10, -1 * 1e10, -1 * 1e10, -1 * 1e10,
                          -1 * 1e10, -1 * 1e10, -1 * 1e10, -1 * 1e10, -1 * 1e10, -1 * 1e10, -1 * 1e10, -1 * 1e10, -1 * 1e10, -1 * 1e10};

      if (_hit_z >= 5 * mm && _hit_z < 1.055 * m)
      {
        n_hits++;
        _h_DC_E->Fill(_ibf_vol, hit_eion * 1e6);
        for (int iz = 0; iz < 30; iz++)
        {
          double bX = _beamxing[iz];
          //double bX_end = _beamxing_end[iz];
          if (_event_bunchXing <= bX)
          {
            if (_fAvg == 1)
            {
              z_prim[iz] = _hit_z - z_bias_avg;
              z_ibf[iz] = 1.055 * m - z_bias_avg;
            }
            else
            {
              z_prim[iz] = _hit_z - (bX - _event_bunchXing) * 106 * vIon * ns;
              z_ibf[iz] = 1.055 * m - (bX - _event_bunchXing) * 106 * vIon * ns;
            }
            if (z_prim[iz] > 0 && z_prim[iz] < 1.055 * m)
            {
              f_fill_prim[iz] = 1;
            }
            if (z_ibf[iz] > 0 && z_ibf[iz] < 1.055 * m)
            {
              f_fill_ibf[iz] = 1;
            }
          }
        }
      }

      if (_hit_z < -5 * mm && _hit_z > -1.055 * m)
      {
        n_hits++;
        _h_DC_E->Fill(_ibf_vol, hit_eion * 1e6);
        for (int iz = 0; iz < 30; iz++)
        {
          double bX = _beamxing[iz];
          //double bX_end = _beamxing_end[iz];
          if (_event_bunchXing <= bX)
          {
            if (_fAvg == 1)
            {
              z_prim[iz] = _hit_z + z_bias_avg;
              z_ibf[iz] = -1.055 * m + z_bias_avg;
            }
            else
            {
              z_prim[iz] = _hit_z + (bX - _event_bunchXing) * 106 * vIon * ns;
              z_ibf[iz] = -1.055 * m + (bX - _event_bunchXing) * 106 * vIon * ns;
            }
            if (z_prim[iz] < 0 && z_prim[iz] > -1.055 * m)
            {
              f_fill_prim[iz] = 1;
            }
            if (z_ibf[iz] < 0 && z_ibf[iz] > -1.055 * m)
            {
              f_fill_ibf[iz] = 1;
            }
          }
        }
      }

      double w_prim = _hit_eion * Tpc_ElectronsPerGeV;
      for (int iz = 0; iz < 30; iz++)
      {
        if (f_fill_prim[iz] == 1)
        {
          _h_SC_prim[iz]->Fill(_hit_phi, _hit_r, z_prim[iz], w_prim);
        }
        if (f_fill_ibf[iz] == 1)
        {
          if (!_isOnPlane)
          {
            //Redistribute charges
            //std::vector<double> newWeights = getNewWeights(_h_SC_ibf[iz], _h_modules_anode, _h_modules_measuredibf, _hit_r, _hit_phi, dr_bin, dphi_bin, _fUseIBFMap);
            std::vector<double> newWeights = getNewWeights(_h_SC_ibf[iz], _h_modules_anode, _h_modules_measuredibf, new_r, new_phi, dr_bin, dphi_bin, _fUseIBFMap);

            double w_ibf_tmp = newWeights[0];
            double w_gain_tmp = newWeights[1];
            _hit_r = newWeights[2];
            _hit_phi = newWeights[3];

            _ibf_vol = N_electrons * w_gain_tmp * _ampGain * w_ibf_tmp * _ampIBFfrac;
            _h_SC_ibf[iz]->Fill(_hit_phi, _hit_r, z_ibf[iz], _ibf_vol);
            
          }
          else
          {
            _h_SC_ibf[iz]->Fill(new_phi, new_r, z_ibf[iz], _ibf_vol);
          }
        }
      }

      if (f_fill_ibf[0] == 1)
      {
        _h_R->Fill(_hit_r);
      }
    }
  }
  else
  {
    if (_fSliming == 1) _rawHits->Fill();
  }
  _h_hits->Fill(n_hits);

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int fillSpaceChargeMaps::End(PHCompositeNode * /*topNode*/)
{
  std::cout << "fillSpaceChargeMaps::End" << std::endl;
  if (_fSliming == 1)
  {
    outfile->cd();
    outfile->Write();
    outfile->Close();
    delete outfile;
    for (int iz = 0; iz < 30; iz++)
    {
      _h_SC_prim[iz]->Sumw2(false);
      _h_SC_ibf[iz]->Sumw2(false);
    }

    _h_hits->Sumw2(false);
    _h_DC_E->Sumw2(false);
    _h_R->Sumw2(false);
    hm->dumpHistos(_filename, "UPDATE");
  }
  else
  {
    std::cout << "Writing File:" << _filename << std::endl;
    hm->dumpHistos(_filename, "RECREATE");
  }

  return 0;
}

void fillSpaceChargeMaps::SetFrequency(int freq)
{
  _freqKhz = freq;
  std::cout << "Frequency is set to: " << _freqKhz << " kHz" << std::endl;
}

void fillSpaceChargeMaps::SetBeamXing(const std::vector<int> &beamXs)
{
  _beamxing = beamXs;
  std::cout << "Initial BeamXing is set to: " << _beamxing[0] << std::endl;
}
//void fillSpaceChargeMaps::SetBeamXingEnd(std::vector<int> beamXs_end){
//  _beamxing_end = beamXs_end;
//  std::cout<<"Last BeamXing to fill TPC is set to: "<<_beamxing_end[0]<<std::endl;
//
//}
void fillSpaceChargeMaps::SetEvtStart(int newEvtStart)
{
  _evtstart = newEvtStart;
  std::cout << "Starting event is set to: " << newEvtStart << std::endl;
}

void fillSpaceChargeMaps::SetUseIBFMap(bool useIBFMap)
{
  _fUseIBFMap = useIBFMap;
  std::cout << "Using IBF and Gain Maps" << std::endl;
}
void fillSpaceChargeMaps::SetGain(float ampGain)
{
  _ampGain = ampGain;
  std::cout << "Gain is set to: " << _ampGain << std::endl;
}
void fillSpaceChargeMaps::SetIBF(float ampIBFfrac)
{
  _ampIBFfrac = ampIBFfrac;
  std::cout << "IBF is set to: " << _ampIBFfrac << std::endl;
}

void fillSpaceChargeMaps::SetCollSyst(int coll_syst)
{
  _collSyst = coll_syst;
  static const std::string s_syst[2] = {"AA", "pp"};
  std::cout << "Collision system is set to: " << s_syst[_collSyst] << std::endl;
}

void fillSpaceChargeMaps::SetAvg(int fAvg)
{
  _fAvg = fAvg;
  static const std::string s_avg[2] = {"OFF", "ON"};
  std::cout << "Averaging is set to: " << s_avg[_fAvg] << std::endl;
}
void fillSpaceChargeMaps::UseSliming(int fSliming)
{
  _fSliming = fSliming;
  static const std::string s_sliming[2] = {"OFF", "ON"};
  std::cout << "Sliming is set to: " << s_sliming[_fSliming] << std::endl;
}

void fillSpaceChargeMaps::UseFieldMaps(int shiftElectrons)
{
  _shiftElectrons = shiftElectrons;
  static const std::string s_shiftElectrons[2] = {"OFF", "ON"};
  std::cout << "Use Field-Maps is set to: " << s_shiftElectrons[_shiftElectrons] << std::endl;
}

std::vector<double> fillSpaceChargeMaps::getNewWeights(TH3 *h_SC_ibf, TH2 *h_modules_anode, TH2 *h_modules_measuredibf, double hit_r, double hit_phi, double dr_bin, double dphi_bin, bool fUseIBFMap)
{
  double w_ibf_tmp = 1.0;
  double w_gain_tmp = 1.0;
  int r_bin = _h_R->GetXaxis()->FindBin(dr_bin);
  double r_bin_c = _h_R->GetXaxis()->GetBinCenter(r_bin);
  double r_bin_r = _h_R->GetXaxis()->GetBinCenter(r_bin+1);
  double r_bin_l = _h_R->GetXaxis()->GetBinCenter(r_bin-1);

  int phi_bin = h_SC_ibf->GetXaxis()->FindBin(dphi_bin);
    
  double phi_bin_c = h_SC_ibf->GetXaxis()->GetBinCenter(phi_bin);
  double phi_bin_r = h_SC_ibf->GetXaxis()->GetBinCenter(phi_bin+1);
  double phi_bin_l = h_SC_ibf->GetXaxis()->GetBinCenter(phi_bin-1);



  if (dr_bin > 0 && dphi_bin < 0)
  {
    if (hit_r - r_bin_c > 0)
    {
      hit_r = r_bin_r;//hit_r + dr_bin;
    }
    else
    {
      hit_r = r_bin_l;//hit_r - dr_bin;
    }
  }

  if (dr_bin < 0 && dphi_bin > 0)
  {
    if (hit_phi - phi_bin_c > 0)
    {
      hit_phi = phi_bin_r; //hit_phi + dphi_bin;
    }
    else
    {
      hit_phi = phi_bin_l; //hit_phi - dphi_bin;
    }
  }
  if (dr_bin > 0 && dphi_bin > 0)
  {
    if (hit_phi - phi_bin_c >= 0 && hit_r - r_bin_c >= 0)
    {
      hit_phi = phi_bin_r; //hit_phi + dphi_bin;
      hit_r = r_bin_r;//hit_r + dr_bin;
    }
    if (hit_phi - phi_bin_c <= 0 && hit_r - r_bin_c <= 0)
    {
      hit_phi = phi_bin_l; //hit_phi - dphi_bin;
      hit_r = r_bin_l;//hit_r - dr_bin;
    }
    if (hit_phi - phi_bin_c >= 0 && hit_r - r_bin_c <= 0)
    {
      hit_phi = phi_bin_r; //hit_phi + dphi_bin;
      hit_r = r_bin_l;//hit_r - dr_bin;
    }
    if (hit_phi - phi_bin_c <= 0 && hit_r - r_bin_c >= 0)
    {
      hit_phi = phi_bin_l; //hit_phi - dphi_bin;
      hit_r = r_bin_r;//hit_r + dr_bin;
    }
  }
  if (fUseIBFMap)
  {
    double x_tmp = hit_r * cos(hit_phi);
    double y_tmp = hit_r * sin(hit_phi);
    int bin_x = h_modules_anode->GetXaxis()->FindBin(x_tmp);
    int bin_y = h_modules_anode->GetYaxis()->FindBin(y_tmp);
    w_ibf_tmp = h_modules_measuredibf->GetBinContent(bin_x, bin_y);
    w_gain_tmp = h_modules_anode->GetBinContent(bin_x, bin_y);
  }
  //for(int p=0;p<24;p++){
  //  if(hit_phi>=phi_dead_bins[p] && hit_phi<=phi_dead_bins[p+1]){
  //    std::cout<< 
  //    " dphi_bin  = "<< dphi_bin << 
  //    " phi_bin_c = " <<phi_bin_c <<
  //    " phi_bin_r = " <<phi_bin_r <<
  //    " phi_bin_l = " <<phi_bin_l
  //    << std::endl;
  //    std::cout<<"end new Weights hit_phi = "<<hit_phi<<std::endl;
  //  }
  //  p++;
  //}
  std::vector<double> newWeights;
  newWeights.push_back(w_ibf_tmp);
  newWeights.push_back(w_gain_tmp);
  newWeights.push_back(hit_r);
  newWeights.push_back(hit_phi);
  return newWeights;
}

bool fillSpaceChargeMaps::IsOverFrame(double r, double phi)
{
  //these parameters are taken from Feb 12 drawings of frames.
  //double tpc_frame_side_gap = 0.8 * mm;    //mm //space between radial line and start of frame
  //double tpc_frame_side_width = 2.6 * mm;  //mm //thickness of frame
  double tpc_margin = 2.0 * mm;            //mm // extra gap between edge of frame and start of GEM holes

  double tpc_frame_r3_outer = 759.11 * mm;  //758.4;//mm inner edge of larger-r frame of r3
  double tpc_frame_r3_inner = 583.67 * mm;  //583.5;//mm outer edge of smaller-r frame of r3

  double tpc_frame_r2_outer = 574.76 * mm;  //574.9;//mm inner edge of larger-r frame of r3
  double tpc_frame_r2_inner = 411.53 * mm;  //411.4;//mm outer edge of smaller-r frame of r3

  double tpc_frame_r1_outer = 402.49 * mm;  //402.6;//mm inner edge of larger-r frame of r3
  double tpc_frame_r1_inner = 217.83 * mm;  //221.0;//mm outer edge of smaller-r frame of r3

  //double tpc_sec0_phi=0.0;//get_double_param("tpc_sec0_phi");

  //if the coordinate is in the radial spaces of the frames, return true:
  if (r <= tpc_frame_r1_inner + tpc_margin)
  {
    return true;
  }
  if (r >= tpc_frame_r1_outer - tpc_margin && r <= tpc_frame_r2_inner + tpc_margin)
  {
    return true;
  }
  if (r >= tpc_frame_r2_outer - tpc_margin && r <= tpc_frame_r3_inner + tpc_margin)
  {
    return true;
  }
  if (r >= tpc_frame_r3_outer - tpc_margin)
  {
    return true;
  }

  for(int p=0;p<24;p=p+2){
    if(phi>=phi_dead_bins[p] && phi<=phi_dead_bins[p+1]){
      return true;
    }
  }
  return false;
}

std::vector<double> fillSpaceChargeMaps::putOnPlane(double r, double phi)
{
  //these parameters are taken from Feb 12 drawings of frames.
  //double tpc_frame_side_gap = 0.8*mm;    //mm //space between radial line and start of frame
  //double tpc_frame_side_width = 2.6*mm;  //mm //thickness of frame
  double tpc_margin = 2.0*mm;            //mm // extra gap between edge of frame and start of GEM holes

  double tpc_frame_r3_outer = 759.11*mm;  //758.4;//mm inner edge of larger-r frame of r3
  double tpc_frame_r3_inner = 583.67*mm;  //583.5;//mm outer edge of smaller-r frame of r3

  double tpc_frame_r2_outer = 574.76*mm;  //574.9;//mm inner edge of larger-r frame of r3
  double tpc_frame_r2_inner = 411.53*mm;  //411.4;//mm outer edge of smaller-r frame of r3

  double tpc_frame_r1_outer = 402.49*mm;  //402.6;//mm inner edge of larger-r frame of r3
  double tpc_frame_r1_inner = 217.83*mm;  //221.0;//mm outer edge of smaller-r frame of r3

  //double tpc_sec0_phi=0.0;//get_double_param("tpc_sec0_phi");
  double r_0_bin = -1;
  double phi_0_bin = -1;
  //if the coordinate is in the radial spaces of the frames, return true:
  if (r >= tpc_frame_r1_inner - 2*tpc_margin && r <= tpc_frame_r1_inner + tpc_margin)
  {
    //Collect charges from inner margin + 1 mm
    r_0_bin = tpc_frame_r1_inner - tpc_margin/2;//- r; 
  }
  if (r >= tpc_frame_r1_outer - tpc_margin && r <= tpc_frame_r2_inner + tpc_margin)
  {
    r_0_bin = tpc_frame_r2_inner-(tpc_frame_r2_inner - tpc_frame_r1_outer)/2;
  }
  if (r >= tpc_frame_r2_outer - tpc_margin && r <= tpc_frame_r3_inner + tpc_margin)
  {
    r_0_bin = tpc_frame_r3_inner-(tpc_frame_r3_inner - tpc_frame_r2_outer)/2;
  }
  if (r >= tpc_frame_r3_outer - tpc_margin)
  {
    //Collect charges from outer margin + 1 mm
    r_0_bin = tpc_frame_r3_outer + tpc_margin/2;
  }

  //if the coordinate is within gap+width of a sector boundary, return true:
  //note that this is not a line of constant radius, but a linear distance from a radius.

  //find the two spokes we're between:
  for(int p=0;p<24;p=p+2){
    if(phi>=phi_dead_bins[p] && phi<=phi_dead_bins[p+1]){
      phi_0_bin = phi_dead_bins[p] + (phi_dead_bins[p+1]-phi_dead_bins[p])/2;
    }
  }

  // // float sectorangle = (M_PI / 6.);
  // // float nsectors = (phi+sectorangle/2) / sectorangle;
  // // if(nsectors>11) nsectors-=11;
  // // int nsec = floor(nsectors);
  // // float reduced_phi = phi - nsec * sectorangle;  //between zero and sixty degrees.
  // // float dist_to_previous = r * sin(reduced_phi );
  // // float dist_to_next = r * sin(sectorangle - reduced_phi );
  // // if (dist_to_previous <= tpc_frame_side_gap + tpc_frame_side_width + tpc_margin)
  // // {
  // //   //phi_0_bin = 0.0136;
  // //   phi_0_bin = asin((tpc_frame_side_gap + tpc_frame_side_width + tpc_margin) / r);
  // // }
  // // if (dist_to_next <= tpc_frame_side_gap + tpc_frame_side_width + tpc_margin)
  // // {
  // //   //phi_0_bin = 0.0136;
  // //   phi_0_bin = asin((tpc_frame_side_gap + tpc_frame_side_width + tpc_margin) / r);
  // // }
  std::vector<double> r_phi_bin;
  r_phi_bin.push_back(r_0_bin);
  r_phi_bin.push_back(phi_0_bin);
  return r_phi_bin;
}
