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
#include <cmath>      // for sin, asin, cos, floor, M_PI
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
  const int nphi = 205;
  double phi_bins[nphi + 1] = {0., 0.0068, 0.038675, 0.07055, 0.102425, 0.1343, 0.166175, 0.19805,
                               0.229925, 0.2618, 0.293675, 0.32555, 0.357425, 0.3893, 0.421175, 0.45305, 0.484925,
                               0.5168, 0.5304, 0.562275, 0.59415, 0.626025, 0.6579, 0.689775, 0.72165, 0.753525, 0.7854,
                               0.817275, 0.84915, 0.881025, 0.9129, 0.944775, 0.97665, 1.008525, 1.0404, 1.054, 1.085875,
                               1.11775, 1.149625, 1.1815, 1.213375, 1.24525, 1.277125, 1.309, 1.340875, 1.37275, 1.404625, 1.4365,
                               1.468375, 1.50025, 1.532125, 1.564, 1.5776, 1.609475, 1.64135, 1.673225, 1.7051, 1.736975, 1.76885,
                               1.800725, 1.8326, 1.864475, 1.89635, 1.928225, 1.9601, 1.991975, 2.02385, 2.055725, 2.0876, 2.1012,
                               2.133075, 2.16495, 2.196825, 2.2287, 2.260575, 2.29245, 2.324325, 2.3562, 2.388075, 2.41995, 2.451825,
                               2.4837, 2.515575, 2.54745, 2.579325, 2.6112, 2.6248, 2.656675, 2.68855, 2.720425, 2.7523, 2.784175, 2.81605,
                               2.847925, 2.8798, 2.911675, 2.94355, 2.975425, 3.0073, 3.039175, 3.07105, 3.102925, 3.1348, 3.1484, 3.180275,
                               3.21215, 3.244025, 3.2759, 3.307775, 3.33965, 3.371525, 3.4034, 3.435275, 3.46715, 3.499025, 3.5309, 3.562775,
                               3.59465, 3.626525, 3.6584, 3.672, 3.703875, 3.73575, 3.767625, 3.7995, 3.831375, 3.86325, 3.895125, 3.927, 3.958875,
                               3.99075, 4.022625, 4.0545, 4.086375, 4.11825, 4.150125, 4.182, 4.1956, 4.227475, 4.25935, 4.291225, 4.3231, 4.354975,
                               4.38685, 4.418725, 4.4506, 4.482475, 4.51435, 4.546225, 4.5781, 4.609975, 4.64185, 4.673725, 4.7056, 4.7192, 4.751075,
                               4.78295, 4.814825, 4.8467, 4.878575, 4.91045, 4.942325, 4.9742, 5.006075, 5.03795, 5.069825, 5.1017, 5.133575, 5.16545,
                               5.197325, 5.2292, 5.2428, 5.274675, 5.30655, 5.338425, 5.3703, 5.402175, 5.43405, 5.465925, 5.4978, 5.529675, 5.56155,
                               5.593425, 5.6253, 5.657175, 5.68905, 5.720925, 5.7528, 5.7664, 5.798275, 5.83015, 5.862025, 5.8939, 5.925775, 5.95765,
                               5.989525, 6.0214, 6.053275, 6.08515, 6.117025, 6.1489, 6.180775, 6.21265, 6.244525, 6.2764, 2 * M_PI};

  double z_bins[2 * nz + 1];
  for (int z = 0; z <= 2 * nz; z++)
  {
    z_bins[z] = -z_rdo + z_rdo / nz * z;
  }

  _h_R = new TH1F("_h_R", "_h_R;R, [m]", r_bins_N, r_bins);
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
        //std::cout<<"dr_bin="<<dr_bin<<"dphi_bin="<<dphi_bin<<std::endl;
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
            std::vector<double> newWeights = getNewWeights(_h_SC_ibf[iz], _h_modules_anode, _h_modules_measuredibf, _hit_r, _hit_phi, dr_bin, dphi_bin, _fUseIBFMap);
            //std::vector<double> newWeights = getNewWeights(_h_SC_ibf[iz], _h_modules_anode, _h_modules_measuredibf, new_r, new_phi, dr_bin, dphi_bin, _fUseIBFMap);

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
  int r_bin = h_SC_ibf->GetYaxis()->FindBin(hit_r);
  double r_bin_c = h_SC_ibf->GetYaxis()->GetBinCenter(r_bin);

  int phi_bin = h_SC_ibf->GetXaxis()->FindBin(hit_phi);
  double phi_bin_c = h_SC_ibf->GetXaxis()->GetBinCenter(phi_bin);

  if (dr_bin > 0 && dphi_bin < 0)
  {
    if (hit_r - r_bin_c > 0)
    {
      hit_r = hit_r + dr_bin;
    }
    else
    {
      hit_r = hit_r - dr_bin;
    }
  }

  if (dr_bin < 0 && dphi_bin > 0)
  {
    if (hit_phi - phi_bin_c > 0)
    {
      hit_phi = hit_phi + dphi_bin;
    }
    else
    {
      hit_phi = hit_phi - dphi_bin;
    }
  }
  if (dr_bin > 0 && dphi_bin > 0)
  {
    if (hit_phi - phi_bin_c > 0 && hit_r - r_bin_c >= 0)
    {
      hit_phi = hit_phi + dphi_bin;
      hit_r = hit_r + dr_bin;
    }
    if (hit_phi - phi_bin_c <= 0 && hit_r - r_bin_c < 0)
    {
      hit_phi = hit_phi - dphi_bin;
      hit_r = hit_r - dr_bin;
    }
    if (hit_phi - phi_bin_c > 0 && hit_r - r_bin_c <= 0)
    {
      hit_phi = hit_phi + dphi_bin;
      hit_r = hit_r - dr_bin;
    }
    if (hit_phi - phi_bin_c <= 0 && hit_r - r_bin_c > 0)
    {
      hit_phi = hit_phi - dphi_bin;
      hit_r = hit_r + dr_bin;
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
  double tpc_frame_side_gap = 0.8;    //mm //space between radial line and start of frame
  double tpc_frame_side_width = 2.6;  //mm //thickness of frame
  double tpc_margin = 2.0;            //mm // extra gap between edge of frame and start of GEM holes

  double tpc_frame_r3_outer = 759.11;  //758.4;//mm inner edge of larger-r frame of r3
  double tpc_frame_r3_inner = 583.67;  //583.5;//mm outer edge of smaller-r frame of r3

  double tpc_frame_r2_outer = 574.76;  //574.9;//mm inner edge of larger-r frame of r3
  double tpc_frame_r2_inner = 411.53;  //411.4;//mm outer edge of smaller-r frame of r3

  double tpc_frame_r1_outer = 402.49;  //402.6;//mm inner edge of larger-r frame of r3
  double tpc_frame_r1_inner = 217.83;  //221.0;//mm outer edge of smaller-r frame of r3

  //double tpc_sec0_phi=0.0;//get_double_param("tpc_sec0_phi");

  //if the coordinate is in the radial spaces of the frames, return true:
  if (r < tpc_frame_r1_inner + tpc_margin)
  {
    return true;
  }
  if (r > tpc_frame_r1_outer - tpc_margin && r < tpc_frame_r2_inner + tpc_margin)
  {
    return true;
  }
  if (r > tpc_frame_r2_outer - tpc_margin && r < tpc_frame_r3_inner + tpc_margin)
  {
    return true;
  }
  if (r > tpc_frame_r3_outer - tpc_margin)
  {
    return true;
  }
  //if the coordinate is within gap+width of a sector boundary, return true:
  //note that this is not a line of constant radius, but a linear distance from a radius.

  //find the two spokes we're between:

  float sectorangle = (M_PI / 6.);
  float nsectors = phi / sectorangle;
  int nsec = floor(nsectors);
  float reduced_phi = phi - nsec * sectorangle;  //between zero and sixty degrees.
  float dist_to_previous = r * sin(reduced_phi);
  float dist_to_next = r * sin(sectorangle - reduced_phi);
  if (dist_to_previous < tpc_frame_side_gap + tpc_frame_side_width + tpc_margin)
  {
    return true;
  }
  if (dist_to_next < tpc_frame_side_gap + tpc_frame_side_width + tpc_margin)
  {
    return true;
  }
  return false;
}

std::vector<double> fillSpaceChargeMaps::putOnPlane(double r, double phi)
{
  //these parameters are taken from Feb 12 drawings of frames.
  double tpc_frame_side_gap = 0.8;    //mm //space between radial line and start of frame
  double tpc_frame_side_width = 2.6;  //mm //thickness of frame
  double tpc_margin = 2.0;            //mm // extra gap between edge of frame and start of GEM holes

  double tpc_frame_r3_outer = 759.11;  //758.4;//mm inner edge of larger-r frame of r3
  double tpc_frame_r3_inner = 583.67;  //583.5;//mm outer edge of smaller-r frame of r3

  double tpc_frame_r2_outer = 574.76;  //574.9;//mm inner edge of larger-r frame of r3
  double tpc_frame_r2_inner = 411.53;  //411.4;//mm outer edge of smaller-r frame of r3

  double tpc_frame_r1_outer = 402.49;  //402.6;//mm inner edge of larger-r frame of r3
  double tpc_frame_r1_inner = 217.83;  //221.0;//mm outer edge of smaller-r frame of r3

  //double tpc_sec0_phi=0.0;//get_double_param("tpc_sec0_phi");
  double r_0_bin = -1;
  double phi_0_bin = -1;
  //if the coordinate is in the radial spaces of the frames, return true:
  if (r < tpc_frame_r1_inner + tpc_margin)
  {
    r_0_bin = tpc_frame_r1_inner - r;
  }
  if (r > tpc_frame_r1_outer - tpc_margin && r < tpc_frame_r2_inner + tpc_margin)
  {
    r_0_bin = tpc_frame_r2_inner - tpc_frame_r1_outer;
  }
  if (r > tpc_frame_r2_outer - tpc_margin && r < tpc_frame_r3_inner + tpc_margin)
  {
    r_0_bin = tpc_frame_r3_inner - tpc_frame_r2_outer;
  }
  if (r > tpc_frame_r3_outer - tpc_margin)
  {
    r_0_bin = r - tpc_frame_r3_outer;
  }

  //if the coordinate is within gap+width of a sector boundary, return true:
  //note that this is not a line of constant radius, but a linear distance from a radius.

  //find the two spokes we're between:

  float sectorangle = (M_PI / 6.);
  float nsectors = phi / sectorangle;
  int nsec = floor(nsectors);
  float reduced_phi = phi - nsec * sectorangle;  //between zero and sixty degrees.
  float dist_to_previous = r * sin(reduced_phi);
  float dist_to_next = r * sin(sectorangle - reduced_phi);
  if (dist_to_previous < tpc_frame_side_gap + tpc_frame_side_width + tpc_margin)
  {
    //phi_0_bin = 0.0136;
    phi_0_bin = asin((tpc_frame_side_gap + tpc_frame_side_width + tpc_margin) / r);
  }
  if (dist_to_next < tpc_frame_side_gap + tpc_frame_side_width + tpc_margin)
  {
    //phi_0_bin = 0.0136;
    phi_0_bin = asin((tpc_frame_side_gap + tpc_frame_side_width + tpc_margin) / r);
  }
  std::vector<double> r_phi_bin;
  r_phi_bin.push_back(r_0_bin / 2);
  r_phi_bin.push_back(phi_0_bin);
  return r_phi_bin;
}
