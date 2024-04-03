#include "readDigitalCurrents.h"

#include <g4detectors/PHG4CylinderCellGeom.h>
#include <g4detectors/PHG4CylinderCellGeomContainer.h>
#include <g4detectors/PHG4TpcCylinderGeom.h>
#include <g4detectors/PHG4TpcCylinderGeomContainer.h>

#include <trackbase/TpcDefs.h>
#include <trackbase/TrkrDefs.h>
#include <trackbase/TrkrHit.h>
#include <trackbase/TrkrHitSet.h>
#include <trackbase/TrkrHitSetContainer.h>

#include <fun4all/Fun4AllHistoManager.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>  // for SubsysReco

#include <phool/getClass.h>
#include <phool/phool.h>  // for PHWHERE

#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>

#include <algorithm>  // for max
#include <cmath>      // for M_PI, sin, cos
#include <fstream>
#include <iostream>
#include <set>
#include <sstream>
#include <string>
#include <utility>  // for pair

bool IsOverFrame(double r, double phi);

bool IsOverFrame(double r, double phi)
{
  // these parameters are taken from Feb 12 drawings of frames.
  double tpc_frame_side_gap = 0.8;    // mm //space between radial line and start of frame
  double tpc_frame_side_width = 2.6;  // mm //thickness of frame
  double tpc_margin = 0.0;            // mm // extra gap between edge of frame and start of GEM holes

  double tpc_frame_r3_outer = 758.4;  // mm inner edge of larger-r frame of r3
  double tpc_frame_r3_inner = 583.5;  // mm outer edge of smaller-r frame of r3

  double tpc_frame_r2_outer = 574.9;  // mm inner edge of larger-r frame of r2
  double tpc_frame_r2_inner = 411.4;  // mm outer edge of smaller-r frame of r2

  double tpc_frame_r1_outer = 402.6;  // mm inner edge of larger-r frame of r1
  double tpc_frame_r1_inner = 221.0;  // mm outer edge of smaller-r frame of r1

  // double tpc_sec0_phi=0.0;//get_double_param("tpc_sec0_phi");

  // if the coordinate is in the radial spaces of the frames, return true:
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

  // if the coordinate is within gap+width of a sector boundary, return true:
  // note that this is not a line of constant radius, but a linear distance from a radius.

  // find the two spokes we're between:

  double sectorangle = (M_PI / 6);
  double nsectors = phi / sectorangle;
  int nsec = std::floor(nsectors);
  double reduced_phi = phi - nsec * sectorangle;  // between zero and sixty degrees.
  double dist_to_previous = r * std::sin(reduced_phi);
  double dist_to_next = r * std::sin(sectorangle - reduced_phi);
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

//____________________________________________________________________________..
readDigitalCurrents::readDigitalCurrents(const std::string &name, const std::string &filename)
  : SubsysReco(name)
  , hm(nullptr)
  , _filename(filename)
  , _ampIBFfrac(0.02)
  , _collSyst(0)
// , outfile(nullptr)
{
  std::cout << "readDigitalCurrents::readDigitalCurrents(const std::string &name) Calling ctor" << std::endl;
}

//____________________________________________________________________________..
readDigitalCurrents::~readDigitalCurrents()
{
  std::cout << "readDigitalCurrents::~readDigitalCurrents() Calling dtor" << std::endl;
  delete hm;
}

//____________________________________________________________________________..
int readDigitalCurrents::Init(PHCompositeNode * /*topNode*/)
{
  std::cout << "readDigitalCurrents::Init(PHCompositeNode *topNode) Initializing" << std::endl;

  int nz = 72;
  double z_rdo = 108 * cm;

  int nr = 159;
  // const int nphi=128*3;
  // const int nz=62*2;
  // double z_rdo=105.5*cm;
  // double rmin=20*cm;
  double rmax = 78 * cm;

  hm = new Fun4AllHistoManager("HITHIST");
  const int r_bins_N = 66;  // 51;
  double r_bins[r_bins_N + 1] = {217.83,
                                 223.83, 229.83, 235.83, 241.83, 247.83, 253.83, 259.83, 265.83, 271.83, 277.83, 283.83, 289.83, 295.83, 301.83, 306.83,
                                 311.05, 317.92, 323.31, 329.27, 334.63, 340.59, 345.95, 351.91, 357.27, 363.23, 368.59, 374.55, 379.91, 385.87, 391.23, 397.19, 402.49,
                                 411.53, 421.70, 431.90, 442.11, 452.32, 462.52, 472.73, 482.94, 493.14, 503.35, 513.56, 523.76, 533.97, 544.18, 554.39, 564.59, 574.76,
                                 583.67, 594.59, 605.57, 616.54, 627.51, 638.48, 649.45, 660.42, 671.39, 682.36, 693.33, 704.30, 715.27, 726.24, 737.21, 748.18, 759.11};
  const int nphi = 205;
  double phi_bins[nphi + 1] = {0., 6.3083 - 2 * M_PI, 6.3401 - 2 * M_PI, 6.372 - 2 * M_PI, 6.4039 - 2 * M_PI, 6.4358 - 2 * M_PI, 6.4676 - 2 * M_PI, 6.4995 - 2 * M_PI, 6.5314 - 2 * M_PI,
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

  // R0 sector 0: 2.36679 < phi < 2.86919
  // R0 sector 1: 1.8432 < phi < 2.3456
  // R0 sector 2: 1.3196 < phi < 1.822
  // R0 sector 3: 0.795998 < phi < 1.2984
  // R0 sector 4: 0.272399 < phi < 0.774799
  // R0 sector 5: -0.2512 < phi < 0.2512
  // R0 sector 6: -0.774799 < phi < -0.272399
  // R0 sector 7: -1.2984 < phi < -0.795998
  // R0 sector 8: -1.822 < phi < -1.3196
  // R0 sector 9: -2.3456 < phi < -1.8432
  // R0 sector 10: -2.86919 < phi < -2.36679
  // R0 sector 11: -3.39279 < phi < -2.89039

  double z_bins[2 * nz + 1];
  for (int z = 0; z <= 2 * nz; z++)
  {
    z_bins[z] = -z_rdo + z_rdo / nz * z;
  }

  _h_R = new TH1F("_h_R", "_h_R;R, [m]", r_bins_N, r_bins);
  _h_hits = new TH1F("_h_hits", "_h_hits;N, [hit]", 1e5, 0 - 0.5, 1e5 - 0.5);
  _h_hit_XY = new TH2F("_h_hit_XY", "_h_hit_XY;X, [m];Y, [m]", 4 * nr, -1 * rmax, rmax, 4 * nr, -1 * rmax, rmax);

  //_h_SC_ibf  = new TH3F("_h_SC_ibf" ,"_h_SC_ibf;#phi, [rad];R, [mm];Z, [mm]"   ,nphi,phi_bins,r_bins_N ,r_bins,2*nz,z_bins);
  _h_DC_SC = new TH3F("_h_DC_SC", "_h_DC_SC;#phi, [rad];R, [mm];Z, [mm]", nphi, phi_bins, r_bins_N, r_bins, 2 * nz, z_bins);
  _h_DC_SC_XY = new TH2F("_h_DC_SC_XY", "_h_DC_SC_XY;X, [mm];Y, [mm];ADC;", 4 * nr, -1 * rmax, rmax, 4 * nr, -1 * rmax, rmax);
  _h_DC_E = new TH2F("_h_DC_E", "_h_DC_E;ADC;E", 200, -100, 2e3 - 100, 500, -100, 5e3 - 100);

  std::string name;
  std::string name_ax;
  for (int iz = 0; iz < nFrames; iz++)
  {
    name = "_h_SC_ibf_" + std::to_string(iz);
    name_ax = "_h_SC_ibf_" + std::to_string(iz) + ";#phi, [rad];R, [mm];Z, [mm]";
    _h_SC_ibf[iz] = new TH3F(name.c_str(), name_ax.c_str(), nphi, phi_bins, r_bins_N, r_bins, 2 * nz, z_bins);

    hm->registerHisto(_h_SC_ibf[iz]);
  }

  hm->registerHisto(_h_R);
  hm->registerHisto(_h_hits);
  hm->registerHisto(_h_DC_SC);
  hm->registerHisto(_h_DC_SC_XY);
  hm->registerHisto(_h_hit_XY);
  hm->registerHisto(_h_DC_E);

  _event_timestamp = 0;

  if (_fillCSVFile)
  {
    myCSVFile.open("./Files/example_1ms_120evts_AA.csv");
    myCSVFile << "Event,"
              << "T,"
              << "Pad,"
              << "Radius,"
              << "ADC"
              << "\n";
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int readDigitalCurrents::InitRun(PHCompositeNode * /*topNode*/)
{
  std::cout << "readDigitalCurrents::InitRun(PHCompositeNode *topNode) Initializing for Run XXX" << std::endl;
  std::string line;
  // AA collisions timestamps
  // std::string txt_file = "/sphenix/user/shulga/Work/IBF/DistortionMap/timestamps_50kHz.txt";
  std::string txt_file = "/sphenix/user/shulga/Work/TpcPadPlane_phi_coresoftware/coresoftware/calibrations/tpc/fillSpaceChargeMaps/data/timestamps_50kHz_1M.txt";
  int start_line = 3;
  if (_collSyst == 1)
  {
    // pp collisions timestamps
    txt_file = "/phenix/u/hpereira/sphenix/work/g4simulations/timestamps_3MHz.txt";
    // txt_file = "/sphenix/user/shulga/Work/IBF/DistortionMap/timestamps_50kHz.txt";
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
  {
    std::cout << "Unable to open file:" << txt_file << std::endl;
  }

  TFile *MapsFile;
  // if(_fUseIBFMap){
  MapsFile = new TFile("/sphenix/user/shulga/Work/IBF/DistortionMap/IBF_Map.root", "READ");
  if (MapsFile->IsOpen())
  {
    std::cout << "Gain/IBF Maps File opened successfully" << std::endl;
  }
  //_h_modules_anode       = (TH2F*)MapsFile ->Get("h_modules_anode")      ->Clone("_h_modules_anode");
  _h_modules_measuredibf = (TH2F *) MapsFile->Get("h_modules_measuredibf")->Clone("_h_modules_measuredibf");
  //}
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int readDigitalCurrents::process_event(PHCompositeNode *topNode)
{
  // double bX = _beamxing;
  // float bX = 1508071;
  // double z_bias_avg = 0;
  // if (_fAvg==1){
  //   z_bias_avg=1.05*(float) rand()/RAND_MAX;
  // }
  int bemxingsInFile = _keys.size();
  if (_evtstart >= bemxingsInFile)
  {
    _evtstart = _evtstart - bemxingsInFile;
  }
  int key = _keys.at(_evtstart);
  _event_timestamp = (float) _timestamps[key] * ns;  // units in seconds
  _event_bunchXing = key;
  if (_evtstart % 100 == 0)
  {
    std::cout << "_evtstart = " << _evtstart << std::endl;
  }
  _evtstart++;

  // std::cout << "readDigitalCurrents::process_event(PHCompositeNode *topNode) Processing Event" << std::endl;
  std::set<std::string>::const_iterator iter;
  // //nodename << "G4HIT_TPC";
  std::string nodename = "TRKR_HITSET";

  // //  SvtxEvaluator
  // SvtxEvaluator *hits = findNode::getClass<SvtxEvaluator>(topNode, nodename.str().c_str());
  // //int n_hits = 0;
  // if (hits){
  //   PHG4HitContainer::ConstRange hit_range = hits->getHits();
  // }
  //===================================
  // get node containing the digitized hits
  TrkrHitSetContainer *_hitmap = findNode::getClass<TrkrHitSetContainer>(topNode, nodename);
  if (!_hitmap)
  {
    std::cout << PHWHERE << "ERROR: Can't find node TRKR_HITSET" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  std::string geo_nodename = "CYLINDERCELLGEOM_SVTX";

  PHG4TpcCylinderGeomContainer *_geom_container_ccgc = nullptr;
  PHG4CylinderCellGeomContainer *_geom_container_cgc = nullptr;
  if (_f_ccgc == 1)
  {
    _geom_container_ccgc = findNode::getClass<PHG4TpcCylinderGeomContainer>(topNode, geo_nodename);
    if (!_geom_container_ccgc)
    {
      std::cout << PHWHERE << "ERROR: Can't find node CYLINDERCELLGEOM_SVTX" << std::endl;
      return Fun4AllReturnCodes::ABORTRUN;
    }
  }
  else
  {
    _geom_container_cgc = findNode::getClass<PHG4CylinderCellGeomContainer>(topNode, geo_nodename);

    if (!_geom_container_cgc)
    {
      std::cout << PHWHERE << "ERROR: Can't find node CYLINDERCELLGEOM_SVTX" << std::endl;
      return Fun4AllReturnCodes::ABORTRUN;
    }
  }

  // loop over all the hits
  // hits are stored in hitsets, so have to get the hitset first
  int n_hits = 0;
  // float _event_bunchXing = 1508071;
  TrkrHitSetContainer::ConstRange all_hitsets = _hitmap->getHitSets();
  for (TrkrHitSetContainer::ConstIterator iter_hitset = all_hitsets.first; iter_hitset != all_hitsets.second; ++iter_hitset)
  {
    // checking that the object is inside TPC
    if (TrkrDefs::getTrkrId(iter_hitset->first) == TrkrDefs::tpcId)
    {
      TrkrDefs::hitsetkey hitsetkey = iter_hitset->first;
      const unsigned int zside = TpcDefs::getSide(hitsetkey);
      TrkrHitSet::ConstRange range = iter_hitset->second->getHits();
      unsigned int layer = TrkrDefs::getLayer(iter_hitset->first);
      // if(layer>6){
      PHG4TpcCylinderGeom *layergeom_ccgc = nullptr;
      PHG4CylinderCellGeom *layergeom_cgc = nullptr;
      double radius = 0;
      if (_f_ccgc == 1)
      {
        layergeom_ccgc = _geom_container_ccgc->GetLayerCellGeom(layer);
        radius = layergeom_ccgc->get_radius() * cm;
      }
      else
      {
        layergeom_cgc = _geom_container_cgc->GetLayerCellGeom(layer);
        radius = layergeom_cgc->get_radius() * cm;
      }
      // PHG4TpcCylinderGeom *layergeom = _geom_container->GetLayerCellGeom(layer);
      // double radius = layergeom->get_radius()*cm;  // returns center of the layer
      int min_phiBin = 1e5;
      int max_phiBin = -1;
      for (TrkrHitSet::ConstIterator hit_iter = range.first; hit_iter != range.second; ++hit_iter)
      {
        // int f_fill_ibf=0;
        int f_fill_ibf[30] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

        unsigned short phibin = TpcDefs::getPad(hit_iter->first);
        unsigned short zbin = TpcDefs::getTBin(hit_iter->first);
        double _drift_velocity = 8.0e-3;  // ActsGeometry::get_drift_velocity();
        double phi_center = 0;
        if (_f_ccgc == 1)
        {
          phi_center = layergeom_ccgc->get_phicenter(phibin);
        }
        else
        {
          phi_center = layergeom_cgc->get_phicenter(phibin);
        }
        if (phi_center < 0)
        {
          phi_center += 2 * M_PI;
        }
        if (phi_center < M_PI / 2 + M_PI / 12 && phi_center > M_PI / 2 - M_PI / 12)
        {
          if (min_phiBin > phibin)
          {
            min_phiBin = phibin;
          }
          if (max_phiBin < phibin)
          {
            max_phiBin = phibin;
          }
        }
        float x = radius * cos(phi_center);
        float y = radius * sin(phi_center);

        double z = 0;  // layergeom->get_zcenter(zbin)*cm;
        if (_f_ccgc == 1)
        {
          z = layergeom_ccgc->get_zcenter(zbin) * _drift_velocity * cm;  //*cm/ns;
          if (zside == 0)
          {
            z = -z;
          }
        }
        else
        {
          z = layergeom_cgc->get_zcenter(zbin) * cm;
        }
        TrkrHit *hit = hit_iter->second;
        int adc = hit->getAdc() - adc_pedestal;
        float E = hit->getEnergy();
        // double z = 0;
        // double z_prim = -1*1e10;
        // double z_ibf =  -1*1e10;
        double z_ibf[30] = {-1 * 1e10, -1 * 1e10, -1 * 1e10, -1 * 1e10, -1 * 1e10, -1 * 1e10, -1 * 1e10, -1 * 1e10, -1 * 1e10, -1 * 1e10,
                            -1 * 1e10, -1 * 1e10, -1 * 1e10, -1 * 1e10, -1 * 1e10, -1 * 1e10, -1 * 1e10, -1 * 1e10, -1 * 1e10, -1 * 1e10,
                            -1 * 1e10, -1 * 1e10, -1 * 1e10, -1 * 1e10, -1 * 1e10, -1 * 1e10, -1 * 1e10, -1 * 1e10, -1 * 1e10, -1 * 1e10};

        int RBin = _h_R->GetXaxis()->FindBin(radius);
        if ((RBin > 33 && RBin < 50) && z > 0)
        {
	  assert(layergeom_ccgc);
          int nRBins = layergeom_ccgc->get_phibins();
          if (phibin < nRBins / 12)
          {
            if (_fillCSVFile)
            {
              myCSVFile << _evtstart << ","
                        << layergeom_ccgc->get_zcenter(zbin) << ","
                        << phibin << ","
                        << RBin - 34 << ","
                        << adc << "\n";
            }
            _h_hit_XY->Fill(x, y);
          }
        }
        // if(!IsOverFrame(radius/mm,phi_center)){
        for (int iz = 0; iz < nFrames; iz++)
        {
          double bX = _beamxing[iz];
          if (_event_bunchXing <= bX)
          {
            if (z >= 0 && z < 1.055 * m)
            {
              z_ibf[iz] = 1.055 * m - (bX - _event_bunchXing) * 106 * vIon * ns;
              if (z_ibf[iz] > 0 && z_ibf[iz] < 1.055 * m)
              {
                f_fill_ibf[iz] = 1;
              }
            }
            if (z < 0 && z > -1.055 * m)
            {
              z_ibf[iz] = -1.055 * m + (bX - _event_bunchXing) * 106 * vIon * ns;
              if (z_ibf[iz] < 0 && z_ibf[iz] > -1.055 * m)
              {
                f_fill_ibf[iz] = 1;
              }
            }
          }
        }
        if (z >= 0 && z < 1.055 * m)
        {
          if (adc >= 0)
          {
            n_hits++;
          }
          if (adc >= 0)
          {
            _h_DC_E->Fill(adc, E);
          }
        }
        if (z < 0 && z > -1.055 * m)
        {
          if (adc >= 0)
          {
            n_hits++;
          }
          if (adc >= 0)
          {
            _h_DC_E->Fill(adc, E);
          }
        }

        // Reading IBF and Gain weights according to X-Y position
        float w_ibf = 1.;
        // float w_gain = 1.;
        // if(_fUseIBFMap){
        int bin_x = _h_modules_measuredibf->GetXaxis()->FindBin(x / mm);
        int bin_y = _h_modules_measuredibf->GetYaxis()->FindBin(y / mm);
        w_ibf = _h_modules_measuredibf->GetBinContent(bin_x, bin_y);
        // w_gain = _h_modules_anode->GetBinContent(bin_x,bin_y);
        w_ibf = 1.;
        //}
        float w_adc = adc * w_ibf;
        _h_DC_SC->Fill(phi_center, radius, z, w_adc);
        _h_DC_SC_XY->Fill(x, y, w_adc);
        if (f_fill_ibf[0] == 1)
        {
          _h_R->Fill(radius);
        }
        for (int iz = 0; iz < nFrames; iz++)
        {
          if (f_fill_ibf[iz] == 1)
          {
            _h_SC_ibf[iz]->Fill(phi_center, radius, z_ibf[iz], w_adc);
          }
        }
        //}
        // if(n_hits%100==0) std::cout<<radius<<"|"<<phi_center<<"|"<<z<<std::endl;
      }
      // std::cout<<" min_phiBin"<< min_phiBin <<" max_phiBin"<< max_phiBin<< std::endl;
    }
  }

  // TrkrHitSetContainer, TrkrHitSet, TrkrHit, and TrkrDefs objects used above in offline/packages/trackbase.
  _h_hits->Fill(n_hits);

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int readDigitalCurrents::ResetEvent(PHCompositeNode * /*topNode*/)
{
  // std::cout << "readDigitalCurrents::ResetEvent(PHCompositeNode *topNode) Resetting internal structures, prepare for next event" << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int readDigitalCurrents::EndRun(const int runnumber)
{
  std::cout << "readDigitalCurrents::EndRun(const int runnumber) Ending Run for Run " << runnumber << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int readDigitalCurrents::End(PHCompositeNode * /*topNode*/)
{
  std::cout << "readDigitalCurrents::End(PHCompositeNode *topNode) This is the End..." << std::endl;
  _h_R->Sumw2(false);
  _h_hits->Sumw2(false);
  _h_DC_E->Sumw2(false);
  _h_DC_SC->Sumw2(false);
  _h_hit_XY->Sumw2(false);
  _h_DC_SC_XY->Sumw2(false);
  for (auto &iz : _h_SC_ibf)
  {
    iz->Sumw2(false);
  }
  hm->dumpHistos(_filename, "RECREATE");
  if (_fillCSVFile)
  {
    myCSVFile.close();
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int readDigitalCurrents::Reset(PHCompositeNode * /*topNode*/)
{
  std::cout << "readDigitalCurrents::Reset(PHCompositeNode *topNode) being Reset" << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
void readDigitalCurrents::Print(const std::string &what) const
{
  std::cout << "readDigitalCurrents::Print(const std::string &what) const Printing info for " << what << std::endl;
}

void readDigitalCurrents::SetEvtStart(int newEvtStart)
{
  _evtstart = newEvtStart;
  std::cout << "Start event is set to: " << newEvtStart << std::endl;
}
void readDigitalCurrents::FillCSV(int fillCSVFile)
{
  _fillCSVFile = fillCSVFile;
  std::cout << "Fill CSV file is set to: " << fillCSVFile << std::endl;
}

void readDigitalCurrents::SetBeamXing(const std::vector<int> &beamXs)
{
  _beamxing = beamXs;
  std::cout << "Initial BeamXing is set to: " << _beamxing[0] << std::endl;
}
void readDigitalCurrents::SetCollSyst(int coll_syst)
{
  _collSyst = coll_syst;
  std::string s_syst[2] = {"AA", "pp"};
  std::cout << "Collision system is set to: " << s_syst[_collSyst] << std::endl;
}

void readDigitalCurrents::SetIBF(double ampIBFfrac)
{
  _ampIBFfrac = ampIBFfrac;
  std::cout << "IBF is set to: " << _ampIBFfrac << std::endl;
}

void readDigitalCurrents::SetCCGC(double f_ccgc)
{
  _f_ccgc = f_ccgc;
  std::cout << "IBF is set to: " << _f_ccgc << std::endl;
}
