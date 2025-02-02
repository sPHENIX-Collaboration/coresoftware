// One-stop header
// Must include first to avoid conflict with "ClassDef" in Rtypes.h
#include <torch/script.h>

#include "TpcClusterizer.h"

#include "LaserEventInfo.h"

#include "TrainingHits.h"
#include "TrainingHitsContainer.h"

#include <trackbase/ClusHitsVerbosev1.h>
#include <trackbase/TpcDefs.h>
#include <trackbase/TrkrClusterContainerv4.h>
#include <trackbase/TrkrClusterHitAssocv3.h>
#include <trackbase/TrkrClusterv3.h>
#include <trackbase/TrkrClusterv4.h>
#include <trackbase/TrkrClusterv5.h>
#include <trackbase/TrkrDefs.h>  // for hitkey, getLayer
#include <trackbase/TrkrHit.h>
#include <trackbase/TrkrHitSet.h>
#include <trackbase/TrkrHitSetContainer.h>
#include <trackbase/alignmentTransformationContainer.h>

#include <trackbase/RawHit.h>
#include <trackbase/RawHitSet.h>
#include <trackbase/RawHitSetContainer.h>
#include <trackbase/RawHitSet.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>  // for SubsysReco

#include <g4detectors/PHG4TpcCylinderGeom.h>
#include <g4detectors/PHG4TpcCylinderGeomContainer.h>

#include <Acts/Definitions/Units.hpp>
#include <Acts/Surfaces/Surface.hpp>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>  // for PHIODataNode
#include <phool/PHNode.h>        // for PHNode
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>  // for PHObject
#include <phool/getClass.h>
#include <phool/phool.h>  // for PHWHERE

#include <TMatrixFfwd.h>    // for TMatrixF
#include <TMatrixT.h>       // for TMatrixT, ope...
#include <TMatrixTUtils.h>  // for TMatrixTRow

#include <TFile.h>

#include <array>
#include <cmath>  // for sqrt, cos, sin
#include <iostream>
#include <limits>
#include <map>  // for _Rb_tree_cons...
#include <string>
#include <utility>  // for pair
#include <vector>
// Terra incognita....
#include <pthread.h>

namespace
{
  template <class T>
  inline constexpr T square(const T &x)
  {
    return x * x;
  }

  using assoc = std::pair<TrkrDefs::cluskey, TrkrDefs::hitkey>;

  struct ihit
  {
    unsigned short iphi = 0;
    unsigned short it = 0;
    unsigned short adc = 0;
    unsigned short edge = 0;
  };

  using vec_dVerbose = std::vector<std::vector<std::pair<int, int>>>;

  // Neural network parameters and modules
  bool gen_hits = false;
  bool use_nn = false;
  const int nd = 5;
  torch::jit::script::Module module_pos;

  struct thread_data
  {
    PHG4TpcCylinderGeom *layergeom = nullptr;
    TrkrHitSet *hitset = nullptr;
    RawHitSet *rawhitset = nullptr;
    ActsGeometry *tGeometry = nullptr;
    unsigned int layer = 0;
    int side = 0;
    unsigned int sector = 0;
    float radius = 0;
    float drift_velocity = 0;
    unsigned short pads_per_sector = 0;
    float phistep = 0;
    float pedestal = 0;
    float seed_threshold = 0;
    float edge_threshold = 0;
    float min_err_squared = 0;
    float min_clus_size = 0;
    float min_adc_sum = 0;
    bool do_assoc = true;
    bool do_wedge_emulation = true;
    bool do_singles = true;
    bool do_split = true;
    int FixedWindow = 0;
    unsigned short phibins = 0;
    unsigned short phioffset = 0;
    unsigned short tbins = 0;
    unsigned short toffset = 0;
    unsigned short maxHalfSizeT = 0;
    unsigned short maxHalfSizePhi = 0;
    double m_tdriftmax = 0;
    double sampa_tbias = 0;
    std::vector<assoc> association_vector;
    std::vector<TrkrCluster *> cluster_vector;
    std::vector<TrainingHits *> v_hits;
    int verbosity = 0;
    bool fillClusHitsVerbose = false;
    vec_dVerbose phivec_ClusHitsVerbose;  // only fill if fillClusHitsVerbose
    vec_dVerbose zvec_ClusHitsVerbose;    // only fill if fillClusHitsVerbose
  };

  pthread_mutex_t mythreadlock;

  void remove_hit(double adc, int phibin, int tbin, int edge, std::multimap<unsigned short, ihit> &all_hit_map, std::vector<std::vector<unsigned short>> &adcval)
  {
    using hit_iterator = std::multimap<unsigned short, ihit>::iterator;
    std::pair<hit_iterator, hit_iterator> iterpair = all_hit_map.equal_range(adc);
    hit_iterator it = iterpair.first;
    for (; it != iterpair.second; ++it)
    {
      if (it->second.iphi == phibin && it->second.it == tbin)
      {
        all_hit_map.erase(it);
        break;
      }
    }
    if (edge)
    {
      adcval[phibin][tbin] = USHRT_MAX;
    }
    else
    {
      adcval[phibin][tbin] = 0;
    }
  }

  void remove_hits(std::vector<ihit> &ihit_list, std::multimap<unsigned short, ihit> &all_hit_map, std::vector<std::vector<unsigned short>> &adcval)
  {
    for (auto &iter : ihit_list)
    {
      unsigned short adc = iter.adc;
      unsigned short phibin = iter.iphi;
      unsigned short tbin = iter.it;
      unsigned short edge = iter.edge;
      remove_hit(adc, phibin, tbin, edge, all_hit_map, adcval);
    }
  }

  void find_t_range(int phibin, int tbin, const thread_data &my_data, const std::vector<std::vector<unsigned short>> &adcval, int &tdown, int &tup, int &touch, int &edge)
  {
    const int FitRangeT = (int) my_data.maxHalfSizeT;
    const int NTBinsMax = (int) my_data.tbins;
    const int FixedWindow = (int) my_data.FixedWindow;
    tup = 0;
    tdown = 0;
    if (FixedWindow != 0)
    {
      tup = FixedWindow;
      tdown = FixedWindow;
      if (tbin + tup >= NTBinsMax)
      {
        tup = NTBinsMax - tbin - 1;
        edge++;
      }
      if ((tbin - tdown) <= 0)
      {
        tdown = tbin;
        edge++;
      }
      return;
    }
    for (int it = 0; it < FitRangeT; it++)
    {
      int ct = tbin + it;

      if (ct <= 0 || ct >= NTBinsMax)
      {
        // tup = it;
        edge++;
        break;  // truncate edge
      }

      if (adcval[phibin][ct] <= 0)
      {
        break;
      }
      if (adcval[phibin][ct] == USHRT_MAX)
      {
        touch++;
        break;
      }
      if (my_data.do_split)
      {
        // check local minima and break at minimum.
        if (ct < NTBinsMax - 4)
        {  // make sure we stay clear from the edge
          if (adcval[phibin][ct] + adcval[phibin][ct + 1] <
              adcval[phibin][ct + 2] + adcval[phibin][ct + 3])
          {  // rising again
            tup = it + 1;
            touch++;
            break;
          }
        }
      }
      tup = it;
    }
    for (int it = 0; it < FitRangeT; it++)
    {
      int ct = tbin - it;
      if (ct <= 0 || ct >= NTBinsMax)
      {
        //      tdown = it;
        edge++;
        break;  // truncate edge
      }
      if (adcval[phibin][ct] <= 0)
      {
        break;
      }
      if (adcval[phibin][ct] == USHRT_MAX)
      {
        touch++;
        break;
      }
      if (my_data.do_split)
      {  // check local minima and break at minimum.
        if (ct > 4)
        {  // make sure we stay clear from the edge
          if (adcval[phibin][ct] + adcval[phibin][ct - 1] <
              adcval[phibin][ct - 2] + adcval[phibin][ct - 3])
          {  // rising again
            tdown = it + 1;
            touch++;
            break;
          }
        }
      }
      tdown = it;
    }
    return;
  }

  void find_phi_range(int phibin, int tbin, const thread_data &my_data, const std::vector<std::vector<unsigned short>> &adcval, int &phidown, int &phiup, int &touch, int &edge)
  {
    int FitRangePHI = (int) my_data.maxHalfSizePhi;
    int NPhiBinsMax = (int) my_data.phibins;
    const int FixedWindow = (int) my_data.FixedWindow;
    phidown = 0;
    phiup = 0;
    if (FixedWindow != 0)
    {
      phiup = FixedWindow;
      phidown = FixedWindow;
      if (phibin + phiup >= NPhiBinsMax)
      {
        phiup = NPhiBinsMax - phibin - 1;
        edge++;
      }
      if (phibin - phidown <= 0)
      {
        phidown = phibin;
        edge++;
      }
      return;
    }
    for (int iphi = 0; iphi < FitRangePHI; iphi++)
    {
      int cphi = phibin + iphi;
      if (cphi < 0 || cphi >= NPhiBinsMax)
      {
        // phiup = iphi;
        edge++;
        break;  // truncate edge
      }

      // break when below minimum
      if (adcval[cphi][tbin] <= 0)
      {
        // phiup = iphi;
        break;
      }
      if (adcval[cphi][tbin] == USHRT_MAX)
      {
        touch++;
        break;
      }
      if (my_data.do_split)
      {  // check local minima and break at minimum.
        if (cphi < NPhiBinsMax - 4)
        {  // make sure we stay clear from the edge
          if (adcval[cphi][tbin] + adcval[cphi + 1][tbin] <
              adcval[cphi + 2][tbin] + adcval[cphi + 3][tbin])
          {  // rising again
            phiup = iphi + 1;
            touch++;
            break;
          }
        }
      }
      phiup = iphi;
    }

    for (int iphi = 0; iphi < FitRangePHI; iphi++)
    {
      int cphi = phibin - iphi;
      if (cphi < 0 || cphi >= NPhiBinsMax)
      {
        // phidown = iphi;
        edge++;
        break;  // truncate edge
      }

      if (adcval[cphi][tbin] <= 0)
      {
        // phidown = iphi;
        break;
      }
      if (adcval[cphi][tbin] == USHRT_MAX)
      {
        touch++;
        break;
      }
      if (my_data.do_split)
      {  // check local minima and break at minimum.
        if (cphi > 4)
        {  // make sure we stay clear from the edge
          if (adcval[cphi][tbin] + adcval[cphi - 1][tbin] <
              adcval[cphi - 2][tbin] + adcval[cphi - 3][tbin])
          {  // rising again
            phidown = iphi + 1;
            touch++;
            break;
          }
        }
      }
      phidown = iphi;
    }
    return;
  }

  int is_hit_isolated(int iphi, int it, int NPhiBinsMax, int NTBinsMax, const std::vector<std::vector<unsigned short>> &adcval)
  {
    // check isolated hits
    //  const int NPhiBinsMax = (int) my_data.phibins;
    // const int NTBinsMax = (int) my_data.tbins;

    int isosum = 0;
    int isophimin = iphi - 1;
    if (isophimin < 0)
    {
      isophimin = 0;
    }
    int isophimax = iphi + 1;
    if (!(isophimax < NPhiBinsMax))
    {
      isophimax = NPhiBinsMax - 1;
    }
    int isotmin = it - 1;
    if (isotmin < 0)
    {
      isotmin = 0;
    }
    int isotmax = it + 1;
    if (!(isotmax < NTBinsMax))
    {
      isotmax = NTBinsMax - 1;
    }
    for (int isophi = isophimin; isophi <= isophimax; isophi++)
    {
      for (int isot = isotmin; isot <= isotmax; isot++)
      {
        if (isophi == iphi && isot == it)
        {
          continue;
        }
        /*
        std::cout <<" isominphi: " << isophimin
                  << "phi: " << isophi
                  <<" isomaxphi: " << isophimax
                  <<" maxphi: " << NPhiBinsMax
                  <<" isomint: " << isotmin
                  << " t: " << isot
                  <<" isomaxt: " << isotmax
                  << " maxt: " << NTBinsMax
                  << std::endl;
        */
        isosum += adcval[isophi][isot];
        /*
        std::cout << "adc " << adcval[isophi][isot]
                  << " isosum " << isosum
                  << " done"
                  << std::endl;
        */
      }
    }
    int isiso = 0;
    if (isosum == 0)
    {
      isiso = 1;
    }
    return isiso;
  }

  void get_cluster(int phibin, int tbin, const thread_data &my_data, const std::vector<std::vector<unsigned short>> &adcval, std::vector<ihit> &ihit_list, int &touch, int &edge)
  {
    // search along phi at the peak in t
    //    const int NPhiBinsMax = (int) my_data.phibins;
    // const int NTBinsMax = (int) my_data.tbins;
    int tup = 0;
    int tdown = 0;
    find_t_range(phibin, tbin, my_data, adcval, tdown, tup, touch, edge);
    // now we have the t extent of the cluster, go find the phi edges
    for (int it = tbin - tdown; it <= (tbin + tup); it++)
    {
      int phiup = 0;
      int phidown = 0;
      find_phi_range(phibin, it, my_data, adcval, phidown, phiup, touch, edge);
      for (int iphi = (phibin - phidown); iphi <= (phibin + phiup); iphi++)
      {
        if (adcval[iphi][it] > 0 && adcval[iphi][it] != USHRT_MAX)
        {
          if (my_data.do_singles)
          {
            if (is_hit_isolated(iphi, it, (int) my_data.phibins, (int) my_data.tbins, adcval))
            {
              continue;
            }
          }
          ihit hit;
          hit.iphi = iphi;
          hit.it = it;

          hit.adc = adcval[iphi][it];
          if (touch > 0)
          {
            if ((iphi == (phibin - phidown)) ||
                (iphi == (phibin + phiup)))
            {
              hit.edge = 1;
            }
          }
          ihit_list.push_back(hit);
        }
      }
    }
    return;
  }

  void calc_cluster_parameter(const int iphi_center, const int it_center,
                              const std::vector<ihit> &ihit_list, thread_data &my_data, int ntouch, int nedge)
  {
    //
    // get z range from layer geometry
    /* these are used for rescaling the drift velocity */
    // const double z_min = -105.5;
    // const double z_max = 105.5;
    //  std::cout << "calc clus" << std::endl;
    //  loop over the hits in this cluster
    double t_sum = 0.0;
    // double phi_sum = 0.0;
    double adc_sum = 0.0;
    double t2_sum = 0.0;
    // double phi2_sum = 0.0;

    double iphi_sum = 0.0;
    double iphi2_sum = 0.0;

    double radius = my_data.layergeom->get_radius();  // returns center of layer

    int phibinhi = -1;
    int phibinlo = 666666;
    int tbinhi = -1;
    int tbinlo = 666666;
    int clus_size = ihit_list.size();
    int max_adc = 0;
    if (clus_size <= my_data.min_clus_size)
    {
      return;
    }

    // training information
    TrainingHits *training_hits = nullptr;
    if (gen_hits)
    {
      training_hits = new TrainingHits;
      assert(training_hits);
      training_hits->radius = radius;
      training_hits->phi = my_data.layergeom->get_phicenter(iphi_center + my_data.phioffset);
      double center_t = my_data.layergeom->get_zcenter(it_center + my_data.toffset) + my_data.sampa_tbias;
      training_hits->z = (my_data.m_tdriftmax - center_t) * my_data.tGeometry->get_drift_velocity();
      if (my_data.side == 0)
      {
        training_hits->z = -training_hits->z;
      }
      training_hits->phistep = my_data.layergeom->get_phistep();
      training_hits->zstep = my_data.layergeom->get_zstep() * my_data.tGeometry->get_drift_velocity();
      training_hits->layer = my_data.layer;
      training_hits->ntouch = ntouch;
      training_hits->nedge = nedge;
      training_hits->v_adc.fill(0);
    }

    //      std::cout << "process list" << std::endl;
    std::vector<TrkrDefs::hitkey> hitkeyvec;

    // keep track of the hit locations in a given cluster
    std::map<int, unsigned int> m_phi{};
    std::map<int, unsigned int> m_z{};

    for (auto iter : ihit_list)
    {
      int iphi = iter.iphi + my_data.phioffset;
      int it = iter.it + my_data.toffset;
      double adc = iter.adc;

      if (adc <= 0)
      {
        continue;
      }

      if (adc > max_adc)
      {
        max_adc = adc;
      }

      if (iphi > phibinhi)
      {
        phibinhi = iphi;
      }

      if (iphi < phibinlo)
      {
        phibinlo = iphi;
      }

      if (it > tbinhi)
      {
        tbinhi = it;
      }

      if (it < tbinlo)
      {
        tbinlo = it;
      }

      // if(it==it_center){ yg_sum += adc; }
      // update phi sums
      //	double phi_center = my_data.layergeom->get_phicenter(iphi);

      // phi_sum += phi_center * adc;
      // phi2_sum += square(phi_center)*adc;
      //	std::cout << "phi_center: " << phi_center << " adc: " << adc <<std::endl;
      iphi_sum += iphi * adc;
      iphi2_sum += square(iphi) * adc;

      // update t sums
      double t = my_data.layergeom->get_zcenter(it);
      t_sum += t * adc;
      t2_sum += square(t) * adc;

      adc_sum += adc;

      if (my_data.fillClusHitsVerbose)
      {
        auto pnew = m_phi.try_emplace(iphi, adc);
        if (!pnew.second)
        {
          pnew.first->second += adc;
        }

        pnew = m_z.try_emplace(it, adc);
        if (!pnew.second)
        {
          pnew.first->second += adc;
        }
      }

      // capture the hitkeys for all adc values above a certain threshold
      TrkrDefs::hitkey hitkey = TpcDefs::genHitKey(iphi, it);
      // if(adc>5)
      hitkeyvec.push_back(hitkey);

      // training adc
      if (gen_hits && training_hits)
      {
        int iphi_diff = iter.iphi - iphi_center;
        int it_diff = iter.it - it_center;
        if (std::abs(iphi_diff) <= nd && std::abs(it_diff) <= nd)
        {
          training_hits->v_adc[(iphi_diff + nd) * (2 * nd + 1) + (it_diff + nd)] = adc;
        }
      }
    }
    //      std::cout << "done process list" << std::endl;
    if (adc_sum < my_data.min_adc_sum)
    {
      hitkeyvec.clear();
      return;  // skip obvious noise "clusters"
    }

    // This is the global position
    double clusiphi = iphi_sum / adc_sum;
    double clusphi = my_data.layergeom->get_phi(clusiphi);

    float clusx = radius * cos(clusphi);
    float clusy = radius * sin(clusphi);
    double clust = t_sum / adc_sum;
    // needed for surface identification
    double zdriftlength = clust * my_data.tGeometry->get_drift_velocity();
    // convert z drift length to z position in the TPC
    double clusz = my_data.m_tdriftmax * my_data.tGeometry->get_drift_velocity() - zdriftlength;
    if (my_data.side == 0)
    {
      clusz = -clusz;
    }
    // std::cout << " side " << my_data.side << " clusz " << clusz << " clust " << clust << " driftmax " << my_data.m_tdriftmax << std::endl;
    const double phi_cov = (iphi2_sum / adc_sum - square(clusiphi)) * pow(my_data.layergeom->get_phistep(), 2);
    const double t_cov = t2_sum / adc_sum - square(clust);

    // Get the surface key to find the surface from the
    TrkrDefs::hitsetkey tpcHitSetKey = TpcDefs::genHitSetKey(my_data.layer, my_data.sector, my_data.side);
    Acts::Vector3 global(clusx, clusy, clusz);
    TrkrDefs::subsurfkey subsurfkey = 0;

    Surface surface = my_data.tGeometry->get_tpc_surface_from_coords(
        tpcHitSetKey,
        global,
        subsurfkey);

    if (!surface)
    {
      /// If the surface can't be found, we can't track with it. So
      /// just return and don't add the cluster to the container
      hitkeyvec.clear();
      return;
    }

    // Estimate the errors
    const double phi_err_square = (phibinhi == phibinlo) ? square(radius * my_data.layergeom->get_phistep()) / 12 : square(radius) * phi_cov / (adc_sum * 0.14);

    const double t_err_square = (tbinhi == tbinlo) ? square(my_data.layergeom->get_zstep()) / 12 : t_cov / (adc_sum * 0.14);

    char tsize = tbinhi - tbinlo + 1;
    char phisize = phibinhi - phibinlo + 1;
    // std::cout << "phisize: "  << (int) phisize << " phibinhi " << phibinhi << " phibinlo " << phibinlo << std::endl;
    // phi_cov = (weighted mean of dphi^2) - (weighted mean of dphi)^2,  which is essentially the weighted mean of dphi^2. The error is then:
    // e_phi = sigma_dphi/sqrt(N) = sqrt( sigma_dphi^2 / N )  -- where N is the number of samples of the distribution with standard deviation sigma_dphi
    //    - N is the number of electrons that drift to the readout plane
    // We have to convert (sum of adc units for all bins in the cluster) to number of ionization electrons N
    // Conversion gain is 20 mV/fC - relates total charge collected on pad to PEAK voltage out of ADC. The GEM gain is assumed to be 2000
    // To get equivalent charge per T bin, so that summing ADC input voltage over all T bins returns total input charge, divide voltages by 2.4 for 80 ns SAMPA
    // Equivalent charge per T bin is then  (ADU x 2200 mV / 1024) / 2.4 x (1/20) fC/mV x (1/1.6e-04) electrons/fC x (1/2000) = ADU x 0.14

    // SAMPA shaping bias correction
    clust = clust + my_data.sampa_tbias;

    /// convert to Acts units
    global *= Acts::UnitConstants::cm;
    // std::cout << "transform" << std::endl;
    Acts::Vector3 local = surface->transform(my_data.tGeometry->geometry().getGeoContext()).inverse() * global;
    local /= Acts::UnitConstants::cm;
    // std::cout << "done transform" << std::endl;
    //  we need the cluster key and all associated hit keys (note: the cluster key includes the hitset key)

    TrkrCluster *clus_base = nullptr;
    bool b_made_cluster{false};

    // std::cout << "ver5" << std::endl;
    //	std::cout << "clus num" << my_data.cluster_vector.size() << " X " << local(0) << " Y " << clust << std::endl;
    if (sqrt(phi_err_square) > my_data.min_err_squared)
    {
      auto clus = new TrkrClusterv5;
      // auto clus = std::make_unique<TrkrClusterv3>();
      clus_base = clus;
      clus->setAdc(adc_sum);
      clus->setMaxAdc(max_adc);
      clus->setEdge(nedge);
      clus->setPhiSize(phisize);
      clus->setZSize(tsize);
      clus->setSubSurfKey(subsurfkey);
      clus->setOverlap(ntouch);
      clus->setLocalX(local(0));
      clus->setLocalY(clust);
      clus->setPhiError(sqrt(phi_err_square));
      clus->setZError(sqrt(t_err_square * pow(my_data.tGeometry->get_drift_velocity(), 2)));
      my_data.cluster_vector.push_back(clus);
      b_made_cluster = true;
    }

    if (use_nn && clus_base && training_hits)
    {
      try
      {
        // Create a vector of inputs
        std::vector<torch::jit::IValue> inputs;
        inputs.emplace_back(torch::stack({torch::from_blob(std::vector<float>(training_hits->v_adc.begin(), training_hits->v_adc.end()).data(), {1, 2 * nd + 1, 2 * nd + 1}, torch::kFloat32),
                                          torch::full({1, 2 * nd + 1, 2 * nd + 1}, std::clamp((training_hits->layer - 7) / 16, 0, 2), torch::kFloat32),
                                          torch::full({1, 2 * nd + 1, 2 * nd + 1}, training_hits->z / radius, torch::kFloat32)},
                                         1));

        // Execute the model and turn its output into a tensor
        at::Tensor ten_pos = module_pos.forward(inputs).toTensor();
        float nn_phi = training_hits->phi + std::clamp(ten_pos[0][0][0].item<float>(), -(float) nd, (float) nd) * training_hits->phistep;
        float nn_z = training_hits->z + std::clamp(ten_pos[0][1][0].item<float>(), -(float) nd, (float) nd) * training_hits->zstep;
        float nn_x = radius * std::cos(nn_phi);
        float nn_y = radius * std::sin(nn_phi);
        Acts::Vector3 nn_global(nn_x, nn_y, nn_z);
        nn_global *= Acts::UnitConstants::cm;
        Acts::Vector3 nn_local = surface->transform(my_data.tGeometry->geometry().geoContext).inverse() * nn_global;
        nn_local /= Acts::UnitConstants::cm;
        float nn_t = my_data.m_tdriftmax - std::fabs(nn_z) / my_data.tGeometry->get_drift_velocity();
        clus_base->setLocalX(nn_local(0));
        clus_base->setLocalY(nn_t);
      }
      catch (const c10::Error &e)
      {
        std::cout << PHWHERE << "Error: Failed to execute NN modules" << std::endl;
      }
    }  // use_nn

    if (my_data.fillClusHitsVerbose && b_made_cluster)
    {
      // push the data back to
      my_data.phivec_ClusHitsVerbose.push_back(std::vector<std::pair<int, int>>{});
      my_data.zvec_ClusHitsVerbose.push_back(std::vector<std::pair<int, int>>{});

      auto &vphi = my_data.phivec_ClusHitsVerbose.back();
      auto &vz = my_data.zvec_ClusHitsVerbose.back();

      for (auto &entry : m_phi)
      {
        vphi.push_back({entry.first, entry.second});
      }
      for (auto &entry : m_z)
      {
        vz.push_back({entry.first, entry.second});
      }
    }

    // std::cout << "end clus out" << std::endl;
    //       if(my_data.do_assoc && my_data.clusterhitassoc){
    if (my_data.do_assoc)
    {
      // get cluster index in vector. It is used to store associations, and build relevant cluster keys when filling the containers
      uint32_t index = my_data.cluster_vector.size() - 1;
      for (unsigned int &i : hitkeyvec)
      {
        my_data.association_vector.emplace_back(index, i);
      }
      if (gen_hits && training_hits)
      {
        training_hits->cluskey = TrkrDefs::genClusKey(tpcHitSetKey, index);
      }
    }
    hitkeyvec.clear();
    if (gen_hits && training_hits)
    {
      my_data.v_hits.emplace_back(training_hits);
    }
    //      std::cout << "done calc" << std::endl;
  }

  void ProcessSectorData(thread_data *my_data)
  {
    const auto &pedestal = my_data->pedestal;
    const auto &phibins = my_data->phibins;
    const auto &phioffset = my_data->phioffset;
    const auto &tbins = my_data->tbins;
    const auto &toffset = my_data->toffset;
    const auto &layer = my_data->layer;
    //    int nhits = 0;
    // for convenience, create a 2D vector to store adc values in and initialize to zero
    std::vector<std::vector<unsigned short>> adcval(phibins, std::vector<unsigned short>(tbins, 0));
    std::multimap<unsigned short, ihit> all_hit_map;
    std::vector<ihit> hit_vect;

    int tbinmax = tbins;
    int tbinmin = 0;
    if (my_data->do_wedge_emulation)
    {
      if (layer >= 7 && layer < 22)
      {
        int etacut = (tbins / 2.) - ((50 + (layer - 7)) / 105.5) * (tbins / 2.);
        tbinmin = etacut;
        tbinmax -= etacut;
      }
      if (layer >= 22 && layer <= 48)
      {
        int etacut = (tbins / 2.) - ((65 + ((40.5 / 26) * (layer - 22))) / 105.5) * (tbins / 2.);
        tbinmin = etacut;
        tbinmax -= etacut;
      }
    }

    if (my_data->hitset != nullptr)
    {
      TrkrHitSet *hitset = my_data->hitset;
      TrkrHitSet::ConstRange hitrangei = hitset->getHits();

      for (TrkrHitSet::ConstIterator hitr = hitrangei.first;
           hitr != hitrangei.second;
           ++hitr)
      {
        if (TpcDefs::getPad(hitr->first) - phioffset < 0)
        {
          // std::cout << "WARNING phibin out of range: " << TpcDefs::getPad(hitr->first) - phioffset << " | " << phibins << std::endl;
          continue;
        }
        if (TpcDefs::getTBin(hitr->first) - toffset < 0)
        {
          // std::cout << "WARNING tbin out of range: " << TpcDefs::getTBin(hitr->first) - toffset  << " | " << tbins <<std::endl;
        }
        unsigned short phibin = TpcDefs::getPad(hitr->first) - phioffset;
        unsigned short tbin = TpcDefs::getTBin(hitr->first) - toffset;
        unsigned short tbinorg = TpcDefs::getTBin(hitr->first);
        if (phibin >= phibins)
        {
          // std::cout << "WARNING phibin out of range: " << phibin << " | " << phibins << std::endl;
          continue;
        }
        if (tbin >= tbins)
        {
          // std::cout << "WARNING z bin out of range: " << tbin << " | " << tbins << std::endl;
          continue;
        }
        if (tbinorg > tbinmax || tbinorg < tbinmin)
        {
          continue;
        }
        float_t fadc = (hitr->second->getAdc()) - pedestal;  // proper int rounding +0.5
        unsigned short adc = 0;
        if (fadc > 0)
        {
          adc = (unsigned short) fadc;
        }
        if (phibin >= phibins)
        {
          continue;
        }
        if (tbin >= tbins)
        {
          continue;  // tbin is unsigned int, <0 cannot happen
        }

        if (adc > 0)
        {
          if (adc > (my_data->seed_threshold))
          {
            ihit thisHit;

            thisHit.iphi = phibin;
            thisHit.it = tbin;
            thisHit.adc = adc;
            thisHit.edge = 0;
            all_hit_map.insert(std::make_pair(adc, thisHit));
          }
          if (adc > my_data->edge_threshold)
          {
            adcval[phibin][tbin] = (unsigned short) adc;
          }
        }
      }
    }
    else if (my_data->rawhitset != nullptr)
    {
      RawHitSet *hitset = my_data->rawhitset;
      /*std::cout << "Layer: " << my_data->layer
                << "Side: " << my_data->side
                << "Sector: " << my_data->sector
                << " nhits:  " << hitset.size()
                << std::endl;
      */
      for (int nphi = 0; nphi < phibins; nphi++)
      {
        if (hitset->size(nphi) == 0)
        {
          continue;
        }

        int pindex = 0;
        for (unsigned int nt = 0; nt < hitset->size(nphi); nt++)
        {
          unsigned short val = (*(hitset->getHits(nphi)))[nt];

          if (val == 0)
          {
            pindex++;
          }
          else
          {
            if (nt == 0)
            {
              if (val > 5)
              {
                ihit thisHit;
                thisHit.iphi = nphi;
                thisHit.it = pindex;
                thisHit.adc = val;
                thisHit.edge = 0;
                all_hit_map.insert(std::make_pair(val, thisHit));
              }
              adcval[nphi][pindex++] = val;
            }
            else
            {
              if (((*(hitset->getHits(nphi)))[nt - 1] == 0) && ((*(hitset->getHits(nphi)))[nt + 1] == 0))
              {  // found zero count
                pindex += val;
              }
              else
              {
                if (val > 5)
                {
                  ihit thisHit;
                  thisHit.iphi = nphi;
                  thisHit.it = pindex;
                  thisHit.adc = val;
                  thisHit.edge = 0;
                  all_hit_map.insert(std::make_pair(val, thisHit));
                }
                adcval[nphi][pindex++] = val;
              }
            }
          }
        }
      }
    }
    /*
    if (my_data->do_singles)
    {
      for (auto ahit : all_hit_map)
      {
        ihit hiHit = ahit.second;
        int iphi = hiHit.iphi;
        int it = hiHit.it;
        unsigned short edge = hiHit.edge;
        double adc = hiHit.adc;
        if (it > 0 && it < tbins)
        {
          if (adcval[iphi][it - 1] == 0 &&
              adcval[iphi][it + 1] == 0)
          {
            remove_hit(adc, iphi, it, edge, all_hit_map, adcval);
          }
        }
      }
    }
    */
    // std::cout << "done filling " << std::endl;
    while (all_hit_map.size() > 0)
    {
      // std::cout << "all hit map size: " << all_hit_map.size() << std::endl;
      auto iter = all_hit_map.rbegin();
      if (iter == all_hit_map.rend())
      {
        break;
      }
      ihit hiHit = iter->second;
      int iphi = hiHit.iphi;
      int it = hiHit.it;
      unsigned short edge = hiHit.edge;
      double adc = hiHit.adc;
      if (my_data->do_singles)
      {
        if (is_hit_isolated(iphi, it, (int) my_data->phibins, (int) my_data->tbins, adcval))
        {
          remove_hit(adc, iphi, it, edge, all_hit_map, adcval);
          continue;
        }
      }

      // put all hits in the all_hit_map (sorted by adc)
      // start with highest adc hit
      //  -> cluster around it and get vector of hits
      std::vector<ihit> ihit_list;
      int ntouch = 0;
      int nedge = 0;
      get_cluster(iphi, it, *my_data, adcval, ihit_list, ntouch, nedge);

      if (my_data->FixedWindow > 0)
      {
        // get cluster dimensions and check if they hit the window boundaries
        int wphibinhi = -1;
        int wphibinlo = 666666;
        int wtbinhi = -1;
        int wtbinlo = 666666;
        for (auto witer : ihit_list)
        {
          int wiphi = witer.iphi + my_data->phioffset;
          int wit = witer.it + my_data->toffset;
          double wadc = witer.adc;
          if (wadc <= 0)
          {
            continue;
          }
          if (wiphi > wphibinhi)
          {
            wphibinhi = wiphi;
          }
          if (wiphi < wphibinlo)
          {
            wphibinlo = wiphi;
          }
          if (wit > wtbinhi)
          {
            wtbinhi = wit;
          }
          if (wit < wtbinlo)
          {
            wtbinlo = wit;
          }
        }
        char wtsize = wtbinhi - wtbinlo + 1;
        char wphisize = wphibinhi - wphibinlo + 1;

        // check if we have a super big cluster and switch from fixed window to step down
        if (ihit_list.size() > (0.5 * pow((2 * my_data->FixedWindow + 1), 2)) ||
            wphisize >= (2 * my_data->FixedWindow + 1) ||
            wtsize >= (2 * my_data->FixedWindow + 1))
        {
          int window_cache = my_data->FixedWindow;
          // std::cout << " fixed size before " << ihit_list.size() << " | " << pow((2*my_data->FixedWindow+1),2) << std::endl;
          my_data->FixedWindow = 0;
          // reset hit list and try again without fixed window
          ihit_list.clear();
          get_cluster(iphi, it, *my_data, adcval, ihit_list, ntouch, nedge);
          // std::cout << " stepdown size after " << ihit_list.size() << std::endl;
          my_data->FixedWindow = window_cache;
        }
      }
      if (ihit_list.size() <= 1)
      {
        remove_hits(ihit_list, all_hit_map, adcval);
        ihit_list.clear();
        remove_hit(adc, iphi, it, edge, all_hit_map, adcval);
      }
      // -> calculate cluster parameters
      // -> add hits to truth association
      // remove hits from all_hit_map
      // repeat untill all_hit_map empty
      calc_cluster_parameter(iphi, it, ihit_list, *my_data, ntouch, nedge);
      remove_hits(ihit_list, all_hit_map, adcval);
      ihit_list.clear();
    }
    /*    if( my_data->rawhitset!=nullptr){
      RawHitSetv1 *hitset = my_data->rawhitset;
      std::cout << "Layer: " << my_data->layer
                << " Side: " << my_data->side
                << " Sector: " << my_data->sector
                << " nhits:  " << hitset->size()
                << " nhits coutn :  " << nhits
                << " nclus: " << my_data->cluster_vector.size()
                << std::endl;
    }
    */
    // pthread_exit(nullptr);
  }
  void *ProcessSector(void *threadarg)
  {
    auto my_data = static_cast<thread_data *>(threadarg);
    ProcessSectorData(my_data);
    pthread_exit(nullptr);
  }
}  // namespace

TpcClusterizer::TpcClusterizer(const std::string &name)
  : SubsysReco(name)
  , m_training(nullptr)
{
}

bool TpcClusterizer::is_in_sector_boundary(int phibin, int sector, PHG4TpcCylinderGeom *layergeom) const
{
  bool reject_it = false;

  // sector boundaries occur every 1/12 of the full phi bin range
  int PhiBins = layergeom->get_phibins();
  int PhiBinsSector = PhiBins / 12;

  double radius = layergeom->get_radius();
  double PhiBinSize = 2.0 * radius * M_PI / (double) PhiBins;

  // sector starts where?
  int sector_lo = sector * PhiBinsSector;
  int sector_hi = sector_lo + PhiBinsSector - 1;

  int sector_fiducial_bins = (int) (SectorFiducialCut / PhiBinSize);

  if (phibin < sector_lo + sector_fiducial_bins || phibin > sector_hi - sector_fiducial_bins)
  {
    reject_it = true;
    /*
    int layer = layergeom->get_layer();
    std::cout << " local maximum is in sector fiducial boundary: layer " << layer << " radius " << radius << " sector " << sector
    << " PhiBins " << PhiBins << " sector_fiducial_bins " << sector_fiducial_bins
    << " PhiBinSize " << PhiBinSize << " phibin " << phibin << " sector_lo " << sector_lo << " sector_hi " << sector_hi << std::endl;
    */
  }

  return reject_it;
}

int TpcClusterizer::InitRun(PHCompositeNode *topNode)
{
  PHNodeIterator iter(topNode);

  // Looking for the DST node
  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    std::cout << PHWHERE << "DST Node missing, doing nothing." << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  // Create the Cluster node if required
  auto trkrclusters = findNode::getClass<TrkrClusterContainer>(dstNode, "TRKR_CLUSTER");
  if (!trkrclusters)
  {
    PHNodeIterator dstiter(dstNode);
    PHCompositeNode *DetNode =
        dynamic_cast<PHCompositeNode *>(dstiter.findFirst("PHCompositeNode", "TRKR"));
    if (!DetNode)
    {
      DetNode = new PHCompositeNode("TRKR");
      dstNode->addNode(DetNode);
    }

    trkrclusters = new TrkrClusterContainerv4;
    PHIODataNode<PHObject> *TrkrClusterContainerNode =
        new PHIODataNode<PHObject>(trkrclusters, "TRKR_CLUSTER", "PHObject");
    DetNode->addNode(TrkrClusterContainerNode);
  }

  auto clusterhitassoc = findNode::getClass<TrkrClusterHitAssoc>(topNode, "TRKR_CLUSTERHITASSOC");
  if (!clusterhitassoc)
  {
    PHNodeIterator dstiter(dstNode);
    PHCompositeNode *DetNode =
        dynamic_cast<PHCompositeNode *>(dstiter.findFirst("PHCompositeNode", "TRKR"));
    if (!DetNode)
    {
      DetNode = new PHCompositeNode("TRKR");
      dstNode->addNode(DetNode);
    }

    clusterhitassoc = new TrkrClusterHitAssocv3;
    PHIODataNode<PHObject> *newNode = new PHIODataNode<PHObject>(clusterhitassoc, "TRKR_CLUSTERHITASSOC", "PHObject");
    DetNode->addNode(newNode);
  }

  auto training_container = findNode::getClass<TrainingHitsContainer>(dstNode, "TRAINING_HITSET");
  if (!training_container)
  {
    PHNodeIterator dstiter(dstNode);
    PHCompositeNode *DetNode =
        dynamic_cast<PHCompositeNode *>(dstiter.findFirst("PHCompositeNode", "TRKR"));
    if (!DetNode)
    {
      DetNode = new PHCompositeNode("TRKR");
      dstNode->addNode(DetNode);
    }

    training_container = new TrainingHitsContainer;
    PHIODataNode<PHObject> *TrainingHitsContainerNode =
        new PHIODataNode<PHObject>(training_container, "TRAINING_HITSET", "PHObject");
    DetNode->addNode(TrainingHitsContainerNode);
  }

  gen_hits = _store_hits || _use_nn;
  use_nn = _use_nn;
  if (use_nn)
  {
    const char *offline_main = std::getenv("OFFLINE_MAIN");
    assert(offline_main);
    std::string net_model = std::string(offline_main) + "/share/tpc/net_model.pt";
    try
    {
      // Deserialize the ScriptModule from a file using torch::jit::load()
      module_pos = torch::jit::load(net_model);
      std::cout << PHWHERE << "Load NN module: " << net_model << std::endl;
    }
    catch (const c10::Error &e)
    {
      std::cout << PHWHERE << "Error: Cannot load module " << net_model << std::endl;
      exit(1);
    }
  }
  else
  {
    std::cout << PHWHERE << "Use traditional clustering" << std::endl;
  }

  if (record_ClusHitsVerbose)
  {
    // get the node
    mClusHitsVerbose = findNode::getClass<ClusHitsVerbosev1>(topNode, "Trkr_SvtxClusHitsVerbose");
    if (!mClusHitsVerbose)
    {
      PHNodeIterator dstiter(dstNode);
      auto DetNode = dynamic_cast<PHCompositeNode *>(dstiter.findFirst("PHCompositeNode", "TRKR"));
      if (!DetNode)
      {
        DetNode = new PHCompositeNode("TRKR");
        dstNode->addNode(DetNode);
      }
      mClusHitsVerbose = new ClusHitsVerbosev1();
      auto newNode = new PHIODataNode<PHObject>(mClusHitsVerbose, "Trkr_SvtxClusHitsVerbose", "PHObject");
      DetNode->addNode(newNode);
    }
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int TpcClusterizer::process_event(PHCompositeNode *topNode)
{
  // The TPC is the only subsystem that clusters in global coordinates. For consistency,
  // we must use the construction transforms to get the local coordinates.
  // Set the flag to use ideal transforms for the duration of this process_event, for thread safety
  alignmentTransformationContainer::use_alignment = false;

  //  int print_layer = 18;

  if (Verbosity() > 1000)
  {
    std::cout << "TpcClusterizer::Process_Event" << std::endl;
  }

  PHNodeIterator iter(topNode);
  PHCompositeNode *dstNode = static_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    std::cout << PHWHERE << "DST Node missing, doing nothing." << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  if (!do_read_raw)
  {
    // get node containing the digitized hits
    m_hits = findNode::getClass<TrkrHitSetContainer>(topNode, "TRKR_HITSET");
    if (!m_hits)
    {
      std::cout << PHWHERE << "ERROR: Can't find node TRKR_HITSET" << std::endl;
      return Fun4AllReturnCodes::ABORTRUN;
    }
  }
  else
  {
    // get node containing the digitized hits
    m_rawhits = findNode::getClass<RawHitSetContainer>(topNode, "TRKR_RAWHITSET");
    if (!m_rawhits)
    {
      std::cout << PHWHERE << "ERROR: Can't find node TRKR_HITSET" << std::endl;
      return Fun4AllReturnCodes::ABORTRUN;
    }
  }

  // get laser event info, if exists and event has laser and rejection is on, don't bother with clustering
  LaserEventInfo *laserInfo = findNode::getClass<LaserEventInfo>(topNode, "LaserEventInfo");
  if (m_rejectEvent && laserInfo && laserInfo->isLaserEvent())
  {
    return Fun4AllReturnCodes::EVENT_OK;
  }

  // get node for clusters
  m_clusterlist = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");
  if (!m_clusterlist)
  {
    std::cout << PHWHERE << " ERROR: Can't find TRKR_CLUSTER." << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  // get node for cluster hit associations
  m_clusterhitassoc = findNode::getClass<TrkrClusterHitAssoc>(topNode, "TRKR_CLUSTERHITASSOC");
  if (!m_clusterhitassoc)
  {
    std::cout << PHWHERE << " ERROR: Can't find TRKR_CLUSTERHITASSOC" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  // get node for training hits
  m_training = findNode::getClass<TrainingHitsContainer>(topNode, "TRAINING_HITSET");
  if (!m_training)
  {
    std::cout << PHWHERE << " ERROR: Can't find TRAINING_HITSET." << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  PHG4TpcCylinderGeomContainer *geom_container =
      findNode::getClass<PHG4TpcCylinderGeomContainer>(topNode, "CYLINDERCELLGEOM_SVTX");
  if (!geom_container)
  {
    std::cout << PHWHERE << "ERROR: Can't find node CYLINDERCELLGEOM_SVTX" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  m_tGeometry = findNode::getClass<ActsGeometry>(topNode,
                                                 "ActsGeometry");
  if (!m_tGeometry)
  {
    std::cout << PHWHERE
              << "ActsGeometry not found on node tree. Exiting"
              << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  // The hits are stored in hitsets, where each hitset contains all hits in a given TPC readout (layer, sector, side), so clusters are confined to a hitset
  // The TPC clustering is more complicated than for the silicon, because we have to deal with overlapping clusters

  TrkrHitSetContainer::ConstRange hitsetrange;
  RawHitSetContainer::ConstRange rawhitsetrange;
  int num_hitsets = 0;

  if (!do_read_raw)
  {
    hitsetrange = m_hits->getHitSets(TrkrDefs::TrkrId::tpcId);
    num_hitsets = std::distance(hitsetrange.first, hitsetrange.second);
  }
  else
  {
    rawhitsetrange = m_rawhits->getHitSets(TrkrDefs::TrkrId::tpcId);
    num_hitsets = std::distance(rawhitsetrange.first, rawhitsetrange.second);
  }

  // create structure to store given thread and associated data
  struct thread_pair_t
  {
    pthread_t thread{};
    thread_data data;
  };

  // create vector of thread pairs and reserve the right size upfront to avoid reallocation
  std::vector<thread_pair_t> threads;
  threads.reserve(num_hitsets);

  pthread_attr_t attr;
  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

  if (pthread_mutex_init(&mythreadlock, nullptr) != 0)
  {
    std::cout << std::endl
              << " mutex init failed" << std::endl;
    return 1;
  }
//  int count = 0;

  if (!do_read_raw)
  {
    for (TrkrHitSetContainer::ConstIterator hitsetitr = hitsetrange.first;
         hitsetitr != hitsetrange.second;
         ++hitsetitr)
    {
      // if(count>0)continue;
      TrkrHitSet *hitset = hitsetitr->second;
      unsigned int layer = TrkrDefs::getLayer(hitsetitr->first);
      int side = TpcDefs::getSide(hitsetitr->first);
      unsigned int sector = TpcDefs::getSectorId(hitsetitr->first);
      PHG4TpcCylinderGeom *layergeom = geom_container->GetLayerCellGeom(layer);

      // instanciate new thread pair, at the end of thread vector
      thread_pair_t &thread_pair = threads.emplace_back();
      if (mClusHitsVerbose)
      {
        thread_pair.data.fillClusHitsVerbose = true;
      };

      thread_pair.data.layergeom = layergeom;
      thread_pair.data.hitset = hitset;
      thread_pair.data.rawhitset = nullptr;
      thread_pair.data.layer = layer;
      thread_pair.data.pedestal = pedestal;
      thread_pair.data.seed_threshold = seed_threshold;
      thread_pair.data.edge_threshold = edge_threshold;
      thread_pair.data.sector = sector;
      thread_pair.data.side = side;
      thread_pair.data.do_assoc = do_hit_assoc;
      thread_pair.data.do_wedge_emulation = do_wedge_emulation;
      thread_pair.data.do_singles = do_singles;
      thread_pair.data.tGeometry = m_tGeometry;
      thread_pair.data.maxHalfSizeT = MaxClusterHalfSizeT;
      thread_pair.data.maxHalfSizePhi = MaxClusterHalfSizePhi;
      thread_pair.data.sampa_tbias = m_sampa_tbias;
      thread_pair.data.verbosity = Verbosity();
      thread_pair.data.do_split = do_split;
      thread_pair.data.FixedWindow = do_fixed_window;
      thread_pair.data.min_err_squared = min_err_squared;
      thread_pair.data.min_clus_size = min_clus_size;
      thread_pair.data.min_adc_sum = min_adc_sum;
      unsigned short NPhiBins = (unsigned short) layergeom->get_phibins();
      unsigned short NPhiBinsSector = NPhiBins / 12;
      unsigned short NTBins = 0;
      if (is_reco)
      {
        NTBins = NZBinsSide;
      }
      else
      {
        NTBins = (unsigned short) layergeom->get_zbins();
      }
      unsigned short NTBinsSide = NTBins;
      unsigned short NTBinsMin = 0;
      unsigned short PhiOffset = NPhiBinsSector * sector;
      unsigned short TOffset = NTBinsMin;

      m_tdriftmax = AdcClockPeriod * NZBinsSide;
      thread_pair.data.m_tdriftmax = m_tdriftmax;

      thread_pair.data.phibins = NPhiBinsSector;
      thread_pair.data.phioffset = PhiOffset;
      thread_pair.data.tbins = NTBinsSide;
      thread_pair.data.toffset = TOffset;

      thread_pair.data.radius = layergeom->get_radius();
      thread_pair.data.drift_velocity = m_tGeometry->get_drift_velocity();
      thread_pair.data.pads_per_sector = 0;
      thread_pair.data.phistep = 0;
      int rc;
      rc = pthread_create(&thread_pair.thread, &attr, ProcessSector, (void *) &thread_pair.data);

      if (rc)
      {
        std::cout << "Error:unable to create thread," << rc << std::endl;
      }
      if (do_sequential)
      {
        int rc2 = pthread_join(thread_pair.thread, nullptr);
        if (rc2)
        {
          std::cout << "Error:unable to join," << rc2 << std::endl;
        }

        // get the hitsetkey from thread data
        const auto &data(thread_pair.data);
        const auto hitsetkey = TpcDefs::genHitSetKey(data.layer, data.sector, data.side);

        // copy clusters to map
        for (uint32_t index = 0; index < data.cluster_vector.size(); ++index)
        {
          // generate cluster key
          const auto ckey = TrkrDefs::genClusKey(hitsetkey, index);

          // get cluster
          auto cluster = data.cluster_vector[index];

          // insert in map
          m_clusterlist->addClusterSpecifyKey(ckey, cluster);

          if (mClusHitsVerbose)
          {
            for (auto &hit : data.phivec_ClusHitsVerbose[index])
            {
              mClusHitsVerbose->addPhiHit(hit.first, hit.second);
            }
            for (auto &hit : data.zvec_ClusHitsVerbose[index])
            {
              mClusHitsVerbose->addZHit(hit.first, hit.second);
            }
            mClusHitsVerbose->push_hits(ckey);
          }
        }

        // copy hit associations to map
        for (const auto &[index, hkey] : thread_pair.data.association_vector)
        {
          // generate cluster key
          const auto ckey = TrkrDefs::genClusKey(hitsetkey, index);

          // add to association table
          m_clusterhitassoc->addAssoc(ckey, hkey);
        }
      }
//      count++;
    }
  }
  else
  {
    for (RawHitSetContainer::ConstIterator hitsetitr = rawhitsetrange.first;
         hitsetitr != rawhitsetrange.second;
         ++hitsetitr)
    {
      //	if(count>0)continue;
      //    const auto hitsetid = hitsetitr->first;
      //	std::cout << " starting thread # " << count << std::endl;

      RawHitSet *hitset = hitsetitr->second;
      unsigned int layer = TrkrDefs::getLayer(hitsetitr->first);
      int side = TpcDefs::getSide(hitsetitr->first);
      unsigned int sector = TpcDefs::getSectorId(hitsetitr->first);
      PHG4TpcCylinderGeom *layergeom = geom_container->GetLayerCellGeom(layer);

      // instanciate new thread pair, at the end of thread vector
      thread_pair_t &thread_pair = threads.emplace_back();

      thread_pair.data.layergeom = layergeom;
      thread_pair.data.hitset = nullptr;
      thread_pair.data.rawhitset = hitset;
      thread_pair.data.layer = layer;
      thread_pair.data.pedestal = pedestal;
      thread_pair.data.sector = sector;
      thread_pair.data.side = side;
      thread_pair.data.do_assoc = do_hit_assoc;
      thread_pair.data.do_wedge_emulation = do_wedge_emulation;
      thread_pair.data.tGeometry = m_tGeometry;
      thread_pair.data.maxHalfSizeT = MaxClusterHalfSizeT;
      thread_pair.data.maxHalfSizePhi = MaxClusterHalfSizePhi;
      thread_pair.data.sampa_tbias = m_sampa_tbias;
      thread_pair.data.verbosity = Verbosity();

      unsigned short NPhiBins = (unsigned short) layergeom->get_phibins();
      unsigned short NPhiBinsSector = NPhiBins / 12;
      unsigned short NTBins = (unsigned short) layergeom->get_zbins();
      unsigned short NTBinsSide = NTBins;
      unsigned short NTBinsMin = 0;
      unsigned short PhiOffset = NPhiBinsSector * sector;
      unsigned short TOffset = NTBinsMin;

      m_tdriftmax = AdcClockPeriod * NZBinsSide;
      thread_pair.data.m_tdriftmax = m_tdriftmax;

      thread_pair.data.phibins = NPhiBinsSector;
      thread_pair.data.phioffset = PhiOffset;
      thread_pair.data.tbins = NTBinsSide;
      thread_pair.data.toffset = TOffset;

      /*
      PHG4TpcCylinderGeom *testlayergeom = geom_container->GetLayerCellGeom(32);
      for( float iphi = 1408; iphi < 1408+ 128;iphi+=0.1){
        double clusiphi = iphi;
        double clusphi = testlayergeom->get_phi(clusiphi);
        double radius = layergeom->get_radius();
        float clusx = radius * cos(clusphi);
        float clusy = radius * sin(clusphi);
        float clusz  = -37.524;

        TrkrDefs::hitsetkey tpcHitSetKey = TpcDefs::genHitSetKey( 32,11, 0 );
        Acts::Vector3 global(clusx, clusy, clusz);
        TrkrDefs::subsurfkey subsurfkey = 0;

        Surface surface = m_tGeometry->get_tpc_surface_from_coords(
                                                                   tpcHitSetKey,
                                                                   global,
                                                                   subsurfkey);
        std::cout << " iphi: " << iphi << " clusphi: " << clusphi << " surfkey " << subsurfkey << std::endl;
        //	std::cout << "surfkey" << subsurfkey << std::endl;
      }
      continue;
      */
      int rc = 0;
      //      if(layer==32)
      rc = pthread_create(&thread_pair.thread, &attr, ProcessSector, (void *) &thread_pair.data);
      //      else
      // continue;

      if (rc)
      {
        std::cout << "Error:unable to create thread," << rc << std::endl;
      }

      if (do_sequential)
      {
        int rc2 = pthread_join(thread_pair.thread, nullptr);
        if (rc2)
        {
          std::cout << "Error:unable to join," << rc2 << std::endl;
        }

        // get the hitsetkey from thread data
        const auto &data(thread_pair.data);
        const auto hitsetkey = TpcDefs::genHitSetKey(data.layer, data.sector, data.side);

        // copy clusters to map
        for (uint32_t index = 0; index < data.cluster_vector.size(); ++index)
        {
          // generate cluster key
          const auto ckey = TrkrDefs::genClusKey(hitsetkey, index);

          // get cluster
          auto cluster = data.cluster_vector[index];

          // insert in map
          m_clusterlist->addClusterSpecifyKey(ckey, cluster);
        }

        // copy hit associations to map
        for (const auto &[index, hkey] : thread_pair.data.association_vector)
        {
          // generate cluster key
          const auto ckey = TrkrDefs::genClusKey(hitsetkey, index);

          // add to association table
          m_clusterhitassoc->addAssoc(ckey, hkey);
        }
      }
//      count++;
    }
  }

  pthread_attr_destroy(&attr);
//  count = 0;
  // wait for completion of all threads
  if (!do_sequential)
  {
    for (const auto &thread_pair : threads)
    {
      int rc2 = pthread_join(thread_pair.thread, nullptr);
      if (rc2)
      {
        std::cout << "Error:unable to join," << rc2 << std::endl;
      }

      // get the hitsetkey from thread data
      const auto &data(thread_pair.data);
      const auto hitsetkey = TpcDefs::genHitSetKey(data.layer, data.sector, data.side);

      // copy clusters to map
      for (uint32_t index = 0; index < data.cluster_vector.size(); ++index)
      {
        // generate cluster key
        const auto ckey = TrkrDefs::genClusKey(hitsetkey, index);

        // get cluster
        auto cluster = data.cluster_vector[index];

        // insert in map
        // std::cout << "X: " << cluster->getLocalX() << "Y: " << cluster->getLocalY() << std::endl;
        m_clusterlist->addClusterSpecifyKey(ckey, cluster);

        if (mClusHitsVerbose)
        {
          for (auto &hit : data.phivec_ClusHitsVerbose[index])
          {
            mClusHitsVerbose->addPhiHit(hit.first, (float) hit.second);
          }
          for (auto &hit : data.zvec_ClusHitsVerbose[index])
          {
            mClusHitsVerbose->addZHit(hit.first, (float) hit.second);
          }
          mClusHitsVerbose->push_hits(ckey);
        }
      }

      // copy hit associations to map
      for (const auto &[index, hkey] : thread_pair.data.association_vector)
      {
        // generate cluster key
        const auto ckey = TrkrDefs::genClusKey(hitsetkey, index);

        // add to association table
        m_clusterhitassoc->addAssoc(ckey, hkey);
      }

      for (auto v_hit : thread_pair.data.v_hits)
      {
        if (_store_hits)
        {
          m_training->v_hits.emplace_back(*v_hit);
        }
        delete v_hit;
      }
    }
  }

  // set the flag to use alignment transformations, needed by the rest of reconstruction
  alignmentTransformationContainer::use_alignment = true;

  if (Verbosity() > 0)
  {
    std::cout << "TPC Clusterizer found " << m_clusterlist->size() << " Clusters " << std::endl;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

int TpcClusterizer::End(PHCompositeNode * /*topNode*/)
{
  return Fun4AllReturnCodes::EVENT_OK;
}
