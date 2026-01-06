#include "PHG4TpcPadPlaneReadout.h"

#include <fun4all/Fun4AllReturnCodes.h>
#include <g4detectors/PHG4CellDefs.h>  // for genkey, keytype
#include <g4detectors/PHG4TpcGeom.h>
#include <g4detectors/PHG4TpcGeomContainer.h>

#include <g4main/PHG4Hit.h>  // for PHG4Hit
#include <g4main/PHG4HitContainer.h>

#include <phool/PHRandomSeed.h>
#include <phool/getClass.h>

#include <cdbobjects/CDBTTree.h>
#include <ffamodules/CDBInterface.h>

// Move to new storage containers
#include <trackbase/TpcDefs.h>
#include <trackbase/TrkrDefs.h>  // for hitkey, hitse...
#include <trackbase/TrkrHit.h>   // for TrkrHit
#include <trackbase/TrkrHitSet.h>
#include <trackbase/TrkrHitSetContainer.h>
#include <trackbase/TrkrHitv2.h>  // for TrkrHit

#include <g4tracking/TrkrTruthTrack.h>
#include <g4tracking/TrkrTruthTrackContainer.h>
#include <g4tracking/TrkrTruthTrackContainerv1.h>
#include <g4tracking/TrkrTruthTrackv1.h>

#include <phool/phool.h>  // for PHWHERE

#include <TF1.h>
#include <TFile.h>
#include <TH2.h>
#include <TSystem.h>

#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>  // for gsl_rng_alloc

#include <cmath>
#include <cstdlib>  // for getenv
#include <format>
#include <iostream>
#include <map>      // for _Rb_tree_cons...
#include <utility>  // for pair

class PHCompositeNode;
class TrkrHitTruthAssoc;

namespace
{
  //! convenient square function
  template <class T>
  constexpr T square(const T &x)
  {
    return x * x;
  }

  template <class T>
  inline T get_r(const T &x, const T &y)
  {
    return std::sqrt(square(x) + square(y));
  }

  //! return normalized gaussian centered on zero and of width sigma
  template <class T>
  inline T gaus(const T &x, const T &sigma)
  {
    return std::exp(-square(x / sigma) / 2) / (sigma * std::sqrt(2 * M_PI));
  }

  constexpr unsigned int print_layer = 18;

}  // namespace

PHG4TpcPadPlaneReadout::PHG4TpcPadPlaneReadout(const std::string &name)
  : PHG4TpcPadPlane(name)
  , RandomGenerator(gsl_rng_alloc(gsl_rng_mt19937))
{
  InitializeParameters();
  // if(m_flagToUseGain==1)
  ReadGain();

  gsl_rng_set(RandomGenerator, PHRandomSeed());  // fixed seed is handled in this funtcion

  return;
}

PHG4TpcPadPlaneReadout::~PHG4TpcPadPlaneReadout()
{
  gsl_rng_free(RandomGenerator);
  for (auto *his : h_gain)
  {
    delete his;
  }
}

//_________________________________________________________
int PHG4TpcPadPlaneReadout::InitRun(PHCompositeNode *topNode)
{
  // base class
  const auto reply = PHG4TpcPadPlane::InitRun(topNode);
  if (reply != Fun4AllReturnCodes::EVENT_OK)
  {
    return reply;
  }
  const std::string seggeonodename = "TPCGEOMCONTAINER";
  GeomContainer = findNode::getClass<PHG4TpcGeomContainer>(topNode, seggeonodename);
  assert(GeomContainer);
  
  PHG4TpcGeom *layergeom =  GeomContainer->GetLayerCellGeom(20);  // z geometry is the same for all layers
  double tpc_adc_clock = layergeom->get_adc_clock();
  double extended_readout_time = layergeom->get_extended_readout_time();
  double maxdriftlength = layergeom->get_max_driftlength();
  double drift_velocity_sim = layergeom->get_drift_velocity_sim();
  const double TBinWidth = tpc_adc_clock;
  const double MaxT = extended_readout_time + 2.0 * maxdriftlength / drift_velocity_sim;  // allows for extended time readout
  const double MinT = 0;
  NTBins = (int) ((MaxT - MinT) / TBinWidth) + 1;

  if (m_use_module_gain_weights)
  {
    int side;
    int region;
    int sector;
    double weight;
    std::ifstream weights_file(m_tpc_module_gain_weights_file);
    if (!weights_file.is_open())
    {
      std::cout << ".In PHG4TpcPadPlaneReadout: Option to use module gain weights enabled, but weights file not found. Aborting." << std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }

    for (int iside = 0; iside < 2; ++iside)
    {
      for (int isec = 0; isec < 12; ++isec)
      {
        for (int ir = 0; ir < 3; ++ir)
        {
          weights_file >> side >> region >> sector >> weight;
          m_module_gain_weight[side][region][sector] = weight;
          std::cout << " iside " << iside << " side " << side << " ir " << ir
                    << " region " << region << " isec " << isec
                    << " sector " << sector << " weight " << weight << std::endl;
        }
      }
    }
  }

  if (m_useLangau)
  {
    int side;
    int region;
    int sector;
    double par0;
    double par1;
    double par2;
    double par3;
    std::ifstream pars_file(m_tpc_langau_pars_file);
    if (!pars_file.is_open())
    {
      std::cout << ".In PHG4TpcPadPlaneReadout: Option to use Langau parameters enabled, but parameter file not found. Aborting." << std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }

    for (int iside = 0; iside < 2; ++iside)
    {
      for (int isec = 0; isec < 12; ++isec)
      {
        for (int ir = 0; ir < 3; ++ir)
        {
          pars_file >> side >> region >> sector >> par0 >> par1 >> par2 >> par3;
          flangau[side][region][sector] = new TF1(std::format("flangau_{}_{}_{}", side, region, sector).c_str(), [](double *x, double *par)
                                                  {
		    Double_t invsq2pi = 0.3989422804014;
		    Double_t mpshift  = -0.22278298;
		    Double_t np = 100.0;
		    Double_t sc =   5.0;
		    Double_t xx;
		    Double_t mpc;
		    Double_t fland;
		    Double_t sum = 0.0;
		    Double_t xlow;
		    Double_t xupp;
		    Double_t step;
		    Double_t i;
		    mpc = par[1] - mpshift * par[0]; 
		    xlow = x[0] - sc * par[3];
		    xupp = x[0] + sc * par[3];
		    step = (xupp-xlow) / np;
		    for(i=1.0; i<=np/2; i++) 
		    {
		      xx = xlow + (i-.5) * step;
		      fland = TMath::Landau(xx,mpc,par[0]) / par[0];
		      sum += fland * TMath::Gaus(x[0],xx,par[3]);
		      
		      xx = xupp - (i-.5) * step;
		      fland = TMath::Landau(xx,mpc,par[0]) / par[0];
		      sum += fland * TMath::Gaus(x[0],xx,par[3]);
		    }
      
		    return (par[2] * step * sum * invsq2pi / par[3]); }, 0, 5000, 4);

          flangau[side][region][sector]->SetParameters(par0, par1, par2, par3);
          // std::cout << " iside " << iside << " side " << side << " ir " << ir
          //	    << " region " << region << " isec " << isec
          //	    << " sector " << sector << " weight " << weight << std::endl;
        }
      }
    }
  }
  if (m_maskDeadChannels)
  {
    makeChannelMask(m_deadChannelMap, m_deadChannelMapName, "TotalDeadChannels");
  }
  if (m_maskHotChannels)
  {
    makeChannelMask(m_hotChannelMap, m_hotChannelMapName, "TotalHotChannels");
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

//_________________________________________________________
double PHG4TpcPadPlaneReadout::getSingleEGEMAmplification()
{
  // Jin H.: For the GEM gain in sPHENIX TPC,
  //         Bob pointed out the PHENIX HBD measured it as the Polya function with theta parameter = 0.8.
  //         Just talked with Tom too, he suggest us to start the TPC modeling with simpler exponential function
  //         with lambda parameter of 1/2000, (i.e. Polya function with theta parameter = 0, q_bar = 2000). Please note, this gain variation need to be sampled for each initial electron individually.
  //         Summing over ~30 initial electrons, the distribution is pushed towards more Gauss like.
  // Bob A.: I like Tom's suggestion to use the exponential distribution as a first approximation
  //         for the single electron gain distribution -
  //         and yes, the parameter you're looking for is of course the slope, which is the inverse gain.
  double nelec = gsl_ran_exponential(RandomGenerator, averageGEMGain);
  if (m_usePolya)
  {
    double y;
    double xmax = 5000;
    double ymax = 0.376;
    while (true)
    {
      nelec = gsl_ran_flat(RandomGenerator, 0, xmax);
      y = gsl_rng_uniform(RandomGenerator) * ymax;
      if (y <= pow((1 + polyaTheta) * (nelec / averageGEMGain), polyaTheta) * exp(-(1 + polyaTheta) * (nelec / averageGEMGain)))
      {
        break;
      }
    }
  }
  // Put gain reading here

  return nelec;
}

//_________________________________________________________
double PHG4TpcPadPlaneReadout::getSingleEGEMAmplification(double weight)
{
  // Jin H.: For the GEM gain in sPHENIX TPC,
  //         Bob pointed out the PHENIX HBD measured it as the Polya function with theta parameter = 0.8.
  //         Just talked with Tom too, he suggest us to start the TPC modeling with simpler exponential function
  //         with lambda parameter of 1/2000, (i.e. Polya function with theta parameter = 0, q_bar = 2000). Please note, this gain variation need to be sampled for each initial electron individually.
  //         Summing over ~30 initial electrons, the distribution is pushed towards more Gauss like.
  // Bob A.: I like Tom's suggestion to use the exponential distribution as a first approximation
  //         for the single electron gain distribution -
  //         and yes, the parameter you're looking for is of course the slope, which is the inverse gain.
  double q_bar = averageGEMGain * weight;
  double nelec = gsl_ran_exponential(RandomGenerator, q_bar);
  if (m_usePolya)
  {
    double y;
    double xmax = 5000;
    double ymax = 0.376;
    while (true)
    {
      nelec = gsl_ran_flat(RandomGenerator, 0, xmax);
      y = gsl_rng_uniform(RandomGenerator) * ymax;
      if (y <= pow((1 + polyaTheta) * (nelec / q_bar), polyaTheta) * exp(-(1 + polyaTheta) * (nelec / q_bar)))
      {
        break;
      }
    }
  }
  // Put gain reading here

  return nelec;
}

//_________________________________________________________
double PHG4TpcPadPlaneReadout::getSingleEGEMAmplification(TF1 *f)
{
  double nelec = f->GetRandom(0, 5000);
  // Put gain reading here

  return nelec;
}

void PHG4TpcPadPlaneReadout::MapToPadPlane(
    TpcClusterBuilder &tpc_truth_clusterer,
    TrkrHitSetContainer *single_hitsetcontainer,
    TrkrHitSetContainer *hitsetcontainer,
    TrkrHitTruthAssoc * /*hittruthassoc*/,
    const double x_gem, const double y_gem, const double t_gem, const unsigned int side,
    PHG4HitContainer::ConstIterator hiter, TNtuple * /*ntpad*/, TNtuple * /*nthit*/)
{
  // One electron per call of this method
  // The x_gem and y_gem values have already been randomized within the transverse drift diffusion width
  // The t_gem value already reflects the drift time of the primary electron from the production point, and is randomized within the longitudinal diffusion witdth

  double phi = atan2(y_gem, x_gem);
  if (phi > +M_PI)
  {
    phi -= 2 * M_PI;
  }
  if (phi < -M_PI)
  {
    phi += 2 * M_PI;
  }

  double rad_gem = get_r(x_gem, y_gem);

  // Moving electrons from dead area to a closest pad
  for (int iregion = 0; iregion < 3; ++iregion)
  {
    double daR = 0;
    if (iregion == 0 || iregion == 2)
    {
      daR = 1.0;  // 1.0cm edge to collect electrons from
    }
    else
    {
      daR = MinRadius[iregion] - MaxRadius[iregion - 1];
    }
    if (rad_gem <= MinRadius[iregion] && rad_gem >= MinRadius[iregion] - daR)
    {
      if (rad_gem <= MinRadius[iregion] - daR / 2)
      {
        rad_gem = MinRadius[iregion] - (1.1 * daR);
      }
      else
      {
        rad_gem = MinRadius[iregion] + 0.1 * daR;
      }
    }
  }

  unsigned int layernum = 0;
  /* TpcClusterBuilder pass_data {}; */

  // Find which readout layer this electron ends up in

  PHG4TpcGeomContainer::ConstRange layerrange = GeomContainer->get_begin_end();
  for (PHG4TpcGeomContainer::ConstIterator layeriter = layerrange.first;
       layeriter != layerrange.second;
       ++layeriter)
  {
    double rad_low = layeriter->second->get_radius() - layeriter->second->get_thickness() / 2.0;
    double rad_high = layeriter->second->get_radius() + layeriter->second->get_thickness() / 2.0;

    if (rad_gem > rad_low && rad_gem < rad_high)
    {
      // capture the layer where this electron hits the gem stack
      LayerGeom = layeriter->second;

      layernum = LayerGeom->get_layer();
      /* pass_data.layerGeom = LayerGeom; */
      /* pass_data.layer = layernum; */
      if (Verbosity() > 1000)
      {
        std::cout << " g4hit id " << hiter->first << " rad_gem " << rad_gem << " rad_low " << rad_low << " rad_high " << rad_high
                  << " layer  " << hiter->second->get_layer() << " want to change to " << layernum << std::endl;
      }
      hiter->second->set_layer(layernum);  // have to set here, since the stepping action knows nothing about layers
    }
  }

  if (layernum == 0)
  {
    return;
  }

  // store phi bins and tbins upfront to avoid repetitive checks on the phi methods
  const auto phibins = LayerGeom->get_phibins();
  /* pass_data.nphibins = phibins; */

  const auto tbins = LayerGeom->get_zbins();

  sector_min_Phi = LayerGeom->get_sector_min_phi();
  sector_max_Phi = LayerGeom->get_sector_max_phi();
  phi_bin_width = LayerGeom->get_phistep();

  phi = check_phi(side, phi, rad_gem);

  // Create the distribution function of charge on the pad plane around the electron position

  // The resolution due to pad readout includes the charge spread during GEM multiplication.
  // this now defaults to 400 microns during construction from Tom (see 8/11 email).
  // Use the setSigmaT(const double) method to update...
  // We use a double gaussian to represent the smearing due to the SAMPA chip shaping time - default values of fShapingLead and fShapingTail are for 80 ns SAMPA

  // amplify the single electron in the gem stack
  //===============================

  double nelec = getSingleEGEMAmplification();
  // Applying weight with respect to the rad_gem and phi after electrons are redistributed
  double phi_gain = phi;
  if (phi < 0)
  {
    phi_gain += 2 * M_PI;
  }
  double gain_weight = 1.0;
  if (m_flagToUseGain == 1)
  {
    gain_weight = h_gain[side]->GetBinContent(h_gain[side]->FindBin(rad_gem * 10, phi_gain));  // rad_gem in cm -> *10 to get mm
    nelec = nelec * gain_weight;
  }

  if (m_use_module_gain_weights)
  {
    double phistep = 30.0;
    int sector = 0;

    if ((phi_gain * 180.0 / M_PI) >= 15 && (phi_gain * 180.0 / M_PI) < 345)
    {
      sector = 1 + (int) ((phi_gain * 180.0 / M_PI - 15) / phistep);
    }
    else
    {
      sector = 0;
    }

    int this_region = -1;
    for (int iregion = 0; iregion < 3; ++iregion)
    {
      if (rad_gem < MaxRadius[iregion] && rad_gem > MinRadius[iregion])
      {
        this_region = iregion;
      }
    }
    if (this_region > -1)
    {
      gain_weight = m_module_gain_weight[side][this_region][sector];
    }
    // regenerate nelec with the new distribution
    //    double original_nelec = nelec;
    nelec = getSingleEGEMAmplification(gain_weight);
    //  std::cout << " side " << side << " this_region " << this_region
    //	<<  " sector " << sector << " original nelec "
    //	<< original_nelec << " new nelec " << nelec << std::endl;
  }

  if (m_useLangau)
  {
    double phistep = 30.0;
    int sector = 0;

    if ((phi_gain * 180.0 / M_PI) >= 15 && (phi_gain * 180.0 / M_PI) < 345)
    {
      sector = 1 + (int) ((phi_gain * 180.0 / M_PI - 15) / phistep);
    }
    else
    {
      sector = 0;
    }

    int this_region = -1;
    for (int iregion = 0; iregion < 3; ++iregion)
    {
      if (rad_gem < MaxRadius[iregion] && rad_gem > MinRadius[iregion])
      {
        this_region = iregion;
      }
    }
    if (this_region > -1)
    {
      nelec = getSingleEGEMAmplification(flangau[side][this_region][sector]);
    }
    else
    {
      nelec = getSingleEGEMAmplification();
    }
  }

  // std::cout<<"PHG4TpcPadPlaneReadout::MapToPadPlane gain_weight = "<<gain_weight<<std::endl;
  /* pass_data.neff_electrons = nelec; */

  // Distribute the charge between the pads in phi
  //====================================

  if (Verbosity() > 200)
  {
    std::cout << "  populate phi bins for "
              << " layernum " << layernum
              << " phi " << phi
              << " sigmaT " << sigmaT
              //<< " zigzag_pads " << zigzag_pads
              << std::endl;
  }

  std::vector<int> pad_phibin;
  std::vector<double> pad_phibin_share;

  populate_zigzag_phibins(side, layernum, phi, sigmaT, pad_phibin, pad_phibin_share);
  /* if (pad_phibin.size() == 0) { */
  /* pass_data.neff_electrons = 0; */
  /* } else { */
  /* pass_data.fillPhiBins(pad_phibin); */
  /* } */

  // Normalize the shares so they add up to 1
  double norm1 = 0.0;
  for (unsigned int ipad = 0; ipad < pad_phibin.size(); ++ipad)
  {
    double pad_share = pad_phibin_share[ipad];
    norm1 += pad_share;
  }
  for (unsigned int iphi = 0; iphi < pad_phibin.size(); ++iphi)
  {
    pad_phibin_share[iphi] /= norm1;
  }

  // Distribute the charge between the pads in t
  //====================================
  if (Verbosity() > 100 && layernum == print_layer)
    {
      std::cout << "  populate t bins for layernum " << layernum
		<< " with t_gem " << t_gem << " SAMPA peaking time  " << Ts << std::endl;
    }

  std::vector<int> adc_tbin;
  std::vector<double> adc_tbin_share;
  sampaTimeDistribution(t_gem, adc_tbin, adc_tbin_share);

  /* if (adc_tbin.size() == 0)  { */
  /* pass_data.neff_electrons = 0; */
  /* } else { */
  /* pass_data.fillTimeBins(adc_tbin); */
  /* } */

  // Normalize the shares so that they add up to 1
  double tnorm = 0.0;
  for (unsigned int it = 0; it < adc_tbin.size(); ++it)
  {
    double bin_share = adc_tbin_share[it];
    tnorm += bin_share;
  }
  for (unsigned int it = 0; it < adc_tbin.size(); ++it)
  {
    adc_tbin_share[it] /= tnorm;
  }

  /*
  if(layernum == print_layer)
    {
      std::cout << "t_gem " << t_gem << std::endl;
      for (unsigned int it = 0; it < adc_tbin.size(); ++it)
	{
	  std::cout << " tbin " << adc_tbin[it] << " share " << adc_tbin_share[it] << std::endl;
	}
    }
  */
  
  // Fill HitSetContainer
  //===============
  // These are used to do a quick clustering for checking
  double phi_integral = 0.0;
  double t_integral = 0.0;
  double weight = 0.0;

  for (unsigned int ipad = 0; ipad < pad_phibin.size(); ++ipad)
  {
    int pad_num = pad_phibin[ipad];
    double pad_share = pad_phibin_share[ipad];

    for (unsigned int it = 0; it < adc_tbin.size(); ++it)
    {
      int tbin_num = adc_tbin[it];
      double adc_bin_share = adc_tbin_share[it];

      // Divide electrons from avalanche between bins
      float neffelectrons = nelec * (pad_share) * (adc_bin_share);
      if (neffelectrons < neffelectrons_threshold)
      {
        continue;  // skip signals that will be below the noise suppression threshold
      }

      if (tbin_num >= tbins)
      {
        std::cout << " Error making key: adc_tbin " << tbin_num << " ntbins " << tbins << std::endl;
      }
      if (pad_num >= phibins)
      {
        std::cout << " Error making key: pad_phibin " << pad_num << " nphibins " << phibins << std::endl;
      }

      // collect information to do simple clustering. Checks operation of PHG4CylinderCellTpcReco, and
      // is also useful for comparison with PHG4TpcClusterizer result when running single track events.
      // The only information written to the cell other than neffelectrons is tbin and pad number, so get those from geometry
      double tcenter = LayerGeom->get_zcenter(tbin_num);
      double phicenter = LayerGeom->get_phicenter(pad_num, side);
      phi_integral += phicenter * neffelectrons;
      t_integral += tcenter * neffelectrons;
      weight += neffelectrons;
      if (Verbosity() > 1 && layernum == print_layer)
      {
        std::cout << "   tbin_num " << tbin_num << " tcenter " << tcenter << " pad_num " << pad_num << " phicenter " << phicenter
                  << " neffelectrons " << neffelectrons << " neffelectrons_threshold " << neffelectrons_threshold << std::endl;
      }

      // new containers
      //============
      // We add the Tpc TrkrHitsets directly to the node using hitsetcontainer
      // We need to create the TrkrHitSet if not already made - each TrkrHitSet should correspond to a Tpc readout module
      // The hitset key includes the layer, sector, side

      // The side is an input parameter

      // get the Tpc readout sector - there are 12 sectors with how many pads each?
      unsigned int pads_per_sector = phibins / 12;
      unsigned int sector = pad_num / pads_per_sector;
      TrkrDefs::hitsetkey hitsetkey = TpcDefs::genHitSetKey(layernum, sector, side);
      // Use existing hitset or add new one if needed
      TrkrHitSetContainer::Iterator hitsetit = hitsetcontainer->findOrAddHitSet(hitsetkey);
      TrkrHitSetContainer::Iterator single_hitsetit = single_hitsetcontainer->findOrAddHitSet(hitsetkey);
      TrkrDefs::hitkey hitkey;

      if (m_maskDeadChannels)
      {
        hitkey = TpcDefs::genHitKey((unsigned int) pad_num, 0);
        if (m_deadChannelMap.contains(hitsetkey) &&
            std::find(m_deadChannelMap[hitsetkey].begin(), m_deadChannelMap[hitsetkey].end(), hitkey) != m_deadChannelMap[hitsetkey].end())
        {
          continue;
        }
      }
      if (m_maskHotChannels)
      {
        hitkey = TpcDefs::genHitKey((unsigned int) pad_num, 0);
        if (m_hotChannelMap.contains(hitsetkey) &&
            std::find(m_hotChannelMap[hitsetkey].begin(), m_hotChannelMap[hitsetkey].end(), hitkey) != m_hotChannelMap[hitsetkey].end())
        {
          continue;
        }
      }
      // generate the key for this hit, requires tbin and phibin
      hitkey = TpcDefs::genHitKey((unsigned int) pad_num, (unsigned int) tbin_num);

      // See if this hit already exists
      TrkrHit *hit = nullptr;
      hit = hitsetit->second->getHit(hitkey);
      if (!hit)
      {
        // create a new one
        hit = new TrkrHitv2();
        hitsetit->second->addHitSpecificKey(hitkey, hit);
      }
      // Either way, add the energy to it  -- adc values will be added at digitization
      hit->addEnergy(neffelectrons);

      tpc_truth_clusterer.addhitset(hitsetkey, hitkey, neffelectrons);

      // repeat for the single_hitsetcontainer
      // See if this hit already exists
      TrkrHit *single_hit = nullptr;
      single_hit = single_hitsetit->second->getHit(hitkey);
      if (!single_hit)
      {
        // create a new one
        single_hit = new TrkrHitv2();
        single_hitsetit->second->addHitSpecificKey(hitkey, single_hit);
      }
      // Either way, add the energy to it  -- adc values will be added at digitization
      single_hit->addEnergy(neffelectrons);

      /*
      if (Verbosity() > 0)
      {
        assert(nthit);
        nthit->Fill(layernum, pad_num, tbin_num, neffelectrons);
      }
      */

    }  // end of loop over adc T bins
  }  // end of loop over zigzag pads
  /* pass_data.phi_integral = phi_integral; */
  /* pass_data.time_integral = t_integral; */

  /*
  // Capture the input values at the gem stack and the quick clustering results, elecron-by-electron
  if (Verbosity() > 0)
  {
    assert(ntpad);
    ntpad->Fill(layernum, phi, phi_integral / weight, t_gem, t_integral / weight);
  }
  */

  if (Verbosity() > 100)
  {
    if (layernum == print_layer)
    {
      std::cout << " hit " << m_NHits << " quick centroid for this electron " << std::endl;
      std::cout << "      phi centroid = " << phi_integral / weight << " phi in " << phi << " phi diff " << phi_integral / weight - phi << std::endl;
      std::cout << "      t centroid = " << t_integral / weight << " t in " << t_gem << " t diff " << t_integral / weight - t_gem << std::endl;
      // For a single track event, this captures the distribution of single electron centroids on the pad plane for layer print_layer.
      // The centroid of that should match the cluster centroid found by PHG4TpcClusterizer for layer print_layer, if everything is working
      //   - matches to < .01 cm for a few cases that I checked

      /*
      assert(nthit);
      nthit->Fill(hit, layernum, phi, phi_integral / weight, t_gem, t_integral / weight, weight);
      */
    }
  }

  m_NHits++;
  /* return pass_data; */
}
double PHG4TpcPadPlaneReadout::check_phi(const unsigned int side, const double phi, const double radius)
{
  double new_phi = phi;
  int p_region = -1;
  for (int iregion = 0; iregion < 3; ++iregion)
  {
    if (radius < MaxRadius[iregion] && radius > MinRadius[iregion])
    {
      p_region = iregion;
    }
  }

  if (p_region >= 0)
  {
    for (int s = 0; s < 12; s++)
    {
      double daPhi = 0;
      if (s == 0)
      {
        daPhi = fabs(sector_min_Phi[side][11] + 2 * M_PI - sector_max_Phi[side][s]);
      }
      else
      {
        daPhi = fabs(sector_min_Phi[side][s - 1] - sector_max_Phi[side][s]);
      }
      double min_phi = sector_max_Phi[side][s];
      double max_phi = sector_max_Phi[side][s] + daPhi;
      if (new_phi <= max_phi && new_phi >= min_phi)
      {
        if (fabs(max_phi - new_phi) > fabs(new_phi - min_phi))
        {
          new_phi = min_phi - phi_bin_width / 5;
        }
        else
        {
          new_phi = max_phi + phi_bin_width / 5;
        }
      }
    }
    if (new_phi < sector_min_Phi[side][11] && new_phi >= -M_PI)
    {
      new_phi += 2 * M_PI;
    }
  }

  return new_phi;
}

void PHG4TpcPadPlaneReadout::populate_zigzag_phibins(const unsigned int side, const unsigned int layernum, const double phi, const double cloud_sig_rp, std::vector<int> &phibin_pad, std::vector<double> &phibin_pad_share)
{
  const double radius = LayerGeom->get_radius();
  const double phistepsize = LayerGeom->get_phistep();
  const auto phibins = LayerGeom->get_phibins();

  // make the charge distribution gaussian
  double rphi = phi * radius;
  if (Verbosity() > 100)
  {
    if (LayerGeom->get_layer() == print_layer)
    {
      std::cout << " populate_zigzag_phibins for layer " << layernum << " with radius " << radius << " phi " << phi
                << " rphi " << rphi << " phistepsize " << phistepsize << std::endl;
      std::cout << " fcharge created: radius " << radius << " rphi " << rphi << " cloud_sig_rp " << cloud_sig_rp << std::endl;
    }
  }

  // Get the range of phi values that completely contains all pads  that touch the charge distribution - (nsigmas + 1/2 pad width) in each direction
  const double philim_low_calc = phi - (_nsigmas * cloud_sig_rp / radius) - phistepsize;
  const double philim_high_calc = phi + (_nsigmas * cloud_sig_rp / radius) + phistepsize;

  // Find the pad range that covers this phi range
  const double philim_low = check_phi(side, philim_low_calc, radius);
  const double philim_high = check_phi(side, philim_high_calc, radius);

  int phibin_low = LayerGeom->get_phibin(philim_high, side);
  int phibin_high = LayerGeom->get_phibin(philim_low, side);
  int npads = phibin_high - phibin_low;

  if (Verbosity() > 1000)
  {
    if (layernum == print_layer)
    {
      std::cout << "           zigzags: phi " << phi << " philim_low " << philim_low << " phibin_low " << phibin_low
                << " philim_high " << philim_high << " phibin_high " << phibin_high << " npads " << npads << std::endl;
    }
  }

  if (npads < 0 || npads > 9)
  {
    npads = 9;  // can happen if phibin_high wraps around. If so, limit to 10 pads and fix below
  }

  // Calculate the maximum extent in r-phi of pads in this layer. Pads are assumed to touch the center of the next phi bin on both sides.
  const double pad_rphi = 2.0 * LayerGeom->get_phistep() * radius;

  // Make a TF1 for each pad in the phi range
  using PadParameterSet = std::array<double, 2>;
  std::array<PadParameterSet, 10> pad_parameters{};
  std::array<int, 10> pad_keep{};

  // Now make a loop that steps through the charge distribution and evaluates the response at that point on each pad
  std::array<double, 10> overlap = {{0, 0, 0, 0, 0, 0, 0, 0, 0, 0}};
  double pads_phi[10] = {0., 0., 0., 0., 0., 0., 0., 0., 0., 0.};
  double sum_of_pads_phi = 0.;
  double sum_of_pads_absphi = 0.;
  for (int ipad = 0; ipad <= npads; ipad++)
  {
    int pad_now = phibin_low + ipad;
    if (pad_now >= phibins)
    {
      pad_now -= phibins;
    }
    pads_phi[ipad] = LayerGeom->get_phicenter(pad_now, side);
    sum_of_pads_phi += pads_phi[ipad];
    sum_of_pads_absphi += fabs(pads_phi[ipad]);
  }

  for (int ipad = 0; ipad <= npads; ipad++)
  {
    int pad_now = phibin_low + ipad;
    // if(phibin_low<0 && phibin_high<0) pad_now = phibin_high + ipad;
    //  check that we do not exceed the maximum number of pads, wrap if necessary
    if (pad_now >= phibins)
    {
      pad_now -= phibins;
    }

    pad_keep[ipad] = pad_now;
    double phi_now = pads_phi[ipad];
    const double rphi_pad_now = phi_now * radius;
    pad_parameters[ipad] = {{pad_rphi / 2.0, rphi_pad_now}};

    if (Verbosity() > 1000)
    {
      if (layernum == print_layer)
      {
        std::cout << " zigzags: make fpad for ipad " << ipad << " pad_now " << pad_now << " pad_rphi/2 " << pad_rphi / 2.0
                  << " rphi_pad_now " << rphi_pad_now << std::endl;
      }
    }
    //}

    // use analytic integral
    // for (int ipad = 0; ipad <= npads; ipad++)
    //{
    // const double pitch = pad_parameters[ipad][0];
    // const double x_loc = pad_parameters[ipad][1] - rphi;
    // const double sigma = cloud_sig_rp;

    const double pitch = pad_rphi / 2.0;     // eshulga
    double x_loc_tmp = rphi_pad_now - rphi;  // eshulga
    const double sigma = cloud_sig_rp;       // eshulga

    // Checking if the pads are on the same side of the TPC in phi
    if (fabs(sum_of_pads_phi) != sum_of_pads_absphi)
    {
      if (phi < -M_PI / 2 && phi_now > 0)
      {
        x_loc_tmp = (phi_now - 2 * M_PI) * radius - rphi;
      }
      if (phi > M_PI / 2 && phi_now < 0)
      {
        x_loc_tmp = (phi_now + 2 * M_PI) * radius - rphi;
      }
      if (phi < 0 && phi_now > 0)
      {
        x_loc_tmp = (phi_now + fabs(phi)) * radius;
      }
      if (phi > 0 && phi_now < 0)
      {
        x_loc_tmp = (2 * M_PI - phi_now + phi) * radius;
      }
    }

    const double x_loc = x_loc_tmp;
    // calculate fraction of the total charge on this strip
    /*
    this corresponds to integrating the charge distribution Gaussian function (centered on rphi and of width cloud_sig_rp),
    convoluted with a strip response function, which is triangular from -pitch to +pitch, with a maximum of 1. at stript center
    */
    overlap[ipad] =
        (pitch - x_loc) * (std::erf(x_loc / (M_SQRT2 * sigma)) - std::erf((x_loc - pitch) / (M_SQRT2 * sigma))) / (pitch * 2) + (pitch + x_loc) * (std::erf((x_loc + pitch) / (M_SQRT2 * sigma)) - std::erf(x_loc / (M_SQRT2 * sigma))) / (pitch * 2) + (gaus(x_loc - pitch, sigma) - gaus(x_loc, sigma)) * square(sigma) / pitch + (gaus(x_loc + pitch, sigma) - gaus(x_loc, sigma)) * square(sigma) / pitch;
  }

  // now we have the overlap for each pad
  for (int ipad = 0; ipad <= npads; ipad++)
  {
    phibin_pad.push_back(pad_keep[ipad]);
    phibin_pad_share.push_back(overlap[ipad]);
  }

  return;
}

void PHG4TpcPadPlaneReadout::UseGain(const int flagToUseGain)
{
  m_flagToUseGain = flagToUseGain;
  if (m_flagToUseGain == 1 && Verbosity() > 0)
  {
    std::cout << "PHG4TpcPadPlaneReadout: UseGain: TRUE " << std::endl;
  }
}

void PHG4TpcPadPlaneReadout::ReadGain()
{
  // Reading TPC Gain Maps from the file
  if (m_flagToUseGain == 1)
  {
    char *calibrationsroot = getenv("CALIBRATIONROOT");
    if (calibrationsroot != nullptr)
    {
      std::string gain_maps_filename = std::string(calibrationsroot) + std::string("/TPC/GainMaps/TPCGainMaps.root");
      TFile *fileGain = TFile::Open(gain_maps_filename.c_str(), "READ");
      h_gain[0] = (TH2 *) fileGain->Get("RadPhiPlot0")->Clone();
      h_gain[1] = (TH2 *) fileGain->Get("RadPhiPlot1")->Clone();
      h_gain[0]->SetDirectory(nullptr);
      h_gain[1]->SetDirectory(nullptr);
      fileGain->Close();
    }
  }
}
void PHG4TpcPadPlaneReadout::SetDefaultParameters()
{
  set_default_double_param("tpc_minradius_inner", 31.105);  // 30.0);  // cm
  set_default_double_param("tpc_minradius_mid", 41.153);    // 40.0);
  set_default_double_param("tpc_minradius_outer", 58.367);  // 60.0);

  set_default_double_param("tpc_maxradius_inner", 40.249);  // 40.0);  // cm
  set_default_double_param("tpc_maxradius_mid", 57.475);    // 60.0);
  set_default_double_param("tpc_maxradius_outer", 75.911);  // 77.0);  // from Tom

  set_default_double_param("neffelectrons_threshold", 1.0);
  set_default_double_param("gem_cloud_sigma", 0.04);     // cm = 400 microns
  set_default_double_param("sampa_shaping_lead", 32.0);  // ns, for 80 ns SAMPA
  set_default_double_param("sampa_shaping_tail", 48.0);  // ns, for 80 ns SAMPA

  set_default_double_param("tpc_sector_phi_inner", 0.5024);  // 2 * M_PI / 12 );//sector size in phi for R1 sector
  set_default_double_param("tpc_sector_phi_mid", 0.5087);    // 2 * M_PI / 12 );//sector size in phi for R2 sector
  set_default_double_param("tpc_sector_phi_outer", 0.5097);  // 2 * M_PI / 12 );//sector size in phi for R3 sector

  set_default_int_param("ntpc_phibins_inner", 1128);  // 94 * 12
  set_default_int_param("ntpc_phibins_mid", 1536);    // 128 * 12
  set_default_int_param("ntpc_phibins_outer", 2304);  // 192 * 12

  // GEM Gain
  /*
  hp (2020/09/04): gain changed from 2000 to 1400, to accomodate gas mixture change
  from Ne/CF4 90/10 to Ne/CF4 50/50, and keep the average charge per particle per pad constant
  */
  set_default_double_param("gem_amplification", 1400);
  set_default_double_param("polya_theta", 0.8);
  return;
}

void PHG4TpcPadPlaneReadout::UpdateInternalParameters()
{
  neffelectrons_threshold = get_double_param("neffelectrons_threshold");

  MinRadius =
      {{get_double_param("tpc_minradius_inner"),
        get_double_param("tpc_minradius_mid"),
        get_double_param("tpc_minradius_outer")}};

  MaxRadius =
      {{get_double_param("tpc_maxradius_inner"),
        get_double_param("tpc_maxradius_mid"),
        get_double_param("tpc_maxradius_outer")}};

  sigmaT = get_double_param("gem_cloud_sigma");
  sigmaL = {{get_double_param("sampa_shaping_lead"),
             get_double_param("sampa_shaping_tail")}};


  averageGEMGain = get_double_param("gem_amplification");
  polyaTheta = get_double_param("polya_theta");
}

void PHG4TpcPadPlaneReadout::makeChannelMask(hitMaskTpc &aMask, const std::string &dbName, const std::string &totalChannelsToMask)
{
  CDBTTree *cdbttree;
  if (m_maskFromFile)
  {
    cdbttree = new CDBTTree(dbName);
  }
  else // mask using CDB TTree, default
  {
    std::string database = CDBInterface::instance()->getUrl(dbName);
    cdbttree = new CDBTTree(database);
  }
  
  std::cout << "Masking TPC Channel Map: " << dbName << std::endl;

  int NChan = -1;
  NChan = cdbttree->GetSingleIntValue(totalChannelsToMask);

  for (int i = 0; i < NChan; i++)
  {
    int Layer = cdbttree->GetIntValue(i, "layer");
    int Sector = cdbttree->GetIntValue(i, "sector");
    int Side = cdbttree->GetIntValue(i, "side");
    int Pad = cdbttree->GetIntValue(i, "pad");
    if (Verbosity() > VERBOSITY_A_LOT)
    {
      std::cout << dbName << ": Will mask layer: " << Layer << ", sector: " << Sector << ", side: " << Side << ", Pad: " << Pad << std::endl;
    }

    TrkrDefs::hitsetkey DeadChannelHitKey = TpcDefs::genHitSetKey(Layer, Sector, Side);
    TrkrDefs::hitkey DeadHitKey = TpcDefs::genHitKey((unsigned int) Pad, 0);
    aMask[DeadChannelHitKey].push_back(DeadHitKey);
  }

  delete cdbttree;
}

void PHG4TpcPadPlaneReadout::sampaTimeDistribution(double tzero,  std::vector<int> &adc_tbin, std::vector<double> &adc_tbin_share)
{
  // tzero is the arrival time of the electron at the GEM
  // Ts is the sampa peaking time
  // Assume the response is over after 8 clock cycles (400 ns)
  int nclocks = 8;

  double tstepsize = LayerGeom->get_zstep();
  int tbinzero = LayerGeom->get_zbin(tzero);

  // the first clock bin is a special case
  double tfirst_end = LayerGeom->get_zcenter(tbinzero) + tstepsize/2.0;
  double vfirst_end =  sampaShapingResponseFunction(tzero, tfirst_end); 
  double first_integral = (vfirst_end / 2.0) * (tfirst_end - tzero);
    
  adc_tbin.push_back(tbinzero);
  adc_tbin_share.push_back(first_integral);

  /*
  if (LayerGeom->get_layer() == print_layer)
    {
      std::cout << "     tzero " << tzero << " tbinzero " << tbinzero << " iclock  0 "  
		<< " tfirst_end " << tfirst_end << " vfirst_end " << vfirst_end << " first_integral " << first_integral << std::endl;      
    }
  */
  
  for(int iclock = 1; iclock < nclocks; ++iclock)
    {
      int tbin = tbinzero + iclock;
      if (tbin < 0 || tbin > LayerGeom->get_zbins())
	{
	  if (Verbosity() > 0)
	    {
	      std::cout << " t bin " << tbin << " is outside range of " << LayerGeom->get_zbins() << " so skip it" << std::endl;
	    }
	  continue;
	}

      // get the beginning and end of this clock bin
      double tcenter = LayerGeom->get_zcenter(tbin);
      double tlow = tcenter - tstepsize/2.0;

      // sample the voltage in this bin at nsamples-1 locations
      int nsamples = 6;
      double sample_step = tstepsize / (double) nsamples;
      double sintegral = 0;
      for(int isample = 0; isample < nsamples; ++isample)
	{
	  double tnow = tlow + (double) isample * sample_step + sample_step / 2.0;	  
	  double vnow = sampaShapingResponseFunction(tzero, tnow);
	  sintegral += vnow * sample_step;

	  /*
	  if (LayerGeom->get_layer() == print_layer)
	    {
	      std::cout << "     tzero " << tzero << " tbinzero " << tbinzero << " iclock " << iclock << " tbin " << tbin << " isample " << isample
			<< " tnow " << tnow << " vnow " << vnow  << " sintegral " << sintegral << std::endl;
	    }
	  */
	}

	  
      adc_tbin.push_back(tbin);
      adc_tbin_share.push_back(sintegral);      
    }
}
  
double PHG4TpcPadPlaneReadout::sampaShapingResponseFunction(double tzero, double t) const
  {
    double v = exp(-4*(t-tzero)/Ts) * pow( (t-tzero)/Ts, 4.0);

    return v;
  }
