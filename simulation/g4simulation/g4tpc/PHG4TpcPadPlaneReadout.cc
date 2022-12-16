#include "PHG4TpcPadPlaneReadout.h"

#include <g4detectors/PHG4CellDefs.h>  // for genkey, keytype
#include <g4detectors/PHG4TpcCylinderGeom.h>
#include <g4detectors/PHG4TpcCylinderGeomContainer.h>

#include <g4main/PHG4Hit.h>  // for PHG4Hit
#include <g4main/PHG4HitContainer.h>


#include <phool/PHRandomSeed.h>

// Move to new storage containers
#include <trackbase/TpcDefs.h>
#include <trackbase/TrkrDefs.h>  // for hitkey, hitse...
#include <trackbase/TrkrHit.h>   // for TrkrHit
#include <trackbase/TrkrHitSet.h>
#include <trackbase/TrkrHitSetContainer.h>
#include <trackbase/TrkrHitv2.h>  // for TrkrHit

#include <trackbase/TrkrTruthTrack.h>
#include <trackbase/TrkrTruthTrackv1.h>
#include <trackbase/TrkrTruthTrackContainer.h>
#include <trackbase/TrkrTruthTrackContainerv1.h>

#include <phool/phool.h>  // for PHWHERE

#include <TSystem.h>

#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>  // for gsl_rng_alloc

#include <cmath>
#include <iostream>
#include <map>      // for _Rb_tree_cons...
#include <utility>  // for pair

class PHCompositeNode;
class TrkrHitTruthAssoc;

namespace
{
  //! convenient square function
  template <class T>
  inline constexpr T square(const T &x)
  {
    return x * x;
  }

  //! return normalized gaussian centered on zero and of width sigma
  template <class T>
  inline T gaus(const T &x, const T &sigma)
  {
    return std::exp(-square(x / sigma) / 2) / (sigma * std::sqrt(2 * M_PI));
  }

}  // namespace

PHG4TpcPadPlaneReadout::PHG4TpcPadPlaneReadout(const std::string &name)
  : PHG4TpcPadPlane(name)
{
  InitializeParameters();

  RandomGenerator = gsl_rng_alloc(gsl_rng_mt19937);
  gsl_rng_set(RandomGenerator, PHRandomSeed());  // fixed seed is handled in this funtcion

  return;
}

PHG4TpcPadPlaneReadout::~PHG4TpcPadPlaneReadout()
{
  gsl_rng_free(RandomGenerator);
}

int PHG4TpcPadPlaneReadout::CreateReadoutGeometry(PHCompositeNode * /*topNode*/, PHG4TpcCylinderGeomContainer *seggeo)
{
  if (Verbosity()) std::cout << "PHG4TpcPadPlaneReadout: CreateReadoutGeometry: " << std::endl;

  for (int iregion = 0; iregion < 3; ++iregion)
  {
    //int zside = 0;
    for (int zside = 0; zside < 2;zside++){
      sector_R_bias[zside].clear();
      sector_Phi_bias[zside].clear();
      sector_min_Phi[zside].clear();
      sector_max_Phi[zside].clear();
      sector_min_Phi_sectors[zside][iregion].clear();
      sector_max_Phi_sectors[zside][iregion].clear();
      //int eff_layer = 0;
      for (int isector = 0; isector < NSectors; ++isector)//12 sectors
      {
        sector_R_bias[zside].push_back(dR[zside][isector][iregion]);
        sector_Phi_bias[zside].push_back(dPhi[zside][isector][iregion]);
  
        double sec_gap = (2*M_PI - SectorPhi[iregion]*12)/12;
        double sec_max_phi = M_PI - SectorPhi[iregion]/2 - sec_gap - 2 * M_PI / 12 * isector;// * (isector+1) ;
        double sec_min_phi = sec_max_phi - SectorPhi[iregion];
        sector_min_Phi[zside].push_back(sec_min_phi);
        sector_max_Phi[zside].push_back(sec_max_phi);
        sector_min_Phi_sectors[zside][iregion].push_back(sec_min_phi);
        sector_max_Phi_sectors[zside][iregion].push_back(sec_max_phi);
  
      }// isector
    }

    double sum_r = 0;
    for (int layer = MinLayer[iregion]; layer < MinLayer[iregion] + NTpcLayers[iregion]; ++layer)
    {          
      double r_length = Thickness[iregion];
      if(iregion == 0 && layer>0){
        if(layer%2==0) r_length = Thickness[4];
        else r_length = Thickness[3];
      }
      sum_r += r_length;
    }    
    double pad_space = (MaxRadius[iregion] - MinRadius[iregion] - sum_r)/(NTpcLayers[iregion]-1);
    double current_r = MinRadius[iregion];

    for (int layer = MinLayer[iregion]; layer < MinLayer[iregion] + NTpcLayers[iregion]; ++layer)
    {
      if (Verbosity())
      {
        std::cout << " layer " << layer << " MinLayer " << MinLayer[iregion] << " region " << iregion
                  << " radius " << MinRadius[iregion] + ((double) (layer - MinLayer[iregion]) + 0.5) * Thickness[iregion]
                  << " thickness " << Thickness[iregion]
                  << " NTBins " << NTBins << " tmin " << MinT << " tstep " << TBinWidth
                  << " phibins " << NPhiBins[iregion] << " phistep " << PhiBinWidth[iregion] << std::endl;
      }

      PHG4TpcCylinderGeom *layerseggeo = new PHG4TpcCylinderGeom();
      layerseggeo->set_layer(layer);

      //layerseggeo->set_radius(MinRadius[iregion] + ((double) (layer - MinLayer[iregion]) + 0.5) * Thickness[iregion]);
      //layerseggeo->set_thickness(Thickness[iregion]);

      double r_length = Thickness[iregion];
      if(iregion == 0 && layer>0){
        if(layer%2==0) r_length = Thickness[4];
        else r_length = Thickness[3];
      }
      layerseggeo->set_thickness(r_length);
      layerseggeo->set_radius(current_r+r_length/2);
      layerseggeo->set_binning(PHG4CellDefs::sizebinning);
      layerseggeo->set_zbins(NTBins);
      layerseggeo->set_zmin(MinT);
      layerseggeo->set_zstep(TBinWidth);
      layerseggeo->set_phibins(NPhiBins[iregion]);
      layerseggeo->set_phistep(PhiBinWidth[iregion]);
      layerseggeo->set_r_bias(sector_R_bias);
      layerseggeo->set_phi_bias(sector_Phi_bias);
      layerseggeo->set_sector_min_phi(sector_min_Phi);
      layerseggeo->set_sector_max_phi(sector_max_Phi);

      // Chris Pinkenburg: greater causes huge memory growth which causes problems
      // on our farm. If you need to increase this - TALK TO ME first
      if (NPhiBins[iregion] * NTBins > 5100000)
      {
        std::cout << "increase Tpc cellsize, number of cells "
                  << NPhiBins[iregion] * NTBins << " for layer " << layer
                  << " exceed 5.1M limit" << std::endl;
        gSystem->Exit(1);
      }
      seggeo->AddLayerCellGeom(layerseggeo);

      current_r += r_length + pad_space;
    }
  }

  GeomContainer = seggeo;

  return 0;
}

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

  return nelec;
}


void PHG4TpcPadPlaneReadout::MapToPadPlane(
    TpcClusterBuilder  *tpc_truth_clusterer,
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
  if (phi > +M_PI) phi -= 2 * M_PI;
  if (phi < -M_PI) phi += 2 * M_PI;

  rad_gem = sqrt(x_gem * x_gem + y_gem * y_gem);

  // Moving electrons from dead area to a closest pad 
  for (int iregion = 0; iregion < 3; ++iregion)
  {
    double daR = 0; 
    if(iregion==0 || iregion==2){
      daR=1.0;//1.0cm edge to collect electrons from
    }else{
      daR = MinRadius[iregion]-MaxRadius[iregion-1];
    }
    if ( rad_gem <= MinRadius[iregion] && rad_gem >= MinRadius[iregion]-daR){
      if( rad_gem <= MinRadius[iregion]-daR/2){
        rad_gem = MinRadius[iregion] - (1.1*daR) ; 
      }else{
        rad_gem = MinRadius[iregion] + 0.1*daR; 
      }

    }

  }

  phi = check_phi(side, phi, rad_gem);
  unsigned int layernum = 0;
  /* TpcClusterBuilder pass_data {}; */

  // Find which readout layer this electron ends up in

  PHG4TpcCylinderGeomContainer::ConstRange layerrange = GeomContainer->get_begin_end();
  for (PHG4TpcCylinderGeomContainer::ConstIterator layeriter = layerrange.first;
       layeriter != layerrange.second;
       ++layeriter)
  {
    double rad_low  = layeriter->second->get_radius() - layeriter->second->get_thickness() / 2.0;
    double rad_high = layeriter->second->get_radius() + layeriter->second->get_thickness() / 2.0;

    if (rad_gem > rad_low && rad_gem < rad_high)
    {
      // capture the layer where this electron hits the gem stack
      LayerGeom = layeriter->second;
      layernum = LayerGeom->get_layer();
      /* pass_data.layerGeom = LayerGeom; */
      /* pass_data.layer = layernum; */
      if (Verbosity() > 1000)
        std::cout << " g4hit id " << hiter->first << " rad_gem " << rad_gem << " rad_low " << rad_low << " rad_high " << rad_high
                  << " layer  " << hiter->second->get_layer() << " want to change to " << layernum << std::endl;
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

  const auto tbins   = LayerGeom->get_zbins();

  // Create the distribution function of charge on the pad plane around the electron position

  // The resolution due to pad readout includes the charge spread during GEM multiplication.
  // this now defaults to 400 microns during construction from Tom (see 8/11 email).
  // Use the setSigmaT(const double) method to update...
  // We use a double gaussian to represent the smearing due to the SAMPA chip shaping time - default values of fShapingLead and fShapingTail are for 80 ns SAMPA

  // amplify the single electron in the gem stack
  //===============================

  double nelec = getSingleEGEMAmplification();
  /* pass_data.neff_electrons = nelec; */

  // Distribute the charge between the pads in phi
  //====================================

  if (Verbosity() > 200)
    std::cout << "  populate phi bins for "
              << " layernum " << layernum
              << " phi " << phi
              << " sigmaT " << sigmaT
              //<< " zigzag_pads " << zigzag_pads
              << std::endl;

  pad_phibin.clear();
  pad_phibin_share.clear();

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
    pad_phibin_share[iphi] /= norm1;

  // Distribute the charge between the pads in t
  //====================================
  if (Verbosity() > 100 && layernum == print_layer)
    std::cout << "  populate t bins for layernum " << layernum
              << " with t_gem " << t_gem << " sigmaL[0] " << sigmaL[0] << " sigmaL[1] " << sigmaL[1] << std::endl;

  adc_tbin.clear();
  adc_tbin_share.clear();
  populate_tbins(t_gem, sigmaL, adc_tbin, adc_tbin_share);
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
    adc_tbin_share[it] /= tnorm;

  // Fill HitSetContainer
  //===============
  // These are used to do a quick clustering for checking
  double phi_integral = 0.0;
  double t_integral = 0.0;
  double weight = 0.0;

  for (unsigned int ipad = 0; ipad < pad_phibin.size(); ++ipad)
  {
    int    pad_num   = pad_phibin[ipad];
    double pad_share = pad_phibin_share[ipad];

    for (unsigned int it = 0; it < adc_tbin.size(); ++it)
    {
      int tbin_num = adc_tbin[it];
      double adc_bin_share = adc_tbin_share[it];

      // Divide electrons from avalanche between bins
      float neffelectrons = nelec * (pad_share) * (adc_bin_share);
      if (neffelectrons < neffelectrons_threshold) continue;  // skip signals that will be below the noise suppression threshold

      if (tbin_num >= tbins)  std::cout << " Error making key: adc_tbin " << tbin_num << " ntbins " << tbins << std::endl;
      if (pad_num >= phibins) std::cout << " Error making key: pad_phibin " << pad_num << " nphibins " << phibins << std::endl;

      // collect information to do simple clustering. Checks operation of PHG4CylinderCellTpcReco, and
      // is also useful for comparison with PHG4TpcClusterizer result when running single track events.
      // The only information written to the cell other than neffelectrons is tbin and pad number, so get those from geometry
      double tcenter   = LayerGeom->get_zcenter(tbin_num);
      double phicenter = LayerGeom->get_phicenter(pad_num);
      phi_integral += phicenter * neffelectrons;
      t_integral   += tcenter   * neffelectrons;
      weight += neffelectrons;
      if (Verbosity() > 1 && layernum == print_layer)
        std::cout << "   tbin_num " << tbin_num << " tcenter " << tcenter << " pad_num " << pad_num << " phicenter " << phicenter
                  << " neffelectrons " << neffelectrons << " neffelectrons_threshold " << neffelectrons_threshold << std::endl;

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

      // generate the key for this hit, requires tbin and phibin
      TrkrDefs::hitkey hitkey = TpcDefs::genHitKey((unsigned int) pad_num, (unsigned int) tbin_num);
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

      tpc_truth_clusterer->addhitset(hitsetkey, hitkey, neffelectrons);

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
  }    // end of loop over zigzag pads
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

  m_NHits++;
  /* return pass_data; */
}
double PHG4TpcPadPlaneReadout::check_phi(const unsigned int side, const double phi, const double radius){
  double new_phi = phi;
  int p_region=-1;
  for (int iregion = 0; iregion < 3; ++iregion)
  {  
    if (radius < MaxRadius[iregion] && radius > MinRadius[iregion]) p_region = iregion;
  }
    if(p_region>0){
      for(int s=0; s<12;s++){
        double daPhi = 0;
        if (s==0){
          daPhi = fabs(sector_min_Phi_sectors[side][p_region][11] + 2*M_PI - sector_max_Phi_sectors[side][p_region][s]);
        }else{
          daPhi = fabs(sector_min_Phi_sectors[side][p_region][s-1] - sector_max_Phi_sectors[side][p_region][s]);
        }
      double min_phi = sector_max_Phi_sectors[side][p_region][s];
      double max_phi = sector_max_Phi_sectors[side][p_region][s]+daPhi;
        if (new_phi<=max_phi && new_phi>=min_phi){
          if(fabs(max_phi - new_phi) > fabs(new_phi - min_phi)){
            new_phi = min_phi-PhiBinWidth[p_region]/5;//to be changed
          }else{
            new_phi = max_phi+PhiBinWidth[p_region]/5;
          }


         }

      }
    }
    if(new_phi<sector_min_Phi_sectors[side][p_region][11] && new_phi>=-M_PI){
      new_phi += 2*M_PI;
    }
  return new_phi;
}
void PHG4TpcPadPlaneReadout::populate_zigzag_phibins(const unsigned int side, const unsigned int layernum, const double phi, const double cloud_sig_rp, std::vector<int> &phibin_pad, std::vector<double> &phibin_pad_share)
{
  const double radius      = LayerGeom->get_radius();
  const double phistepsize = LayerGeom->get_phistep();
  const auto   phibins     = LayerGeom->get_phibins();

  // make the charge distribution gaussian
  double rphi = phi * radius;
  if (Verbosity() > 100)
    if (LayerGeom->get_layer() == print_layer)
    {
      std::cout << " populate_zigzag_phibins for layer " << layernum << " with radius " << radius << " phi " << phi
                << " rphi " << rphi << " phistepsize " << phistepsize << std::endl;
      std::cout << " fcharge created: radius " << radius << " rphi " << rphi << " cloud_sig_rp " << cloud_sig_rp << std::endl;
    }

  // Get the range of phi values that completely contains all pads  that touch the charge distribution - (nsigmas + 1/2 pad width) in each direction
  const double philim_low_calc  = phi - (_nsigmas * cloud_sig_rp / radius) - phistepsize;
  const double philim_high_calc = phi + (_nsigmas * cloud_sig_rp / radius) + phistepsize;

  // Find the pad range that covers this phi range
  const double philim_low = check_phi(side, philim_low_calc, radius);
  const double philim_high = check_phi(side, philim_high_calc, radius);

  int phibin_low  = LayerGeom->get_phibin(philim_high);
  int phibin_high = LayerGeom->get_phibin(philim_low);
  int npads       = phibin_high - phibin_low;


  if (Verbosity() > 1000)
    if (layernum == print_layer)
      std::cout << "           zigzags: phi " << phi << " philim_low " << philim_low << " phibin_low " << phibin_low
                << " philim_high " << philim_high << " phibin_high " << phibin_high << " npads " << npads << std::endl;

  if (npads < 0 || npads > 9) npads = 9;  // can happen if phibin_high wraps around. If so, limit to 10 pads and fix below

  // Calculate the maximum extent in r-phi of pads in this layer. Pads are assumed to touch the center of the next phi bin on both sides.
  const double pad_rphi = 2.0 * LayerGeom->get_phistep() * radius;

  // Make a TF1 for each pad in the phi range
  using PadParameterSet = std::array<double, 2>;
  std::array<PadParameterSet, 10> pad_parameters;
  std::array<int, 10> pad_keep;


  // Now make a loop that steps through the charge distribution and evaluates the response at that point on each pad
  std::array<double, 10> overlap = {{0, 0, 0, 0, 0, 0, 0, 0, 0, 0}};
  double pads_phi[10] = {0., 0., 0., 0., 0., 0., 0., 0., 0., 0.};
  double sum_of_pads_phi = 0.;
  double sum_of_pads_absphi = 0.;
  for (int ipad = 0; ipad <= npads; ipad++)
  {
    int pad_now = phibin_low + ipad;
    if (pad_now >= phibins) pad_now -= phibins;
    pads_phi[ipad] = LayerGeom->get_phicenter(pad_now);
    sum_of_pads_phi += pads_phi[ipad];
    sum_of_pads_absphi += fabs(pads_phi[ipad]);
  }

  for (int ipad = 0; ipad <= npads; ipad++)
  {
    int pad_now = phibin_low + ipad;
    //if(phibin_low<0 && phibin_high<0) pad_now = phibin_high + ipad;
    // check that we do not exceed the maximum number of pads, wrap if necessary
    if (pad_now >= phibins) pad_now -= phibins;

    pad_keep[ipad] = pad_now;
    double phi_now = pads_phi[ipad];
    const double rphi_pad_now = phi_now * radius;
    pad_parameters[ipad] = {{pad_rphi / 2.0, rphi_pad_now}};

    if (Verbosity() > 1000)
      if (layernum == print_layer)
        std::cout << " zigzags: make fpad for ipad " << ipad << " pad_now " << pad_now << " pad_rphi/2 " << pad_rphi / 2.0
                  << " rphi_pad_now " << rphi_pad_now << std::endl;
  //}

  // use analytic integral
  //for (int ipad = 0; ipad <= npads; ipad++)
  //{
    //const double pitch = pad_parameters[ipad][0];
    //const double x_loc = pad_parameters[ipad][1] - rphi;
    //const double sigma = cloud_sig_rp;

    const double pitch = pad_rphi / 2.0;//eshulga
    double x_loc_tmp = rphi_pad_now - rphi;//eshulga
    const double sigma = cloud_sig_rp;//eshulga

    // Checking if the pads are on the same side of the TPC in phi 
    if(fabs(sum_of_pads_phi)!= sum_of_pads_absphi){
      if(phi<-M_PI/2 && phi_now>0){
        x_loc_tmp = (phi_now - 2*M_PI) * radius - rphi;
      }
      if(phi>M_PI/2 && phi_now<0){
        x_loc_tmp = (phi_now + 2*M_PI) * radius - rphi;
      }
      if(phi<0 && phi_now>0){
        x_loc_tmp = (phi_now+fabs(phi)) * radius;
      }
      if(phi>0 && phi_now<0){
        x_loc_tmp = (2*M_PI - phi_now + phi) * radius;
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
    if (rad_gem < output_radius) std::cout << "         zigzags: for pad " << ipad << " integral is " << overlap[ipad] << std::endl;
   }

  return;
}

void PHG4TpcPadPlaneReadout::populate_tbins(const double t, const std::array<double, 2> &cloud_sig_tt, std::vector<int> &tbin_adc, std::vector<double> &tbin_adc_share)
{
  int tbin = LayerGeom->get_zbin(t);
  if (tbin < 0 || tbin > LayerGeom->get_zbins())
  {
    //std::cout << " t bin is outside range, return" << std::endl;
    return;
  }

  double tstepsize = LayerGeom->get_zstep();
  double tdisp = t - LayerGeom->get_zcenter(tbin);

  if (Verbosity() > 1000)
    std::cout << "     input:  t " << t << " tbin " << tbin << " tstepsize " << tstepsize << " t center " << LayerGeom->get_zcenter(tbin) << " tdisp " << tdisp << std::endl;

  // Because of diffusion, hits can be shared across the membrane, so we allow all t bins
  int min_cell_tbin = 0;
  int max_cell_tbin = NTBins - 1;

  double cloud_sig_tt_inv[2];
  cloud_sig_tt_inv[0] = 1. / cloud_sig_tt[0];
  cloud_sig_tt_inv[1] = 1. / cloud_sig_tt[1];

  int zsect = 0;
  if (t < 0)
    zsect = -1;
  else
    zsect = 1;

  int n_zz = int(3 * (cloud_sig_tt[0] + cloud_sig_tt[1]) / (2.0 * tstepsize) + 1);
  if (Verbosity() > 1000) std::cout << " n_zz " << n_zz << " cloud_sigzz[0] " << cloud_sig_tt[0] << " cloud_sig_tt[1] " << cloud_sig_tt[1] << std::endl;
  for (int it = -n_zz; it != n_zz + 1; ++it)
  {
    int cur_t_bin = tbin + it;
    if ((cur_t_bin < min_cell_tbin) || (cur_t_bin > max_cell_tbin)) continue;

    if (Verbosity() > 1000)
      std::cout << " it " << it << " cur_t_bin " << cur_t_bin << " min_cell_tbin " << min_cell_tbin << " max_cell_tbin " << max_cell_tbin << std::endl;

    double t_integral = 0.0;
    if (it == 0)
    {
      // the crossover between lead and tail shaping occurs in this bin
      int index1 = -1;
      int index2 = -1;
      if (zsect == -1)
      {
        index1 = 0;
        index2 = 1;
      }
      else
      {
        index1 = 1;
        index2 = 0;
      }

      double tLim1 = 0.0;
      double tLim2 = 0.5 * M_SQRT2 * (-0.5 * tstepsize - tdisp) * cloud_sig_tt_inv[index1];
      // 1/2 * the erf is the integral probability from the argument Z value to zero, so this is the integral probability between the Z limits
      double t_integral1 = 0.5 * (erf(tLim1) - erf(tLim2));

      if (Verbosity() > 1000)
        if (LayerGeom->get_layer() == print_layer)
          std::cout << "   populate_tbins:  cur_t_bin " << cur_t_bin << "  center t " << LayerGeom->get_zcenter(cur_t_bin)
                    << " index1 " << index1 << "  tLim1 " << tLim1 << " tLim2 " << tLim2 << " t_integral1 " << t_integral1 << std::endl;

      tLim2 = 0.0;
      tLim1 = 0.5 * M_SQRT2 * (0.5 * tstepsize - tdisp) * cloud_sig_tt_inv[index2];
      double t_integral2 = 0.5 * (erf(tLim1) - erf(tLim2));

      if (Verbosity() > 1000)
        if (LayerGeom->get_layer() == print_layer)
          std::cout << "   populate_tbins:  cur_t_bin " << cur_t_bin << "  center t " << LayerGeom->get_zcenter(cur_t_bin)
                    << " index2 " << index2 << "  tLim1 " << tLim1 << " tLim2 " << tLim2 << " t_integral2 " << t_integral2 << std::endl;

      t_integral = t_integral1 + t_integral2;
    }
    else
    {
      // The non zero bins are entirely in the lead or tail region
      // lead or tail depends on which side of the membrane
      int index = 0;
      if (it < 0)
      {
        if (zsect == -1)
          index = 0;
        else
          index = 1;
      }
      else
      {
        if (zsect == -1)
          index = 1;
        else
          index = 0;
      }
      double tLim1 = 0.5 * M_SQRT2 * ((it + 0.5) * tstepsize - tdisp) * cloud_sig_tt_inv[index];
      double tLim2 = 0.5 * M_SQRT2 * ((it - 0.5) * tstepsize - tdisp) * cloud_sig_tt_inv[index];
      t_integral = 0.5 * (erf(tLim1) - erf(tLim2));

      if (Verbosity() > 1000)
        if (LayerGeom->get_layer() == print_layer)
          std::cout << "   populate_tbins:  t_bin " << cur_t_bin << "  center t " << LayerGeom->get_zcenter(cur_t_bin)
                    << " index " << index << "  tLim1 " << tLim1 << " tLim2 " << tLim2 << " t_integral " << t_integral << std::endl;
    }

    tbin_adc.push_back(cur_t_bin);
    tbin_adc_share.push_back(t_integral);
  }

  return;
}

void PHG4TpcPadPlaneReadout::SetDefaultParameters()
{
  set_default_int_param("ntpc_layers_inner", 16);
  set_default_int_param("ntpc_layers_mid", 16);
  set_default_int_param("ntpc_layers_outer", 16);

  set_default_int_param("tpc_minlayer_inner", 7);

  //set_default_double_param("tpc_minradius_inner", 30.0);  // cm
  //set_default_double_param("tpc_minradius_mid", 40.0);
  //set_default_double_param("tpc_minradius_outer", 60.0);
//
  //set_default_double_param("tpc_maxradius_inner", 40.0);  // cm
  //set_default_double_param("tpc_maxradius_mid", 60.0);
  //set_default_double_param("tpc_maxradius_outer", 77.0);  // from Tom

  set_default_double_param("tpc_minradius_inner", 31.105);//30.0);  // cm
  set_default_double_param("tpc_minradius_mid", 41.153);//40.0);
  set_default_double_param("tpc_minradius_outer", 58.367);//60.0);


  set_default_double_param("tpc_maxradius_inner", 40.249);//40.0);  // cm
  set_default_double_param("tpc_maxradius_mid", 57.475);//60.0);
  set_default_double_param("tpc_maxradius_outer", 75.911);//77.0);  // from Tom


  set_default_double_param("neffelectrons_threshold", 1.0);
  set_default_double_param("maxdriftlength", 105.5);     // cm
  set_default_double_param("tpc_adc_clock", 53.0);       // ns, for 18.8 MHz clock
  set_default_double_param("gem_cloud_sigma", 0.04);     // cm = 400 microns
  set_default_double_param("sampa_shaping_lead", 32.0);  // ns, for 80 ns SAMPA
  set_default_double_param("sampa_shaping_tail", 48.0);  // ns, for 80 ns SAMPA

  set_default_double_param("tpc_sector_phi_inner", 0.5024);//2 * M_PI / 12 );//sector size in phi for R1 sector
  set_default_double_param("tpc_sector_phi_mid",   0.5087);//2 * M_PI / 12 );//sector size in phi for R2 sector
  set_default_double_param("tpc_sector_phi_outer", 0.5097);//2 * M_PI / 12 );//sector size in phi for R3 sector

  set_default_int_param("ntpc_phibins_inner", 1152);
  set_default_int_param("ntpc_phibins_mid", 1536);
  set_default_int_param("ntpc_phibins_outer", 2304);

  // GEM Gain
  /*
  hp (2020/09/04): gain changed from 2000 to 1400, to accomodate gas mixture change 
  from Ne/CF4 90/10 to Ne/CF4 50/50, and keep the average charge per particle per pad constant
  */
  set_default_double_param("gem_amplification", 1400);
  return;
}

void PHG4TpcPadPlaneReadout::UpdateInternalParameters()
{
  NTpcLayers[0] = get_int_param("ntpc_layers_inner");
  NTpcLayers[1] = get_int_param("ntpc_layers_mid");
  NTpcLayers[2] = get_int_param("ntpc_layers_outer");

  MinLayer[0] = get_int_param("tpc_minlayer_inner");
  MinLayer[1] = MinLayer[0] + NTpcLayers[0];
  MinLayer[2] = MinLayer[1] + NTpcLayers[1];

  neffelectrons_threshold = get_double_param("neffelectrons_threshold");

  MinRadius[0] = get_double_param("tpc_minradius_inner");
  MinRadius[1] = get_double_param("tpc_minradius_mid");
  MinRadius[2] = get_double_param("tpc_minradius_outer");

  MaxRadius[0] = get_double_param("tpc_maxradius_inner");
  MaxRadius[1] = get_double_param("tpc_maxradius_mid");
  MaxRadius[2] = get_double_param("tpc_maxradius_outer");

  //Thickness[0] = NTpcLayers[0] <= 0 ? 0 : (MaxRadius[0] - MinRadius[0]) / NTpcLayers[0];
  //Thickness[1] = NTpcLayers[1] <= 0 ? 0 : (MaxRadius[1] - MinRadius[1]) / NTpcLayers[1];
  //Thickness[2] = NTpcLayers[2] <= 0 ? 0 : (MaxRadius[2] - MinRadius[2]) / NTpcLayers[2];

  MaxRadius[0] = get_double_param("tpc_maxradius_inner");
  MaxRadius[1] = get_double_param("tpc_maxradius_mid");
  MaxRadius[2] = get_double_param("tpc_maxradius_outer");

  Thickness[0] = 0.687;
  Thickness[1] = 1.012;
  Thickness[2] = 1.088;
  Thickness[3] = 0.534;
  Thickness[4] = 0.595;

  sigmaT = get_double_param("gem_cloud_sigma");
  sigmaL[0] = get_double_param("sampa_shaping_lead");
  sigmaL[1] = get_double_param("sampa_shaping_tail");

  tpc_adc_clock = get_double_param("tpc_adc_clock");

  MaxZ = get_double_param("maxdriftlength");
  TBinWidth = tpc_adc_clock;
  ZBinWidth = TBinWidth * drift_velocity;
  MaxT = 2.0 * MaxZ / drift_velocity;  // allows for extended time readout
  MinT = 0;
  NTBins = (int) ((MaxT - MinT) / TBinWidth) + 1;

  std::cout << PHWHERE << "MaxT " << MaxT << " NTBins = " << NTBins << " drift velocity " << drift_velocity << std::endl;

  SectorPhi[0] = get_double_param("tpc_sector_phi_inner");
  SectorPhi[1] = get_double_param("tpc_sector_phi_mid");
  SectorPhi[2] = get_double_param("tpc_sector_phi_outer");

  NPhiBins[0] = get_int_param("ntpc_phibins_inner");
  NPhiBins[1] = get_int_param("ntpc_phibins_mid");
  NPhiBins[2] = get_int_param("ntpc_phibins_outer");

  PhiBinWidth[0] = SectorPhi[0] * 12 / (double) NPhiBins[0];
  PhiBinWidth[1] = SectorPhi[1] * 12 / (double) NPhiBins[1];
  PhiBinWidth[2] = SectorPhi[2] * 12 / (double) NPhiBins[2];

  averageGEMGain = get_double_param("gem_amplification");

  //The modules are segmented in [-M_PI;M_PI] interval
  std::array< std::array< std::array< float,NRSectors >,NSectors >,NSides > dR_tmp = {{
    {{
      {0.,0.,0.},
      {0.,0.,0.},
      {0.,0.,0.},
      {0.,0.,0.},
      {0.,0.,0.},
      {0.,0.,0.},
      {0.,0.,0.},
      {0.,0.,0.},
      {0.,0.,0.},
      {0.,0.,0.},
      {0.,0.,0.},
      {0.,0.,0.}
    }},
    {{
      {0.,0.,0.},
      {0.,0.,0.},
      {0.,0.,0.},
      {0.,0.,0.},
      {0.,0.,0.},
      {0.,0.,0.},
      {0.,0.,0.},
      {0.,0.,0.},
      {0.,0.,0.},
      {0.,0.,0.},
      {0.,0.,0.},
      {0.,0.,0.}
    }}
  }};

  std::array< std::array< std::array< float,NRSectors >,NSectors >,NSides > dPhi_tmp = {{
    {{
      {0.,0.,0.},
      {0.,0.,0.},
      {0.,0.,0.},
      {0.,0.,0.},
      {0.,0.,0.},
      {0.,0.,0.},
      {0.,0.,0.},
      {0.,0.,0.},
      {0.,0.,0.},
      {0.,0.,0.},
      {0.,0.,0.},
      {0.,0.,0.}
    }},
    {{
      {0.,0.,0.},
      {0.,0.,0.},
      {0.,0.,0.},
      {0.,0.,0.},
      {0.,0.,0.},
      {0.,0.,0.},
      {0.,0.,0.},
      {0.,0.,0.},
      {0.,0.,0.},
      {0.,0.,0.},
      {0.,0.,0.},
      {0.,0.,0.}
    }}
  }};

  dR = dR_tmp;
  dPhi = dPhi_tmp;
}
