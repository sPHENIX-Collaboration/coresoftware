#include "PHG4TpcPadPlaneReadout.h"

#include <g4detectors/PHG4Cell.h>                       // for PHG4Cell
#include <g4detectors/PHG4CellDefs.h>                   // for genkey, keytype
#include <g4detectors/PHG4CellContainer.h>
#include <g4detectors/PHG4Cellv1.h>
#include <g4detectors/PHG4CylinderCellGeom.h>
#include <g4detectors/PHG4CylinderCellGeomContainer.h>

#include <g4main/PHG4Hit.h>                             // for PHG4Hit
#include <g4main/PHG4HitContainer.h>
#include <phool/PHRandomSeed.h>

// Move to new storage containers
#include <trackbase/TrkrDefs.h>                         // for hitkey, hitse...
#include <trackbase/TrkrHit.h>                          // for TrkrHit
#include <trackbase/TrkrHitSet.h>
#include <trackbase/TrkrHitSetContainer.h>

#include <tpc/TpcDefs.h>
#include <tpc/TpcHit.h>

#include <TF1.h>
#include <TNtuple.h>
#include <TSystem.h>

#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>                                // for gsl_rng_alloc

#include <cmath>
#include <climits>                                     // for INT_MAX
#include <cstdio>                                      // for sprintf
#include <iostream>
#include <map>                                          // for _Rb_tree_cons...
#include <utility>                                      // for pair

class PHCompositeNode;
class TrkrHitTruthAssoc;

using namespace std;

namespace
{

  //______________________________________________________________________
  template<class T> constexpr T square( const T& x ) { return x*x; }

  //______________________________________________________________________
  double get_pad_response( double position, const std::array<double,2>& par )
  { return std::get<0>(par)-std::abs(position - std::get<1>(par)); }

}

PHG4TpcPadPlaneReadout::PHG4TpcPadPlaneReadout(const string &name)
  : PHG4TpcPadPlane(name)
{
  InitializeParameters();

  // initialize gaussian weights
  for( int i = 0; i < _ngauss_steps; ++i )
  {
    const double x = -4.5 + 2.0*_nsigmas*i/_ngauss_steps;
    _gauss_weights[i] = std::exp( -square( x )/2 );
  }

  RandomGenerator = gsl_rng_alloc(gsl_rng_mt19937);
  gsl_rng_set(RandomGenerator, PHRandomSeed());  // fixed seed is handled in this funtcion

  return;
}

PHG4TpcPadPlaneReadout::~PHG4TpcPadPlaneReadout()
{
  gsl_rng_free(RandomGenerator);
}

int PHG4TpcPadPlaneReadout::CreateReadoutGeometry(PHCompositeNode *topNode, PHG4CylinderCellGeomContainer *seggeo)
{
  if (Verbosity()) cout << "PHG4TpcPadPlaneReadout: CreateReadoutGeometry: " << endl;

  for (int iregion = 0; iregion < 3; ++iregion)
  {
    for (int layer = MinLayer[iregion]; layer < MinLayer[iregion] + NTpcLayers[iregion]; ++layer)
    {
      if (Verbosity())
      {
        cout << " layer " << layer << " MinLayer " << MinLayer[iregion] << " region " << iregion
             << " radius " << MinRadius[iregion] + ((double) (layer - MinLayer[iregion]) + 0.5) * Thickness[iregion]
             << " thickness " << Thickness[iregion]
             << " NZbins " << NZBins << " zmin " << MinZ << " zstep " << ZBinWidth
             << " phibins " << NPhiBins[iregion] << " phistep " << PhiBinWidth[iregion] << endl;
      }

      PHG4CylinderCellGeom *layerseggeo = new PHG4CylinderCellGeom();
      layerseggeo->set_layer(layer);
      layerseggeo->set_radius(MinRadius[iregion] + ((double) (layer - MinLayer[iregion]) + 0.5) * Thickness[iregion]);
      layerseggeo->set_thickness(Thickness[iregion]);
      layerseggeo->set_binning(PHG4CellDefs::sizebinning);
      layerseggeo->set_zbins(NZBins);
      layerseggeo->set_zmin(MinZ);
      layerseggeo->set_zstep(ZBinWidth);
      layerseggeo->set_phibins(NPhiBins[iregion]);
      layerseggeo->set_phistep(PhiBinWidth[iregion]);
      // Chris Pinkenburg: greater causes huge memory growth which causes problems
      // on our farm. If you need to increase this - TALK TO ME first
      if (NPhiBins[iregion] * NZBins > 5100000)
      {
        cout << "increase Tpc cellsize, number of cells "
             << NPhiBins[iregion] * NZBins << " for layer " << layer
             << " exceed 5.1M limit" << endl;
        gSystem->Exit(1);
      }
      seggeo->AddLayerCellGeom(layerseggeo);
    }
  }

  GeomContainer = seggeo;

  return 0;
}

// This is obsolete, it uses the old PHG4Cell containers
void PHG4TpcPadPlaneReadout::MapToPadPlane(PHG4CellContainer *g4cells, const double x_gem, const double y_gem, const double z_gem, PHG4HitContainer::ConstIterator hiter, TNtuple *ntpad, TNtuple *nthit)
{
  // One electron per call of this method
  // The x_gem and y_gem values have already been randomized within the transverse drift diffusion width
  // The z_gem value already reflects the drift time of the primary electron from the production point, and is randomized within the longitudinal diffusion witdth

  double phi = atan2(y_gem, x_gem);
  if (phi > +M_PI) phi -= 2 * M_PI;
  if (phi < -M_PI) phi += 2 * M_PI;

  rad_gem = sqrt(x_gem * x_gem + y_gem * y_gem);
  //cout << "Enter old MapToPadPlane with rad_gem " << rad_gem << endl;

  unsigned int layernum = 0;

  // Find which readout layer this electron ends up in

  PHG4CylinderCellGeomContainer::ConstRange layerrange = GeomContainer->get_begin_end();
  for (PHG4CylinderCellGeomContainer::ConstIterator layeriter = layerrange.first;
       layeriter != layerrange.second;
       ++layeriter)
  {
    double rad_low = layeriter->second->get_radius() - layeriter->second->get_thickness() / 2.0;
    double rad_high = layeriter->second->get_radius() + layeriter->second->get_thickness() / 2.0;

    if (rad_gem > rad_low && rad_gem < rad_high)
    {
      // capture the layer where this electron hits sthe gem stack
      LayerGeom = layeriter->second;
      layernum = LayerGeom->get_layer();
      if (Verbosity() > 1000 && layernum == print_layer)
        cout << " g4hit id " << hiter->first << " rad_gem " << rad_gem << " rad_low " << rad_low << " rad_high " << rad_high
             << " layer  " << hiter->second->get_layer() << " want to change to " << layernum << endl;
      hiter->second->set_layer(layernum);  // have to set here, since the stepping action knows nothing about layers
    }
  }

  if (layernum == 0)
  {
    return;
  }

  // store phi bins and zbins upfront to avoid repetitive checks on the phi methods
  const auto phibins = LayerGeom->get_phibins();
  const auto zbins = LayerGeom->get_zbins();

  // Create the distribution function of charge on the pad plane around the electron position

  // The resolution due to pad readout includes the charge spread during GEM multiplication.
  // this now defaults to 400 microns during construction from Tom (see 8/11 email).
  // Use the setSigmaT(const double) method to update...
  // We use a double gaussian to represent the smearing due to the SAMPA chip shaping time - default values of fShapingLead and fShapingTail are for 80 ns SAMPA

  // amplify the single electron in the gem stack
  //===============================

  double nelec = getSingleEGEMAmplification();


  // Distribute the charge between the pads in phi
  //====================================

  if (Verbosity() > 200)
    cout << "  populate phi bins for "
         << " layernum " << layernum
         << " phi " << phi
         << " sigmaT " << sigmaT
         << " zigzag_pads " << zigzag_pads
         << endl;

  pad_phibin.clear();
  pad_phibin_share.clear();
  if (zigzag_pads)
    populate_zigzag_phibins(layernum, phi, sigmaT, pad_phibin, pad_phibin_share);
  else
    populate_rectangular_phibins(layernum, phi, sigmaT, pad_phibin, pad_phibin_share);

  // Normalize the shares so they add up to 1
  double norm1 = 0.0;
  for (unsigned int ipad = 0; ipad < pad_phibin.size(); ++ipad)
  {
    double pad_share = pad_phibin_share[ipad];
    norm1 += pad_share;
  }
  for (unsigned int iphi = 0; iphi < pad_phibin.size(); ++iphi)
    pad_phibin_share[iphi] /= norm1;

  // Distribute the charge between the pads in z
  //====================================
  if (Verbosity() > 100 && layernum == print_layer)
    cout << "  populate z bins for layernum " << layernum
         << " with z_gem " << z_gem << " sigmaL[0] " << sigmaL[0] << " sigmaL[1] " << sigmaL[1] << endl;

  adc_zbin.clear();
  adc_zbin_share.clear();
  populate_zbins(z_gem, sigmaL, adc_zbin, adc_zbin_share);

  // Normalize the shares so that they add up to 1
  double znorm = 0.0;
  for (unsigned int iz = 0; iz < adc_zbin.size(); ++iz)
  {
    double bin_share = adc_zbin_share[iz];
    znorm += bin_share;
  }
  for (unsigned int iz = 0; iz < adc_zbin.size(); ++iz)
    adc_zbin_share[iz] /= znorm;

  // Fill cells
  //========
  // These are used to do a quick clustering for checking
  double phi_integral = 0.0;
  double z_integral = 0.0;
  double weight = 0.0;

  for (unsigned int ipad = 0; ipad < pad_phibin.size(); ++ipad)
  {
    int pad_num = pad_phibin[ipad];
    double pad_share = pad_phibin_share[ipad];

    for (unsigned int iz = 0; iz < adc_zbin.size(); ++iz)
    {
      int zbin_num = adc_zbin[iz];
      double adc_bin_share = adc_zbin_share[iz];

      // Divide electrons from avalanche between bins
      float neffelectrons = nelec * (pad_share) * (adc_bin_share);
      if (neffelectrons < neffelectrons_threshold) continue;  // skip signals that will be below the noise suppression threshold

      if (zbin_num >= zbins) cout << " Error making key: adc_zbin " << zbin_num << " nzbins " << zbins << endl;
      if (pad_num >= phibins) cout << " Error making key: pad_phibin " << pad_num << " nphibins " << phibins << endl;

      // collect information to do simple clustering. Checks operation of PHG4CylinderCellTpcReco, and
      // is also useful for comparison with PHG4TpcClusterizer result when running single track events.
      // The only information written to the cell other than neffelectrons is zbin and pad number, so get those from geometry
      double zcenter = LayerGeom->get_zcenter(zbin_num);
      double phicenter = LayerGeom->get_phicenter(pad_num);
      phi_integral += phicenter * neffelectrons;
      z_integral += zcenter * neffelectrons;
      weight += neffelectrons;
      if (Verbosity() > 1000 && layernum == print_layer)
        cout << "   zbin_num " << zbin_num << " zcenter " << zcenter << " pad_num " << pad_num << " phicenter " << phicenter
             << " neffelectrons " << neffelectrons << " neffelectrons_threshold " << neffelectrons_threshold << endl;

      // Add edep from this electron for this zbin and phi bin combination to the appropriate cell
      PHG4CellDefs::keytype key = PHG4CellDefs::SizeBinning::genkey(layernum, zbin_num, pad_num);
      PHG4Cell *cell = g4cells->findCell(key);
      if (!cell)
      {
        cell = new PHG4Cellv1(key);
        g4cells->AddCell(cell);
      }
      cell->add_edep(neffelectrons);
      cell->add_edep(hiter->first, neffelectrons);  // associates g4hit with this edep
      if (Verbosity() > 100 && layernum == 50) cell->identify();
    }  // end of loop over adc Z bins
  }    // end of loop over zigzag pads

  // Capture the input values at the gem stack and the quick clustering results, elecron-by-electron
  if (Verbosity() > 0)
    ntpad->Fill(layernum, phi, phi_integral / weight, z_gem, z_integral / weight);

  if (Verbosity() > 100)
    if (layernum == print_layer)
    {
      cout << " hit " << hit << " quick centroid for this electron " << endl;
      cout << "      phi centroid = " << phi_integral / weight << " phi in " << phi << " phi diff " << phi_integral / weight - phi << endl;
      cout << "      z centroid = " << z_integral / weight << " z in " << z_gem << " z diff " << z_integral / weight - z_gem << endl;
      // For a single track event, this captures the distribution of single electron centroids on the pad plane for layer print_layer.
      // The centroid of that should match the cluster centroid found by PHG4TpcClusterizer for layer print_layer, if everything is working
      //   - matches to < .01 cm for a few cases that I checked
      nthit->Fill(hit, layernum, phi, phi_integral / weight, z_gem, z_integral / weight, weight);
    }

  hit++;

  return;
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

void PHG4TpcPadPlaneReadout::MapToPadPlane(TrkrHitSetContainer *hitsetcontainer, TrkrHitTruthAssoc *hittruthassoc, const double x_gem, const double y_gem, const double z_gem, PHG4HitContainer::ConstIterator hiter, TNtuple *ntpad, TNtuple *nthit)
{
  // One electron per call of this method
  // The x_gem and y_gem values have already been randomized within the transverse drift diffusion width
  // The z_gem value already reflects the drift time of the primary electron from the production point, and is randomized within the longitudinal diffusion witdth

  double phi = atan2(y_gem, x_gem);
  if (phi > +M_PI) phi -= 2 * M_PI;
  if (phi < -M_PI) phi += 2 * M_PI;

  rad_gem = sqrt(x_gem * x_gem + y_gem * y_gem);
  //cout << "Enter new MapToPadPlane with rad_gem " << rad_gem << endl;

  unsigned int layernum = 0;

  // Find which readout layer this electron ends up in

  PHG4CylinderCellGeomContainer::ConstRange layerrange = GeomContainer->get_begin_end();
  for (PHG4CylinderCellGeomContainer::ConstIterator layeriter = layerrange.first;
       layeriter != layerrange.second;
       ++layeriter)
  {
    double rad_low = layeriter->second->get_radius() - layeriter->second->get_thickness() / 2.0;
    double rad_high = layeriter->second->get_radius() + layeriter->second->get_thickness() / 2.0;

    if (rad_gem > rad_low && rad_gem < rad_high)
    {
      // capture the layer where this electron hits sthe gem stack
      LayerGeom = layeriter->second;
      layernum = LayerGeom->get_layer();
      if (Verbosity() > 1000)
        cout << " g4hit id " << hiter->first << " rad_gem " << rad_gem << " rad_low " << rad_low << " rad_high " << rad_high
             << " layer  " << hiter->second->get_layer() << " want to change to " << layernum << endl;
      hiter->second->set_layer(layernum);  // have to set here, since the stepping action knows nothing about layers
    }
  }

  if (layernum == 0)
  {
    return;
  }

  // store phi bins and zbins upfront to avoid repetitive checks on the phi methods
  const auto phibins = LayerGeom->get_phibins();
  const auto zbins = LayerGeom->get_zbins();

  // Create the distribution function of charge on the pad plane around the electron position

  // The resolution due to pad readout includes the charge spread during GEM multiplication.
  // this now defaults to 400 microns during construction from Tom (see 8/11 email).
  // Use the setSigmaT(const double) method to update...
  // We use a double gaussian to represent the smearing due to the SAMPA chip shaping time - default values of fShapingLead and fShapingTail are for 80 ns SAMPA

  // amplify the single electron in the gem stack
  //===============================

  double nelec = getSingleEGEMAmplification();

  // Distribute the charge between the pads in phi
  //====================================

  if (Verbosity() > 200)
    cout << "  populate phi bins for "
         << " layernum " << layernum
         << " phi " << phi
         << " sigmaT " << sigmaT
         << " zigzag_pads " << zigzag_pads
         << endl;

  pad_phibin.clear();
  pad_phibin_share.clear();
  if (zigzag_pads)
    populate_zigzag_phibins(layernum, phi, sigmaT, pad_phibin, pad_phibin_share);
  else
    populate_rectangular_phibins(layernum, phi, sigmaT, pad_phibin, pad_phibin_share);

  // Normalize the shares so they add up to 1
  double norm1 = 0.0;
  for (unsigned int ipad = 0; ipad < pad_phibin.size(); ++ipad)
  {
    double pad_share = pad_phibin_share[ipad];
    norm1 += pad_share;
  }
  for (unsigned int iphi = 0; iphi < pad_phibin.size(); ++iphi)
    pad_phibin_share[iphi] /= norm1;

  // Distribute the charge between the pads in z
  //====================================
  if (Verbosity() > 100 && layernum == print_layer)
    cout << "  populate z bins for layernum " << layernum
         << " with z_gem " << z_gem << " sigmaL[0] " << sigmaL[0] << " sigmaL[1] " << sigmaL[1] << endl;

  adc_zbin.clear();
  adc_zbin_share.clear();
  populate_zbins(z_gem, sigmaL, adc_zbin, adc_zbin_share);

  // Normalize the shares so that they add up to 1
  double znorm = 0.0;
  for (unsigned int iz = 0; iz < adc_zbin.size(); ++iz)
  {
    double bin_share = adc_zbin_share[iz];
    znorm += bin_share;
  }
  for (unsigned int iz = 0; iz < adc_zbin.size(); ++iz)
    adc_zbin_share[iz] /= znorm;

  // Fill HitSetContainer
  //===============
  // These are used to do a quick clustering for checking
  double phi_integral = 0.0;
  double z_integral = 0.0;
  double weight = 0.0;

  for (unsigned int ipad = 0; ipad < pad_phibin.size(); ++ipad)
  {
    int pad_num = pad_phibin[ipad];
    double pad_share = pad_phibin_share[ipad];

    for (unsigned int iz = 0; iz < adc_zbin.size(); ++iz)
    {
      int zbin_num = adc_zbin[iz];
      double adc_bin_share = adc_zbin_share[iz];

      // Divide electrons from avalanche between bins
      float neffelectrons = nelec * (pad_share) * (adc_bin_share);
      if (neffelectrons < neffelectrons_threshold) continue;  // skip signals that will be below the noise suppression threshold

      if (zbin_num >= zbins) cout << " Error making key: adc_zbin " << zbin_num << " nzbins " << zbins << endl;
      if (pad_num >= phibins) cout << " Error making key: pad_phibin " << pad_num << " nphibins " << phibins << endl;

      // collect information to do simple clustering. Checks operation of PHG4CylinderCellTpcReco, and
      // is also useful for comparison with PHG4TpcClusterizer result when running single track events.
      // The only information written to the cell other than neffelectrons is zbin and pad number, so get those from geometry
      double zcenter = LayerGeom->get_zcenter(zbin_num);
      double phicenter = LayerGeom->get_phicenter(pad_num);
      phi_integral += phicenter * neffelectrons;
      z_integral += zcenter * neffelectrons;
      weight += neffelectrons;
      if (Verbosity() > 1000 && layernum == print_layer)
        cout << "   zbin_num " << zbin_num << " zcenter " << zcenter << " pad_num " << pad_num << " phicenter " << phicenter
             << " neffelectrons " << neffelectrons << " neffelectrons_threshold " << neffelectrons_threshold << endl;

      // new containers
      //============
      // We add the Tpc TrkrHitsets directly to the node using hitsetcontainer
      // We need to create the TrkrHitSet if not already made - each TrkrHitSet should correspond to a Tpc readout module
      // The hitset key includes the layer, sector, side

      // Get the side - 0 for negative z, 1 for positive z
      unsigned int side = 0;
      if (zcenter > 0) side = 1;
      // get the Tpc readout sector - there are 12 sectors with how many pads each?
      unsigned int pads_per_sector = phibins / 12;
      unsigned int sector = pad_num / pads_per_sector;
      TrkrDefs::hitsetkey hitsetkey = TpcDefs::genHitSetKey(layernum, sector, side);
      // Use existing hitset or add new one if needed
      TrkrHitSetContainer::Iterator hitsetit = hitsetcontainer->findOrAddHitSet(hitsetkey);

      // generate the key for this hit, requires zbin and phibin
      TrkrDefs::hitkey hitkey = TpcDefs::genHitKey((unsigned int) pad_num, (unsigned int) zbin_num);
      // See if this hit already exists
      TrkrHit *hit = nullptr;
      hit = hitsetit->second->getHit(hitkey);
      if (!hit)
      {
        // create a new one
        hit = new TpcHit();
        hitsetit->second->addHitSpecificKey(hitkey, hit);
      }

      // Either way, add the energy to it  -- adc values will be added at digitization
      hit->addEnergy(neffelectrons);

      nthit->Fill(layernum, pad_num, zbin_num, neffelectrons);

    }  // end of loop over adc Z bins
  }    // end of loop over zigzag pads

  // Capture the input values at the gem stack and the quick clustering results, elecron-by-electron
  if (Verbosity() > 0)
    ntpad->Fill(layernum, phi, phi_integral / weight, z_gem, z_integral / weight);

  if (Verbosity() > 100)
    if (layernum == print_layer)
    {
      cout << " hit " << hit << " quick centroid for this electron " << endl;
      cout << "      phi centroid = " << phi_integral / weight << " phi in " << phi << " phi diff " << phi_integral / weight - phi << endl;
      cout << "      z centroid = " << z_integral / weight << " z in " << z_gem << " z diff " << z_integral / weight - z_gem << endl;
      // For a single track event, this captures the distribution of single electron centroids on the pad plane for layer print_layer.
      // The centroid of that should match the cluster centroid found by PHG4TpcClusterizer for layer print_layer, if everything is working
      //   - matches to < .01 cm for a few cases that I checked
      nthit->Fill(hit, layernum, phi, phi_integral / weight, z_gem, z_integral / weight, weight);
    }

  hit++;

  return;
}

void PHG4TpcPadPlaneReadout::populate_rectangular_phibins(const unsigned int layernum, const double phi, const double cloud_sig_rp, std::vector<int> &pad_phibin, std::vector<double> &pad_phibin_share)
{
  double cloud_sig_rp_inv = 1. / cloud_sig_rp;

  const int phibin = LayerGeom->get_phibin(phi);
  const int nphibins = LayerGeom->get_phibins();

  double radius = LayerGeom->get_radius();
  double phidisp = phi - LayerGeom->get_phicenter(phibin);
  double phistepsize = LayerGeom->get_phistep();

  // bin the charge in phi - consider phi bins up and down 3 sigma in r-phi
  int n_rp = int(3 * cloud_sig_rp / (radius * phistepsize) + 1);
  for (int iphi = -n_rp; iphi != n_rp + 1; ++iphi)
  {
    int cur_phi_bin = phibin + iphi;
    // correcting for continuity in phi
    if (cur_phi_bin < 0)
      cur_phi_bin += nphibins;
    else if (cur_phi_bin >= nphibins)
      cur_phi_bin -= nphibins;
    if ((cur_phi_bin < 0) || (cur_phi_bin >= nphibins))
    {
      std::cout << "PHG4CylinderCellTpcReco => error in phi continuity. Skipping" << std::endl;
      continue;
    }
    // Get the integral of the charge probability distribution in phi inside the current phi step
    double phiLim1 = 0.5 * M_SQRT2 * ((iphi + 0.5) * phistepsize * radius - phidisp * radius) * cloud_sig_rp_inv;
    double phiLim2 = 0.5 * M_SQRT2 * ((iphi - 0.5) * phistepsize * radius - phidisp * radius) * cloud_sig_rp_inv;
    double phi_integral = 0.5 * (erf(phiLim1) - erf(phiLim2));

    pad_phibin.push_back(cur_phi_bin);
    pad_phibin_share.push_back(phi_integral);
  }

  return;
}

void PHG4TpcPadPlaneReadout::populate_zigzag_phibins(const unsigned int layernum, const double phi, const double cloud_sig_rp, std::vector<int> &pad_phibin, std::vector<double> &pad_phibin_share)
{
  const double radius = LayerGeom->get_radius();
  const double phistepsize = LayerGeom->get_phistep();
  const auto phibins = LayerGeom->get_phibins();

  // make the charge distribution gaussian
  double rphi = phi * radius;
  if (Verbosity() > 100)
    if (LayerGeom->get_layer() == print_layer)
    {
      cout << " populate_zigzag_phibins for layer " << layernum << " with radius " << radius << " phi " << phi
           << " rphi " << rphi << " phistepsize " << phistepsize << endl;
      cout << " fcharge created: radius " << radius << " rphi " << rphi << " cloud_sig_rp " << cloud_sig_rp << endl;
    }

  // Get the range of phi values that completely contains all pads  that touch the charge distribution - (nsigmas + 1/2 pad width) in each direction
  const double philim_low = phi - (_nsigmas * cloud_sig_rp / radius) - phistepsize;
  const double philim_high = phi + (_nsigmas * cloud_sig_rp / radius) + phistepsize;

  // Find the pad range that covers this phi range
  int phibin_low = LayerGeom->get_phibin(philim_low);
  int phibin_high = LayerGeom->get_phibin(philim_high);
  int npads = phibin_high - phibin_low;

  if (Verbosity() > 1000)
    if (layernum == print_layer)
      cout << "           zigzags: phi " << phi << " philim_low " << philim_low << " phibin_low " << phibin_low
           << " philim_high " << philim_high << " phibin_high " << phibin_high << " npads " << npads << endl;

  if (npads < 0 || npads > 9) npads = 9;  // can happen if phibin_high wraps around. If so, limit to 10 pads and fix below

  // Calculate the maximum extent in r-phi of pads in this layer. Pads are assumed to touch the center of the next phi bin on both sides.
  const double pad_rphi = 2.0 * LayerGeom->get_phistep() * radius;

  // Make a TF1 for each pad in the phi range
  using PadParameterSet = std::array<double,2>;
  std::array<PadParameterSet,10> pad_parameters;
  std::array<int,10> pad_keep;
  for (int ipad = 0; ipad <= npads; ipad++)
  {
    int pad_now = phibin_low + ipad;
    // check that we do not exceed the maximum number of pads, wrap if necessary
    if (pad_now >= phibins) pad_now -= phibins;

    pad_keep[ipad] = pad_now;
    const double rphi_pad_now = LayerGeom->get_phicenter(pad_now) * radius;
    pad_parameters[ipad] = {{ pad_rphi / 2.0, rphi_pad_now }};

    if (Verbosity() > 1000)
      if (layernum == print_layer)
        cout << " zigzags: make fpad for ipad " << ipad << " pad_now " << pad_now << " pad_rphi/2 " << pad_rphi / 2.0
             << " rphi_pad_now " << rphi_pad_now << endl;
  }

  // Now make a loop that steps through the charge distribution and evaluates the response at that point on each pad

  double overlap[10];
  for (int i = 0; i < 10; i++)
    overlap[i] = 0;

  double xstep = 2.0 * _nsigmas * cloud_sig_rp / (double) _ngauss_steps;
  for (int i = 0; i < _ngauss_steps; i++)
  {
    const double x = rphi - 4.5 * cloud_sig_rp + (double) i * xstep;
    const double charge = _gauss_weights[i];
    for (int ipad = 0; ipad <= npads; ipad++)
    {
      const double pad_response = get_pad_response( x, pad_parameters[ipad] );
      if( pad_response > 0 )
      {
        const double prod = charge * pad_response;
        overlap[ipad] += prod;
      }
    }  // pads

  }  // steps

  // now we have the overlap for each pad
  for (int ipad = 0; ipad <= npads; ipad++)
  {
    pad_phibin.push_back(pad_keep[ipad]);
    pad_phibin_share.push_back(overlap[ipad]);
    if (rad_gem < output_radius) cout << "         zigzags: for pad " << ipad << " integral is " << overlap[ipad] << endl;
  }

  return;
}

void PHG4TpcPadPlaneReadout::populate_zbins(const double z, const std::array<double,2>& cloud_sig_zz, std::vector<int> &adc_zbin, std::vector<double> &adc_zbin_share)
{
  int zbin = LayerGeom->get_zbin(z);
  if (zbin < 0 || zbin > LayerGeom->get_zbins())
  {
    //cout << " z bin is outside range, return" << endl;
    return;
  }

  double zstepsize = LayerGeom->get_zstep();
  double zdisp = z - LayerGeom->get_zcenter(zbin);

  if (Verbosity() > 1000)
    cout << "     input:  z " << z << " zbin " << zbin << " zstepsize " << zstepsize << " z center " << LayerGeom->get_zcenter(zbin) << " zdisp " << zdisp << endl;
  
  // Because of diffusion, hits can be shared across the membrane, so we allow all z bins
  int min_cell_zbin = 0;
  int max_cell_zbin = NZBins - 1;

  double cloud_sig_zz_inv[2];
  cloud_sig_zz_inv[0] = 1. / cloud_sig_zz[0];
  cloud_sig_zz_inv[1] = 1. / cloud_sig_zz[1];

  int zsect = 0;
  if (z < 0)
    zsect = -1;
  else
    zsect = 1;

  int n_zz = int(3 * (cloud_sig_zz[0] + cloud_sig_zz[1]) / (2.0 * zstepsize) + 1);
  if (Verbosity() > 1000) cout << " n_zz " << n_zz << " cloud_sigzz[0] " << cloud_sig_zz[0] << " cloud_sig_zz[1] " << cloud_sig_zz[1] << endl;
  for (int iz = -n_zz; iz != n_zz + 1; ++iz)
  {
    int cur_z_bin = zbin + iz;
    if ((cur_z_bin < min_cell_zbin) || (cur_z_bin > max_cell_zbin)) continue;

    if (Verbosity() > 1000)
      cout << " iz " << iz << " cur_z_bin " << cur_z_bin << " min_cell_zbin " << min_cell_zbin << " max_cell_zbin " << max_cell_zbin << endl;

    double z_integral = 0.0;
    if (iz == 0)
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

      double zLim1 = 0.0;
      double zLim2 = 0.5 * M_SQRT2 * (-0.5 * zstepsize - zdisp) * cloud_sig_zz_inv[index1];
      // 1/2 * the erf is the integral probability from the argument Z value to zero, so this is the integral probability between the Z limits
      double z_integral1 = 0.5 * (erf(zLim1) - erf(zLim2));

      if (Verbosity() > 1000)
        if (LayerGeom->get_layer() == print_layer)
          cout << "   populate_zbins:  cur_z_bin " << cur_z_bin << "  center z " << LayerGeom->get_zcenter(cur_z_bin)
               << " index1 " << index1 << "  zLim1 " << zLim1 << " zLim2 " << zLim2 << " z_integral1 " << z_integral1 << endl;

      zLim2 = 0.0;
      zLim1 = 0.5 * M_SQRT2 * (0.5 * zstepsize - zdisp) * cloud_sig_zz_inv[index2];
      double z_integral2 = 0.5 * (erf(zLim1) - erf(zLim2));

      if (Verbosity() > 1000)
        if (LayerGeom->get_layer() == print_layer)
          cout << "   populate_zbins:  cur_z_bin " << cur_z_bin << "  center z " << LayerGeom->get_zcenter(cur_z_bin)
               << " index2 " << index2 << "  zLim1 " << zLim1 << " zLim2 " << zLim2 << " z_integral2 " << z_integral2 << endl;

      z_integral = z_integral1 + z_integral2;
    }
    else
    {
      // The non zero bins are entirely in the lead or tail region
      // lead or tail depends on which side of the membrane
      int index = 0;
      if (iz < 0)
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
      double zLim1 = 0.5 * M_SQRT2 * ((iz + 0.5) * zstepsize - zdisp) * cloud_sig_zz_inv[index];
      double zLim2 = 0.5 * M_SQRT2 * ((iz - 0.5) * zstepsize - zdisp) * cloud_sig_zz_inv[index];
      z_integral = 0.5 * (erf(zLim1) - erf(zLim2));

      if (Verbosity() > 1000)
        if (LayerGeom->get_layer() == print_layer)
          cout << "   populate_zbins:  z_bin " << cur_z_bin << "  center z " << LayerGeom->get_zcenter(cur_z_bin)
               << " index " << index << "  zLim1 " << zLim1 << " zLim2 " << zLim2 << " z_integral " << z_integral << endl;
    }

    adc_zbin.push_back(cur_z_bin);
    adc_zbin_share.push_back(z_integral);
  }

  return;
}

void PHG4TpcPadPlaneReadout::SetDefaultParameters()
{
  set_default_int_param("ntpc_layers_inner", 16);
  set_default_int_param("ntpc_layers_mid", 16);
  set_default_int_param("ntpc_layers_outer", 16);

  set_default_int_param("tpc_minlayer_inner", 7);

  set_default_double_param("tpc_minradius_inner", 30.0);  // cm
  set_default_double_param("tpc_minradius_mid", 40.0);
  set_default_double_param("tpc_minradius_outer", 60.0);

  set_default_double_param("tpc_maxradius_inner", 40.0);  // cm
  set_default_double_param("tpc_maxradius_mid", 60.0);
  set_default_double_param("tpc_maxradius_outer", 77.0);  // from Tom

  set_default_double_param("neffelectrons_threshold", 1.0);
  set_default_double_param("maxdriftlength", 105.5);         // cm
  set_default_double_param("drift_velocity", 8.0 / 1000.0);  // cm/ns
  set_default_double_param("tpc_adc_clock", 53.0);           // ns, for 18.8 MHz clock

  set_default_double_param("gem_cloud_sigma", 0.04);     // cm = 400 microns
  set_default_double_param("sampa_shaping_lead", 32.0);  // ns, for 80 ns SAMPA
  set_default_double_param("sampa_shaping_tail", 48.0);  // ns, for 80 ns SAMPA

  set_default_int_param("ntpc_phibins_inner", 1152);
  set_default_int_param("ntpc_phibins_mid", 1536);
  set_default_int_param("ntpc_phibins_outer", 2304);

  set_default_int_param("zigzag_pads", 1);

  set_default_double_param("gem_amplification", 2000); // GEM Gain

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

  Thickness[0] = NTpcLayers[0]<=0 ? 0 : (MaxRadius[0] - MinRadius[0]) / NTpcLayers[0];
  Thickness[1] = NTpcLayers[1]<=0 ? 0 : (MaxRadius[1] - MinRadius[1]) / NTpcLayers[1];
  Thickness[2] = NTpcLayers[2]<=0 ? 0 : (MaxRadius[2] - MinRadius[2]) / NTpcLayers[2];

  MaxRadius[0] = get_double_param("tpc_maxradius_inner");
  MaxRadius[1] = get_double_param("tpc_maxradius_mid");
  MaxRadius[2] = get_double_param("tpc_maxradius_outer");

  sigmaT = get_double_param("gem_cloud_sigma");
  sigmaL[0] = get_double_param("sampa_shaping_lead") * get_double_param("drift_velocity");
  sigmaL[1] = get_double_param("sampa_shaping_tail") * get_double_param("drift_velocity");

  tpc_drift_velocity = get_double_param("drift_velocity");
  tpc_adc_clock = get_double_param("tpc_adc_clock");
  ZBinWidth = tpc_adc_clock * tpc_drift_velocity;
  MaxZ = get_double_param("maxdriftlength");
  MinZ = -MaxZ;
  NZBins = (int) ((MaxZ - MinZ) / ZBinWidth) + 1;

  NPhiBins[0] = get_int_param("ntpc_phibins_inner");
  NPhiBins[1] = get_int_param("ntpc_phibins_mid");
  NPhiBins[2] = get_int_param("ntpc_phibins_outer");

  PhiBinWidth[0] = 2.0 * M_PI / (double) NPhiBins[0];
  PhiBinWidth[1] = 2.0 * M_PI / (double) NPhiBins[1];
  PhiBinWidth[2] = 2.0 * M_PI / (double) NPhiBins[2];

  zigzag_pads = get_int_param("zigzag_pads");

  averageGEMGain = get_double_param("gem_amplification");
}
