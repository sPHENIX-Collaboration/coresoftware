#include "PHG4TpcDigitizer.h"

//#include <trackbase/TrkrHit.h>
#include <trackbase/TrkrHitv2.h>
#include <trackbase/TrkrHitSet.h>
#include <trackbase/TrkrHitSetContainer.h>
#include <trackbase/TrkrHitTruthAssoc.h>
#include <trackbase/TrkrDefs.h>

#include <tpc/TpcDefs.h>

#include <g4detectors/PHG4CylinderCellGeom.h>
#include <g4detectors/PHG4CylinderCellGeomContainer.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>                         // for SubsysReco
#
#include <phool/PHCompositeNode.h>
#include <phool/PHNode.h>                               // for PHNode
#include <phool/PHNodeIterator.h>
#include <phool/PHRandomSeed.h>
#include <phool/getClass.h>
#include <phool/phool.h>                                // for PHWHERE

#include <gsl/gsl_rng.h>                                // for gsl_rng_alloc
#include <gsl/gsl_randist.h>

#include <cstdlib>                                     // for exit
#include <iostream>
#include <limits>
#include <memory>                                       // for allocator_tra...

using namespace std;

PHG4TpcDigitizer::PHG4TpcDigitizer(const string &name)
  : SubsysReco(name)
  , TpcMinLayer(7)
  , ADCThreshold(2700)
  ,  // electrons
  TpcEnc(670)
  ,  // electrons
  Pedestal(50000)
  ,  // electrons
  ChargeToPeakVolts(20)
  ,  // mV/fC
  ADCSignalConversionGain(numeric_limits<float>::signaling_NaN())
  ,  // will be assigned in PHG4TpcDigitizer::InitRun
  ADCNoiseConversionGain(numeric_limits<float>::signaling_NaN())  // will be assigned in PHG4TpcDigitizer::InitRun
{
  unsigned int seed = PHRandomSeed();  // fixed seed is handled in this funtcion
  cout << Name() << " random seed: " << seed << endl;
  RandomGenerator = gsl_rng_alloc(gsl_rng_mt19937);
  gsl_rng_set(RandomGenerator, seed);

  if (Verbosity() > 0)
    cout << "Creating PHG4TpcDigitizer with name = " << name << endl;
}

PHG4TpcDigitizer::~PHG4TpcDigitizer()
{
  gsl_rng_free(RandomGenerator);
}

int PHG4TpcDigitizer::InitRun(PHCompositeNode *topNode)
{
  ADCThreshold += Pedestal;

  // Factor that convertes the charge in each z bin into a voltage in each z bin
  // ChargeToPeakVolts relates TOTAL charge collected to peak voltage, while the cell maker actually distributes the signal
  // GEM output charge in Z bins using the shaper time response. For 80 ns shaping, the scaleup factor of 2.4 gets the peak voltage right.
  ADCSignalConversionGain = ChargeToPeakVolts * 1.60e-04 * 2.4;  // 20 (or 30) mV/fC * fC/electron * scaleup factor
  // The noise is by definition the RMS noise width voltage divided by ChargeToPeakVolts
  ADCNoiseConversionGain = ChargeToPeakVolts * 1.60e-04;  // 20 (or 30) mV/fC * fC/electron

  //-------------
  // Add Hit Node
  //-------------
  PHNodeIterator iter(topNode);

  // Looking for the DST node
  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    cout << PHWHERE << "DST Node missing, doing nothing." << endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  PHNodeIterator iter_dst(dstNode);

  CalculateCylinderCellADCScale(topNode);

  //----------------
  // Report Settings
  //----------------

  if (Verbosity() > 0)
  {
    cout << "====================== PHG4TpcDigitizer::InitRun() =====================" << endl;
    for (std::map<int, unsigned int>::iterator iter = _max_adc.begin();
         iter != _max_adc.end();
         ++iter)
    {
      cout << " Max ADC in Layer #" << iter->first << " = " << iter->second << endl;
    }
    for (std::map<int, float>::iterator iter = _energy_scale.begin();
         iter != _energy_scale.end();
         ++iter)
    {
      cout << " Energy per ADC in Layer #" << iter->first << " = " << 1.0e6 * iter->second << " keV" << endl;
    }
    cout << "===========================================================================" << endl;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int PHG4TpcDigitizer::process_event(PHCompositeNode *topNode)
{

  DigitizeCylinderCells(topNode);

  return Fun4AllReturnCodes::EVENT_OK;
}

void PHG4TpcDigitizer::CalculateCylinderCellADCScale(PHCompositeNode *topNode)
{
  // defaults to 8-bit ADC, short-axis MIP placed at 1/4 dynamic range

  PHG4CylinderCellGeomContainer *geom_container = findNode::getClass<PHG4CylinderCellGeomContainer>(topNode, "CYLINDERCELLGEOM_SVTX");

  if (!geom_container) return;

  PHG4CylinderCellGeomContainer::ConstRange layerrange = geom_container->get_begin_end();
  for (PHG4CylinderCellGeomContainer::ConstIterator layeriter = layerrange.first;
       layeriter != layerrange.second;
       ++layeriter)
  {
    int layer = layeriter->second->get_layer();
    float thickness = (layeriter->second)->get_thickness();
    float pitch = (layeriter->second)->get_phistep() * (layeriter->second)->get_radius();
    float length = (layeriter->second)->get_zstep();

    float minpath = pitch;
    if (length < minpath) minpath = length;
    if (thickness < minpath) minpath = thickness;
    float mip_e = 0.003876 * minpath;

    if (_max_adc.find(layer) == _max_adc.end())
    {
      _max_adc[layer] = 255;
      _energy_scale[layer] = mip_e / 64;
    }
  }

  return;
}

void PHG4TpcDigitizer::DigitizeCylinderCells(PHCompositeNode *topNode)
{
  unsigned int print_layer = 18;  // to print diagnostic output for layer 47

  // Digitizes the Tpc cells that were created in PHG4CylinderCellTpcReco
  // These contain as edep the number of electrons out of the GEM stack, distributed between Z bins by shaper response and ADC clock window
  // - i.e. all of the phi and Z bins in a cluster have edep values that add up to the primary charge in the layer times 2000

  // NOTES:
  // Modified by ADF June 2018 to do the following:
  //   Add noise to cells before digitizing
  //   Digitize the first adc time bin to exceed the threshold, and the 4 bins after that
  //   If the adc value is still above the threshold after 5 bins, repeat for the next 5 bins

  // Electron production:
  // A MIP produces 32 electrons in 1.25 cm of Ne:CF4 gas
  // The nominal GEM gain is 2000 => 64,000 electrons per MIP through 1.25 cm gas
  // Thus a MIP produces a charge value out of the GEM stack of 64000/6.242x10^18 = 10.2 fC

  // SAMPA:
  // See https://indico.cern.ch/event/489996/timetable/#all.detailed "SAMPA Chip: the New ASIC for the ALICE Tpc and MCH Upgrades", M Bregant
  // The SAMPA has a maximum output voltage of 2200 mV (but the pedestal is about 200 mV)
  // The SAMPA shaper is set to 80 ns peaking time
  // The ADC Digitizes the SAMPA shaper output into 1024 channels
  // Conversion gains of 20 mV/fC or 30 mV/fC are possible - 1 fC charge input produces a peak volatge out of the shaper of 20 or 30 mV
  //   At 30 mV/fC, the input signal saturates at 2.2 V / 30 mV/fC = 73 fC (say 67 with pedestal not at zero)
  //   At 20 mV/fC, the input signal saturates at 2.2 V / 20 mV/fC = 110 fC (say 100 fC with pedestal not at zero) - assume 20 mV/fC
  // The equivalent noise charge RMS at 20 mV/fC was measured (w/o detector capacitance) at 490 electrons
  //      - note: this appears to be just the pedestal RMS voltage spread divided by the conversion gain, so it is a bit of a funny number (see below)
  //      - it is better to think of noise and signal in terms of voltage at the input of the ADC
  // Bregant's slides say 670 electrons ENC for the full chip with 18 pf detector, as in ALICE - should use that number

  // Signal:
  // To normalize the signal in each cell, take the entire charge on the pad and multiply by 20 mV/fC to get the adc input AT THE PEAK of the shaper
  // The cell contents should thus be multipied by the normalization given by:
  // V_peak = Q_pad (electrons) * 1.6e-04 fC/electron * 20 mV/fC
  // From the sims, for 80 ns and 18.8 MHz, if we take the input charge and spread it a across the shaping time (which is how it has always been done, and is
  // not the best way to think about it, because the ADC does not see charge it sees voltage out of a charge integrating preamp followed by a shaper), to get
  // the voltage at the ADC input right, then the values of Q_(pad,z) have to be scaled up by 2.4
  // V_(pad,z) = 2.4 * Q_(pad,z) (electrons) * 1.6e-04 fC/electron * 20 mV/fC = Q_(pad,z) * 7.68e-03 (in mV)
  // ADU_(pad,z) = V_(pad,z) * (1024 ADU / 2200 mV) = V_(pad,z) * 0.465
  // Remember that Q_(pad,z) is the GEM output charge

  // Noise:
  // The ENC is defined as the RMS spread of the ADC pedestal distribution coming out from the SAMPA divided by the corresponding conversion gain.
  // The full range of the ADC input is 2.2V (which will become 1024 adc counts, i.e. 1024 ADU's).
  // If you see the RMS of the pedestal in adc counts as 1 at the gain of 20mV/fC, the ENC would be defined by:
  //             (2200 [mV]) * (1/1024) / (20[mV/fC]) / (1.6*10^-4 [fC]) = 671 [electrons]
  // The RMS noise voltage would be: V_RMS = ENC (electrons) *1.6e-04 fC/electron * 20 mV/fC = ENC (electrons) * 3.2e-03 (in mV)
  // The ADC readout would be:  ADU = V_RMS * (1024 ADU / 2200 mV) = V_RMS * 0.465

  // The cells that we need to digitize here contain as the energy "edep", which is the number of electrons out of the GEM stack
  // distributed over the pads and ADC time bins according to the output time distribution of the SAMPA shaper - not what really happens, see above
  // We convert to volts at the input to the ADC and add noise generated with the RMS value of the noise voltage at the ADC input
  // We assume the pedestal is zero, for simplicity, so the noise fluctuates above and below zero

  // Note that zbin = 0 corresponds to -100.5 cm, zbin 248 corresponds to 0 cm, and zbin 497 corresponds to +100.5 cm
  // increasing time should be bins (497 -> 249) and (0 -> 248)

  //----------
  // Get Nodes
  //----------

  PHG4CylinderCellGeomContainer *geom_container =
    findNode::getClass<PHG4CylinderCellGeomContainer>(topNode, "CYLINDERCELLGEOM_SVTX");
  if (!geom_container)
    {
      std::cout << PHWHERE << "ERROR: Can't find node CYLINDERCELLGEOM_SVTX" << std::endl;
    }


  // Get the TrkrHitSetContainer node
  TrkrHitSetContainer *trkrhitsetcontainer = findNode::getClass<TrkrHitSetContainer>(topNode, "TRKR_HITSET");
  if(!trkrhitsetcontainer)
    {
      cout << "Could not locate TRKR_HITSET node, quit! " << endl;
      exit(1);
    }

  TrkrHitTruthAssoc *hittruthassoc = findNode::getClass<TrkrHitTruthAssoc>(topNode, "TRKR_HITTRUTHASSOC");
  if(!hittruthassoc)
    {
      cout << PHWHERE << " Could not locate TRKR_HITTRUTHASSOC node, quit! " << endl;
      exit(1);
    }


  //-------------
  // Digitization
  //-------------

  // Loop over all hitsets for the Tpc
  TrkrHitSetContainer::ConstRange hitset_range = trkrhitsetcontainer->getHitSets(TrkrDefs::TrkrId::tpcId);
  for (TrkrHitSetContainer::ConstIterator hitset_iter = hitset_range.first;
       hitset_iter != hitset_range.second;
       ++hitset_iter)
    {
      // we have an itrator to one TrkrHitSet for the Tpc from the trkrHitSetContainer
      // get the hitset key so we can find the layer
      TrkrDefs::hitsetkey hitsetkey = hitset_iter->first;
      const unsigned int layer = TrkrDefs::getLayer(hitsetkey);

      if(Verbosity() > 2)
	if(layer == print_layer) 
	  cout << "new: PHG4TpcDigitizer:  processing hits for layer " << layer << " hitsetkey " << hitsetkey << endl;

      // we need the geometry object for this layer
      PHG4CylinderCellGeom *layergeom = geom_container->GetLayerCellGeom(layer);
      if (!layergeom)
	exit(1);

      // for this hitset, make a vector of a vector of cells for each phibin
      phi_sorted_hits.clear();
      
      // start with an empty vector of cells for each phibin      
      int nphibins = layergeom->get_phibins();
      for (int iphi = 0; iphi < nphibins; iphi++)
	{
	  phi_sorted_hits.push_back(std::vector<TrkrHitSet::ConstIterator>());
	}      
      
      // get all of the hits from this hitset      
      TrkrHitSet *hitset = hitset_iter->second;
      TrkrHitSet::ConstRange hit_range = hitset->getHits();
      for(TrkrHitSet::ConstIterator hit_iter = hit_range.first;
	  hit_iter != hit_range.second;
	  ++hit_iter)
	{
	  // Fill the vector of hits for each phibin
	  unsigned int phibin = TpcDefs::getPad(hit_iter->first);
	  phi_sorted_hits[phibin].push_back(hit_iter);
	}
	  	  
      // For this hitset we have the hits sorted into vectors for each phi
      // process these vectors one phi bin at a time
      for (unsigned int iphi = 0; iphi < phi_sorted_hits.size(); iphi++)
	{
	  if (phi_sorted_hits[iphi].size() == 0)
	    continue;
	  
	  // Populate a vector of hits ordered by Z for each phibin
	  int nzbins = layergeom->get_zbins();
	  is_populated.clear();
	  is_populated.assign(nzbins, 2);  // mark all as noise only for now
	  z_sorted_hits.clear();
	      
	  // add an empty vector for each z bin
	  for (int iz = 0; iz < nzbins; iz++)
	    z_sorted_hits.push_back(std::vector<TrkrHitSet::ConstIterator>());

	  if(Verbosity() > 2)
	    if(layer == print_layer) cout << endl;	  

	  // add a hit for each z bin that has one
	  for (unsigned int iz = 0; iz < phi_sorted_hits[iphi].size(); iz++)
	    {
	      int zbin = TpcDefs::getTBin(phi_sorted_hits[iphi][iz]->first);
	      is_populated[zbin] = 1;  // this bin is a associated with a hit
	      z_sorted_hits[zbin].push_back(phi_sorted_hits[iphi][iz]);

	      if(Verbosity() > 2)
		if(layer == print_layer)  
		  {
		    TrkrDefs::hitkey hitkey =  phi_sorted_hits[iphi][iz]->first ;
		    cout << "iphi " << iphi << " adding existing hit to z vector for layer " << layer << " zbin " << zbin << "  hitkey " << hitkey << " pad " << TpcDefs::getPad(hitkey) 
			 << " z bin " << TpcDefs::getTBin(hitkey)   << "  energy " <<( phi_sorted_hits[iphi][iz]->second)->getEnergy() 
			 << endl;
		  }
	    }
	      
	  adc_input.clear();
	  adc_hitid.clear();
	  // Now for this phibin we process all bins ordered by Z into hits with noise
	  //======================================================
	  // For this step we take the edep value and convert it to mV at the ADC input
	  // See comments above for how to do this for signal and noise
	      
	  for (int iz = 0; iz < nzbins; iz++)
	    {
	      if (is_populated[iz] == 1)
		{
		  // This zbin has a hit, add noise
		  float noise = added_noise();                                                                                       // in electrons
		  float noise_voltage = (Pedestal + noise) * ADCNoiseConversionGain;                      // mV - from definition of noise charge and pedestal charge
		  float adc_input_voltage = (z_sorted_hits[iz][0]->second)->getEnergy() * ADCSignalConversionGain;  // mV, see comments above
		  
		  adc_input.push_back(adc_input_voltage + noise_voltage);
		  adc_hitid.push_back(z_sorted_hits[iz][0]->first);

		  if(Verbosity() > 2)
		  if(layer == print_layer) 
		    cout << "existing hit: layer " << layer << " iphi " << iphi  << " iz " << iz << " edep " <<  (z_sorted_hits[iz][0]->second)->getEnergy() 
			   << " adc gain " << ADCSignalConversionGain << " adc_input_voltage " << adc_input_voltage << " noise voltage " << noise_voltage 
			   <<  " adc_input " << adc_input[iz] << endl;
		}
	      else if (is_populated[iz] == 2)
		{
		  // This z bin does not have a filled cell, add noise
		  float noise = added_noise();                                        // returns "electrons"
		  float noise_voltage = (Pedestal + noise) * ADCNoiseConversionGain;  // mV - noise "electrons" x conversion gain
		  
		  adc_input.push_back(noise_voltage);  // mV
		  adc_hitid.push_back(0);             // there is no hit, just add a placeholder in the vector for now, replace it later

		  if(Verbosity() > 2)
		  if(layer == print_layer) 
		    cout << "noise hit: layer " << layer << " iphi " << iphi  << " iz " << iz
			   << " adc gain " << ADCSignalConversionGain << " noise voltage " << noise_voltage 
			   <<  " adc_input " << adc_input[iz] << endl;
		}
	      else
		{
		  // Cannot happen
		  cout << "Impossible value of is_populated, iz = " << iz << " is_populated = " << is_populated[iz] << endl;
		  exit(-1);
		}
	    }
	  
	  // Now we can digitize the entire stream of z bins for this phi bin
	      
	  // start with negative z, the first to arrive is bin 0
	  //================================
	  int binpointer = 0;
	  for (int iz = 0; iz < nzbins / 2; iz++)  // 0-247
	    {
	      if (iz < binpointer) continue;

	      if (adc_input[iz] > ADCThreshold * ADCNoiseConversionGain)  // convert threshold in "equivalent electrons" to mV
		{
		  // digitize this bin and the following 4 bins
		  
		  if (Verbosity() > 2)
		    if (layer == print_layer) cout << endl
						   << "new:   (neg z) Above threshold of " << ADCThreshold * ADCNoiseConversionGain << " for phibin " << iphi
						   << " iz " << iz << " with adc_input " << adc_input[iz] << " digitize this and 4 following bins: " << endl;
		  
		  for (int izup = 0; izup < 5; izup++)
		    {
		      if (iz + izup < nzbins / 2 && iz + izup >= 0)   // stay within the bin limits for negative z
			{
			  unsigned int adc_output = (unsigned int) (adc_input[iz + izup] * 1024.0 / 2200.0);  // input voltage x 1024 channels over 2200 mV max range
			  if (adc_input[iz + izup] < 0) adc_output = 0;
			  if (adc_output > 1023) adc_output = 1023;
			  
			  if (Verbosity() > 2)
			  if (layer == print_layer) cout << "new:  iphi " << iphi << "  (neg z) iz+izup " << iz + izup << " adc_hitid " << adc_hitid[iz + izup]
							   << "  adc_input " << adc_input[iz + izup] << " ADCThreshold " << ADCThreshold * ADCNoiseConversionGain
							   << " adc_output " << adc_output << endl;
			  
			  // Get the hitkey
			  TrkrDefs::hitkey hitkey = TpcDefs::genHitKey(iphi, iz+izup);
			  TrkrHit *hit = nullptr;
			  hit = hitset_iter->second->getHit(hitkey);
			  
			  // noise bins do not have TrkrHits associated with them, have to make one
			  if(!hit)
			    {
			      hit = new TrkrHitv2();
			      hitset_iter->second->addHitSpecificKey(hitkey, hit);
			      //hit->addEnergy(adc_input[iz+izup]);
			      
			      if (Verbosity() > 2)
				if (layer == print_layer) cout << "new:  adding noise hit for iphi " << iphi << " zbin " << iz + izup
							       << " created new hit with hitkey " << hitkey
							       << " energy " << adc_input[iz + izup] << " adc " << adc_output << endl;
			    }

			  hit->setAdc(adc_output);
			  
			}  // end boundary check
		      binpointer++;   // skip this bin in future
		    }    // end izup loop
		  
		}  //  adc threshold if
	      else
		{
		  // set adc value to zero if there is a hit
		  // Get the hitkey
		  TrkrDefs::hitkey hitkey = TpcDefs::genHitKey(iphi, iz);
		  TrkrHit *hit = nullptr;
		  hit = hitset_iter->second->getHit(hitkey);
		  if(hit)
		    {
		      hit->setAdc(0);
		    }
		  // below threshold, move on
		  binpointer++;
		}  // end adc threshold if/else
	      
	    }  // end iz loop for negative z
	  
	  // now positive z, the first to arrive is bin 497
	  //===============================
	  binpointer = nzbins - 1;
	  for (int iz = nzbins - 1; iz >= nzbins / 2; iz--)  // 495 - 248
	    {
	      if (iz > binpointer) continue;
	      
	      if (adc_input[iz] > ADCThreshold * ADCNoiseConversionGain)  // convert threshold in electrons to mV
		{
		  // digitize this bin and the following 4 bins
		  
		  if (Verbosity() > 2)
		    if (layer == print_layer) cout << endl
						   << "new:  (pos z) Above threshold  of " << ADCThreshold * ADCNoiseConversionGain << " for iphi " << iphi << "  iz " << iz
						   << " with adc_input " << adc_input[iz] << " digitize this and 4 following bins: " << endl;
		  
		  for (int izup = 0; izup < 5; izup++)
		    {
		      if (iz - izup < nzbins && iz - izup >= nzbins / 2)
			{
			  unsigned int adc_output = (unsigned int) (adc_input[iz - izup] * 1024.0 / 2200.0);  // input voltage x 1024 channels over 2200 mV max range
			  if (adc_input[iz - izup] < 0) adc_output = 0;
			  if (adc_output > 1023) adc_output = 1023;
			  
			  if (Verbosity() > 2)
			  if (layer == print_layer) cout << "new:  iphi " << iphi << "  (pos z) iz-izup " << iz - izup << " adc_hitid " << adc_hitid[iz - izup]
							   << "  adc_input " << adc_input[iz - izup] << " ADCThreshold " << ADCThreshold * ADCNoiseConversionGain
							   << " adc_output " << adc_output << endl;

			  // Get the hitkey
			  TrkrDefs::hitkey hitkey = TpcDefs::genHitKey(iphi, iz-izup);
			  TrkrHit *hit = nullptr;
			  hit = hitset_iter->second->getHit(hitkey);
			  
			  // noise bins do not have TrkrHits associated with them, have to make one
			  if(!hit)
			    {
			      hit = new TrkrHitv2();
			      hitset_iter->second->addHitSpecificKey(hitkey, hit);
			      //hit->addEnergy(adc_input[iz-izup]);
			      
			      if (Verbosity() > 2)
				if (layer == print_layer) cout << "new:  adding noise hit for iphi " << iphi << " zbin " << iz - izup
							       << " created new hit with hitkey " << hitkey
							       << " energy " << adc_input[iz - izup] << " adc " << adc_output << endl;
			    }

			  hit->setAdc(adc_output);
			} // end boundary check
		      binpointer--;
		    } // end izup loop
		  
		}  //  adc threshold if
	      else
		{
		  // set adc value to zero if there is a hit
		  // Get the hitkey
		  TrkrDefs::hitkey hitkey = TpcDefs::genHitKey(iphi, iz);
		  TrkrHit *hit = nullptr;
		  hit = hitset_iter->second->getHit(hitkey);
		  if(hit)
		    {
		      hit->setAdc(0);
		    }
		  // below threshold, move on
		  binpointer--;
		}  // end adc threshold if/else
	      
	    }  // end iz loop for positive z
			      
	}  // end phibins loop
      
    }  // end loop over hitsets

  //======================================================  
  if(Verbosity() > 2) 
  {
    cout << "From PHG4TpcDigitizer: hitsetcontainer dump at end before cleaning:" << endl;
  }
  std::vector<std::pair<TrkrDefs::hitsetkey, TrkrDefs::hitkey>> delete_hitkey_list;

  // Clean up undigitized hits - we want all hitsets for the Tpc
  TrkrHitSetContainer::ConstRange hitset_range_now = trkrhitsetcontainer->getHitSets(TrkrDefs::TrkrId::tpcId);
  for (TrkrHitSetContainer::ConstIterator hitset_iter = hitset_range_now.first;
       hitset_iter != hitset_range_now.second;
       ++hitset_iter)
    {
     // we have an iterator to one TrkrHitSet for the Tpc from the trkrHitSetContainer
      TrkrDefs::hitsetkey hitsetkey = hitset_iter->first;
      const unsigned int layer = TrkrDefs::getLayer(hitsetkey);
      const int sector = TpcDefs::getSectorId(hitsetkey);
      const int side = TpcDefs::getSide(hitsetkey);
      if(Verbosity() > 2) 
	cout << "PHG4TpcDigitizer: hitset with key: " << hitsetkey << " in layer " << layer << " with sector " << sector << " side " << side << endl;

      // get all of the hits from this hitset      
      TrkrHitSet *hitset = hitset_iter->second;
      TrkrHitSet::ConstRange hit_range = hitset->getHits();
      for(TrkrHitSet::ConstIterator hit_iter = hit_range.first;
	  hit_iter != hit_range.second;
	  ++hit_iter)
	{
	  TrkrDefs::hitkey hitkey = hit_iter->first;
	  TrkrHit *tpchit = hit_iter->second;

	  if(Verbosity() > 2)
	    cout << "    layer " << layer << "  hitkey " << hitkey << " pad " << TpcDefs::getPad(hitkey) << " z bin " << TpcDefs::getTBin(hitkey) 
		 << " adc " << tpchit->getAdc() << endl;

	  if(tpchit->getAdc() == 0)
	    {
	      if(Verbosity() > 20) 
		cout << "                       --   this hit not digitized - delete it" << endl;
	      // screws up the iterator to delete it here, store the hitkey for later deletion
		delete_hitkey_list.push_back(std::make_pair(hitsetkey, hitkey));
	    }
	}
    }

  // delete all undigitized hits
  for(unsigned int i=0; i<delete_hitkey_list.size(); i++)
    {
      TrkrHitSet *hitset = trkrhitsetcontainer->findHitSet(delete_hitkey_list[i].first);
      const unsigned int layer = TrkrDefs::getLayer(delete_hitkey_list[i].first);
      hitset->removeHit(delete_hitkey_list[i].second);
      if(Verbosity() > 20) 
	if (layer == print_layer)
	  cout << "removed hit with hitsetkey " << delete_hitkey_list[i].first << " and hitkey " << delete_hitkey_list[i].second << endl; 

      // should also delete all entries with this hitkey from the TrkrHitTruthAssoc map
      hittruthassoc->removeAssoc(delete_hitkey_list[i].first, delete_hitkey_list[i].second);
    }


  // Final hitset dump
  if(Verbosity() > 2) 
    cout << "From PHG4TpcDigitizer: hitsetcontainer dump at end after cleaning:" << endl;
 // We want all hitsets for the Tpc
  TrkrHitSetContainer::ConstRange hitset_range_final = trkrhitsetcontainer->getHitSets(TrkrDefs::TrkrId::tpcId);
  for (TrkrHitSetContainer::ConstIterator hitset_iter = hitset_range_final.first;
       hitset_iter != hitset_range_now.second;
       ++hitset_iter)
    {
     // we have an itrator to one TrkrHitSet for the Tpc from the trkrHitSetContainer
      TrkrDefs::hitsetkey hitsetkey = hitset_iter->first;
      const unsigned int layer = TrkrDefs::getLayer(hitsetkey);
      if (layer != print_layer)  continue;
      const int sector = TpcDefs::getSectorId(hitsetkey);
      const int side = TpcDefs::getSide(hitsetkey);
      if(Verbosity() > 2 && layer == print_layer) 
	cout << "PHG4TpcDigitizer: hitset with key: " << hitsetkey << " in layer " << layer << " with sector " << sector << " side " << side << endl;

      // get all of the hits from this hitset      
      TrkrHitSet *hitset = hitset_iter->second;
      TrkrHitSet::ConstRange hit_range = hitset->getHits();
      for(TrkrHitSet::ConstIterator hit_iter = hit_range.first;
	  hit_iter != hit_range.second;
	  ++hit_iter)
	{
	  TrkrDefs::hitkey hitkey = hit_iter->first;
	  TrkrHit *tpchit = hit_iter->second;
	  if(Verbosity() > 2)
	  cout << "      LAYER " << layer << " hitkey " << hitkey << " pad " << TpcDefs::getPad(hitkey) << " z bin " << TpcDefs::getTBin(hitkey) 
		 << " adc " << tpchit->getAdc() << endl;

	  if(tpchit->getAdc() == 0)
	    {
	      cout << "   Oops!                    --   this hit not digitized and not deleted!" << endl;
	    }
	}
    }

  //hittruthassoc->identify();

  return;
}

float PHG4TpcDigitizer::added_noise()
{
  float noise = gsl_ran_gaussian(RandomGenerator, TpcEnc);

  return noise;
}
