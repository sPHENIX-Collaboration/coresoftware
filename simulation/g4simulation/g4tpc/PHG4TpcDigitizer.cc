#include "PHG4TpcDigitizer.h"

#include <trackbase/TpcDefs.h>
#include <trackbase/TrkrDefs.h>
#include <trackbase/TrkrHit.h>
#include <trackbase/TrkrHitSet.h>
#include <trackbase/TrkrHitSetContainer.h>
#include <trackbase/TrkrHitSetTpc.h>
#include <trackbase/TrkrHitTruthAssoc.h>
#include <trackbase/TrkrHitv2.h>

#include <g4detectors/PHG4TpcCylinderGeom.h>
#include <g4detectors/PHG4TpcCylinderGeomContainer.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>  // for SubsysReco
#
#include <phool/PHCompositeNode.h>
#include <phool/PHNode.h>  // for PHNode
#include <phool/PHNodeIterator.h>
#include <phool/PHRandomSeed.h>
#include <phool/getClass.h>
#include <phool/phool.h>  // for PHWHERE

#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>  // for gsl_rng_alloc

#include <cassert>
#include <cstdlib>  // for exit
#include <iostream>
#include <limits>
#include <memory>  // for allocator_tra...

PHG4TpcDigitizer::PHG4TpcDigitizer(const std::string &name)
  : SubsysReco(name)
  , ADCThreshold(2700)                                                    // electrons
  , TpcEnc(670)                                                           // electrons
  , Pedestal(50000)                                                       // electrons
  , ChargeToPeakVolts(20)                                                 // mV/fC
  , ADCSignalConversionGain(std::numeric_limits<float>::signaling_NaN())  // will be assigned in PHG4TpcDigitizer::InitRun
  , ADCNoiseConversionGain(std::numeric_limits<float>::signaling_NaN())   // will be assigned in PHG4TpcDigitizer::InitRun
{
  unsigned int seed = PHRandomSeed();  // fixed seed is handled in this funtcion
  std::cout << Name() << " random seed: " << seed << std::endl;
  RandomGenerator = gsl_rng_alloc(gsl_rng_mt19937);
  gsl_rng_set(RandomGenerator, seed);

  if (Verbosity() > 0)
  {
    std::cout << "Creating PHG4TpcDigitizer with name = " << name << std::endl;
  }
}

PHG4TpcDigitizer::~PHG4TpcDigitizer()
{
  gsl_rng_free(RandomGenerator);
}

int PHG4TpcDigitizer::InitRun(PHCompositeNode *topNode)
{
  ADCThreshold += Pedestal;

  // Factor that converts the charge in each z bin into a voltage in each z bin
  // ChargeToPeakVolts relates TOTAL charge collected to peak voltage, while the cell maker actually distributes the signal
  // GEM output charge in Z bins using the shaper time response. For 80 ns shaping, the scaleup factor of 2.4 gets the peak voltage right.
  ADCSignalConversionGain = ChargeToPeakVolts * 1.60e-04 * 2.4;  // 20 (or 30) mV/fC * fC/electron * scaleup factor
  // The noise is by definition the RMS noise width voltage divided by ChargeToPeakVolts
  ADCNoiseConversionGain = ChargeToPeakVolts * 1.60e-04;  // 20 (or 30) mV/fC * fC/electron

  ADCThreshold_mV = ADCThreshold * ADCNoiseConversionGain;

  //-------------
  // Add Hit Node
  //-------------
  PHNodeIterator iter(topNode);

  // Looking for the DST node
  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    std::cout << PHWHERE << "DST Node missing, doing nothing." << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  PHNodeIterator iter_dst(dstNode);

  CalculateCylinderCellADCScale(topNode);

  //----------------
  // Report Settings
  //----------------

  if (Verbosity() > 0)
  {
    std::cout << "====================== PHG4TpcDigitizer::InitRun() =====================" << std::endl;
    for (std::map<int, unsigned int>::iterator tpiter = _max_adc.begin();
         tpiter != _max_adc.end();
         ++tpiter)
    {
      std::cout << " Max ADC in Layer #" << tpiter->first << " = " << tpiter->second << std::endl;
    }
    for (std::map<int, float>::iterator tpiter = _energy_scale.begin();
         tpiter != _energy_scale.end();
         ++tpiter)
    {
      std::cout << " Energy per ADC in Layer #" << tpiter->first << " = " << 1.0e6 * tpiter->second << " keV" << std::endl;
    }
    std::cout << "===========================================================================" << std::endl;
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

  PHG4TpcCylinderGeomContainer *geom_container = findNode::getClass<PHG4TpcCylinderGeomContainer>(topNode, "CYLINDERCELLGEOM_SVTX");

  if (!geom_container) return;

  PHG4TpcCylinderGeomContainer::ConstRange layerrange = geom_container->get_begin_end();
  for (PHG4TpcCylinderGeomContainer::ConstIterator layeriter = layerrange.first;
       layeriter != layerrange.second;
       ++layeriter)
  {
    int layer = layeriter->second->get_layer();
    float thickness = (layeriter->second)->get_thickness();
    float pitch = (layeriter->second)->get_phistep() * (layeriter->second)->get_radius();
    float length = (layeriter->second)->get_zstep() * _drift_velocity;

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
  // See https://indico.cern.ch/event/489996/timetable/#all.detailed
  //      "SAMPA Chip: the New ASIC for the ALICE Tpc and MCH Upgrades", M Bregant
  // The SAMPA has a maximum output voltage of 2200 mV (but the pedestal is about 200 mV)
  // The SAMPA shaper is set to 80 ns peaking time
  // The ADC Digitizes the SAMPA shaper output into 1024 channels
  // Conversion gains of 20 mV/fC or 30 mV/fC are possible - 1 fC charge input produces a peak voltage out of
  // the shaper of 20 or 30 mV
  //   At 30 mV/fC, the input signal saturates at 2.2 V / 30 mV/fC = 73 fC (say 67 with pedestal not at zero)
  //   At 20 mV/fC, the input signal saturates at 2.2 V / 20 mV/fC = 110 fC (say 100 fC with pedestal not at zero) - assume 20 mV/fC
  // The equivalent noise charge RMS at 20 mV/fC was measured (w/o detector capacitance) at 490 electrons
  //      - note: this appears to be just the pedestal RMS voltage spread divided by the conversion gain, so it is a bit of a
  //         funny number (see below)
  //      - it is better to think of noise and signal in terms of voltage at the input of the ADC
  // Bregant's slides say 670 electrons ENC for the full chip with 18 pf detector, as in ALICE - should use that number

  // Signal:
  // To normalize the signal in each cell, take the entire charge on the pad and multiply by 20 mV/fC to get the adc
  //       input AT THE PEAK of the shaper
  // The cell contents should thus be multipied by the normalization given by:
  // V_peak = Q_pad (electrons) * 1.6e-04 fC/electron * 20 mV/fC
  // From the sims, for 80 ns and 18.8 MHz, if we take the input charge and spread it a across the shaping time (which is how it has always been done, and is
  // not the best way to think about it, because the ADC does not see charge it sees voltage out of a charge integrating
  //     preamp followed by a shaper), to get
  // the voltage at the ADC input right, then the values of Q_(pad,z) have to be scaled up by 2.4
  // V_(pad,z) = 2.4 * Q_(pad,z) (electrons) * 1.6e-04 fC/electron * 20 mV/fC = Q_(pad,z) * 7.68e-03 (in mV)
  // ADU_(pad,z) = V_(pad,z) * (1024 ADU / 2200 mV) = V_(pad,z) * 0.465
  // Remember that Q_(pad,z) is the GEM output charge

  // Noise:
  // The ENC is defined as the RMS spread of the ADC pedestal distribution coming out from the SAMPA
  //      divided by the corresponding conversion gain.
  // The full range of the ADC input is 2.2V (which will become 1024 adc counts, i.e. 1024 ADU's).
  // If you see the RMS of the pedestal in adc counts as 1 at the gain of 20mV/fC, the ENC would be defined by:
  //             (2200 [mV]) * (1/1024) / (20[mV/fC]) / (1.6*10^-4 [fC]) = 671 [electrons]
  // The RMS noise voltage would be:
  //     V_RMS = ENC (electrons) *1.6e-04 fC/electron * 20 mV/fC = ENC (electrons) * 3.2e-03 (in mV)
  // The ADC readout would be:  ADU = V_RMS * (1024 ADU / 2200 mV) = V_RMS * 0.465

  // The cells that we need to digitize here contain as the energy "edep", which is the number of electrons out of the GEM stack
  // distributed over the pads and ADC time bins according to the output time distribution of the SAMPA shaper
  //      - not what really happens, see above
  // We convert to volts at the input to the ADC and add noise generated with the RMS value of the noise voltage at the ADC input
  // We assume the pedestal is zero, for simplicity, so the noise fluctuates above and below zero

  // Note that tbin = 0 corresponds to -105.5 cm, tbin 248 corresponds to 0 cm, and tbin 497 corresponds to +105.5 cm
  // increasing time should be bins (497 -> 249) and (0 -> 248)

  //----------
  // Get Nodes
  //----------

  PHG4TpcCylinderGeomContainer *geom_container =
      findNode::getClass<PHG4TpcCylinderGeomContainer>(topNode, "CYLINDERCELLGEOM_SVTX");
  if (!geom_container)
  {
    std::cout << PHWHERE << "ERROR: Can't find node CYLINDERCELLGEOM_SVTX" << std::endl;
  }

  // Get the TrkrHitSetContainer node
  TrkrHitSetContainer *trkrhitsetcontainer = findNode::getClass<TrkrHitSetContainer>(topNode, "TRKR_HITSET_TPC");
  if (!trkrhitsetcontainer)
  {
    std::cout << "Could not locate TRKR_HITSET_TPC node, quit! " << std::endl;
    exit(1);
  }
  //-------------
  // Digitization
  //-------------

  TrkrHitSetContainer::ConstRange hitset_range = trkrhitsetcontainer->getHitSets();
  for (TrkrHitSetContainer::ConstIterator hitset_iter = hitset_range.first;
       hitset_iter != hitset_range.second;
       ++hitset_iter)
  {
    // we have an iterator to one TrkrHitSet for the Tpc from the trkrHitSetContainer
    // get the hitset key
    TrkrDefs::hitsetkey hitsetkey = hitset_iter->first;

    // get all of the hits from this hitset
    TrkrHitSetTpc *hitset = dynamic_cast<TrkrHitSetTpc *>(hitset_iter->second);

    if (hitset == nullptr)
    {
      std::cout << __PRETTY_FUNCTION__ << " : fatal error : hitset received is not a TrkrHitSetTpc with key " << hitsetkey << " and content ";
      hitset_iter->second->identify();
    }
    assert(hitset);

    std::vector<bool> pass_zero_suppression;

    for (auto & t_sorted_hits : hitset->getTimeFrameAdcData())
    {
      // add noise
      for (auto & adc_bin : t_sorted_hits)
      {
        adc_bin = static_cast<TpcDefs::ADCDataType>(add_noise_to_bin((float) adc_bin));
      }

      // zero suppression
      assert(m_nPostSample >= 1);
      pass_zero_suppression.resize(t_sorted_hits.size(), false);
      for (int i = 0; i < (int) t_sorted_hits.size(); ++i)
      {
        if (t_sorted_hits[i] > ADCThreshold_mV)
        {
          for (int j = i - m_nPreSample; j < i + (int)m_nPostSample; ++j)
          {
            if (j < 0) continue;
            if (j >= (int)pass_zero_suppression.size()) continue;

            pass_zero_suppression[j] = true;
          }
          i += m_nPostSample - 1;
        }
      }
      for (unsigned int i = 0; i < t_sorted_hits.size(); ++i)
      {
        if (not pass_zero_suppression[i]) t_sorted_hits[i] = 0;
      }

      // mV -> ADC
      for (auto &adc_bin : t_sorted_hits)
      {
        // input voltage x 1024 channels over 2200 mV max range
        adc_bin = static_cast<TpcDefs::ADCDataType>(adc_bin * 1024.0 / 2200.0);

        if (adc_bin < 0) adc_bin = 0;
        if (adc_bin > 1023) adc_bin = 1023;
      }  //       for (auto & adc_bin : t_sorted_hits)

    }  //  for (auto & t_sorted_hits : hitset->getTimeFrameAdcData())

  }  //   for (TrkrHitSetContainer::ConstIterator hitset_iter = hitset_range.first;

  return;
}

float PHG4TpcDigitizer::add_noise_to_bin(float signal)
{
  // add noise to the signal and return adc input voltage
  float adc_input_voltage = signal * ADCSignalConversionGain;                 // mV, see comments above
  float noise_voltage = (Pedestal + added_noise()) * ADCNoiseConversionGain;  // mV - from definition of noise charge and pedestal charge
  adc_input_voltage += noise_voltage;

  return adc_input_voltage;
}

float PHG4TpcDigitizer::added_noise()
{
  float noise = gsl_ran_gaussian(RandomGenerator, TpcEnc);

  return noise;
}

void PHG4TpcDigitizer::
SetTpcMinLayer(const int )
{
  std::cout <<__PRETTY_FUNCTION__<<" is deprecated. Ignore this call"<<std::endl;
}
