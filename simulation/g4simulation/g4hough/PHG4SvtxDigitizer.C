#include "PHG4SvtxDigitizer.h"

#include <g4main/PHG4Hit.h>

#include "SvtxHitMap.h"
#include "SvtxHitMap_v1.h"
#include "SvtxHit.h"
#include "SvtxHit_v1.h"

#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/getClass.h>
#include <g4detectors/PHG4CylinderCellGeomContainer.h>
#include <g4detectors/PHG4CylinderCellGeom.h>
#include <g4detectors/PHG4CylinderGeomContainer.h>
#include <g4detectors/PHG4CylinderGeom.h>

#include <g4detectors/PHG4Cell.h>
#include <g4detectors/PHG4Cellv1.h>
#include <g4detectors/PHG4Cellv2.h>
#include <g4detectors/PHG4CellContainer.h>
#include <g4detectors/PHG4CellDefs.h>
#include <phool/PHRandomSeed.h>

#include <gsl/gsl_randist.h>

#include <iostream>
#include <cmath>
#include <limits>

using namespace std;

PHG4SvtxDigitizer::PHG4SvtxDigitizer(const string &name) :
  SubsysReco(name),
  TPCMinLayer(7),
  ADCThreshold(2700),   // electrons
  TPCEnc(670),        // electrons
  Pedestal(50000),    // electrons
  ChargeToPeakVolts(20),    // mV/fC
  ADCSignalConversionGain(numeric_limits<float>::signaling_NaN()), // will be assigned in PHG4SvtxDigitizer::InitRun
  ADCNoiseConversionGain(numeric_limits<float>::signaling_NaN()), // will be assigned in PHG4SvtxDigitizer::InitRun
  _hitmap(NULL),
  _timer(PHTimeServer::get()->insert_new(name)) {

  unsigned int seed = PHRandomSeed();  // fixed seed is handled in this funtcion
  cout << Name() << " random seed: " << seed << endl;
  RandomGenerator = gsl_rng_alloc(gsl_rng_mt19937);
  gsl_rng_set(RandomGenerator, seed);

  if(Verbosity() > 0)
    cout << "Creating PHG4SvtxDigitizer with name = " << name << endl;
}

int PHG4SvtxDigitizer::InitRun(PHCompositeNode* topNode) {

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
  PHCompositeNode *dstNode 
    = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode","DST"));
  if (!dstNode) {
    cout << PHWHERE << "DST Node missing, doing nothing." << endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  PHNodeIterator iter_dst(dstNode);
    
  // Create the SVX node if required
  PHCompositeNode* svxNode 
    = dynamic_cast<PHCompositeNode*>(iter_dst.findFirst("PHCompositeNode","SVTX"));
  if (!svxNode) {
    svxNode = new PHCompositeNode("SVTX");
    dstNode->addNode(svxNode);
  }
  
  // Create the Hit node if required
  SvtxHitMap *svxhits = findNode::getClass<SvtxHitMap>(dstNode,"SvtxHitMap");
  if (!svxhits) {
    svxhits = new SvtxHitMap_v1();
    PHIODataNode<PHObject> *SvtxHitMapNode =
      new PHIODataNode<PHObject>(svxhits, "SvtxHitMap", "PHObject");
    svxNode->addNode(SvtxHitMapNode);
  }

  CalculateCylinderCellADCScale(topNode);
//  CalculateLadderCellADCScale(topNode); // obsolete, use PHG4SiliconTrackerDigitizer
  CalculateMapsLadderCellADCScale(topNode);
  
  //----------------
  // Report Settings
  //----------------
  
  if (Verbosity() > 0) {
    cout << "====================== PHG4SvtxDigitizer::InitRun() =====================" << endl;
    for (std::map<int,unsigned int>::iterator iter = _max_adc.begin();
	 iter != _max_adc.end();
	 ++iter) {
      cout << " Max ADC in Layer #" << iter->first << " = " << iter->second << endl;
    }
    for (std::map<int,float>::iterator iter = _energy_scale.begin();
	 iter != _energy_scale.end();
	 ++iter) {
      cout << " Energy per ADC in Layer #" << iter->first << " = " << 1.0e6*iter->second << " keV" << endl;
    }
    cout << "===========================================================================" << endl;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int PHG4SvtxDigitizer::process_event(PHCompositeNode *topNode) {

  _timer.get()->restart();

  _hitmap = findNode::getClass<SvtxHitMap>(topNode,"SvtxHitMap");
  if (!_hitmap) 
    {
      cout << PHWHERE << " ERROR: Can't find SvtxHitMap." << endl;
      return Fun4AllReturnCodes::ABORTRUN;
    }

  //Jin: don't clear up node. Fun4all server does that. Extra cleaning usually cause problems
//  _hitmap->Reset();
  
  DigitizeCylinderCells(topNode);
//  DigitizeLadderCells(topNode);  // obsolete, use PHG4SiliconTrackerDigitizer
  DigitizeMapsLadderCells(topNode);

  PrintHits(topNode);
  
  _timer.get()->stop();
  return Fun4AllReturnCodes::EVENT_OK;
}

void PHG4SvtxDigitizer::CalculateCylinderCellADCScale(PHCompositeNode *topNode) {

  // defaults to 8-bit ADC, short-axis MIP placed at 1/4 dynamic range

  PHG4CylinderCellGeomContainer *geom_container = findNode::getClass<PHG4CylinderCellGeomContainer>(topNode,"CYLINDERCELLGEOM_SVTX");
    

  if (!geom_container) return;
  
  PHG4CylinderCellGeomContainer::ConstRange layerrange = geom_container->get_begin_end();
  for(PHG4CylinderCellGeomContainer::ConstIterator layeriter = layerrange.first;
      layeriter != layerrange.second;
      ++layeriter) {

    int layer = layeriter->second->get_layer();
    float thickness = (layeriter->second)->get_thickness();
    float pitch = (layeriter->second)->get_phistep()*(layeriter->second)->get_radius();
    float length = (layeriter->second)->get_zstep();
   
    float minpath = pitch;
    if (length < minpath) minpath = length;
    if (thickness < minpath) minpath = thickness;
    float mip_e = 0.003876*minpath;  

    if (_max_adc.find(layer) == _max_adc.end()) {
      _max_adc[layer] = 255;
      _energy_scale[layer] = mip_e / 64;
    }
  }    

  return;
}

//void PHG4SvtxDigitizer::CalculateLadderCellADCScale(PHCompositeNode *topNode) {
//
//  // defaults to 8-bit ADC, short-axis MIP placed at 1/4 dynamic range
//
//  PHG4CylinderGeomContainer *geom_container = findNode::getClass<PHG4CylinderGeomContainer>(topNode,"CYLINDERGEOM_SILICON_TRACKER");
//
//  if (!geom_container) return;
//
//  PHG4CylinderGeomContainer::ConstRange layerrange = geom_container->get_begin_end();
//  for(PHG4CylinderGeomContainer::ConstIterator layeriter = layerrange.first;
//      layeriter != layerrange.second;
//      ++layeriter) {
//
//    int layer = layeriter->second->get_layer();
//    float thickness = (layeriter->second)->get_thickness();
//    float pitch = (layeriter->second)->get_strip_y_spacing();
//    float length = (layeriter->second)->get_strip_z_spacing();
//
//    float minpath = pitch;
//    if (length < minpath) minpath = length;
//    if (thickness < minpath) minpath = thickness;
//    float mip_e = 0.003876*minpath;
//
//    if (_max_adc.find(layer) == _max_adc.end()) {
//      _max_adc[layer] = 255;
//      _energy_scale[layer] = mip_e / 64;
//    }
//  }
//
//  return;
//}

void PHG4SvtxDigitizer::CalculateMapsLadderCellADCScale(PHCompositeNode *topNode) {

  // defaults to 8-bit ADC, short-axis MIP placed at 1/4 dynamic range

  PHG4CylinderGeomContainer *geom_container = findNode::getClass<PHG4CylinderGeomContainer>(topNode,"CYLINDERGEOM_MAPS");
    
  if (!geom_container) return;

  if(Verbosity())
    cout << "Found CYLINDERGEOM_MAPS node" << endl;
  
  PHG4CylinderGeomContainer::ConstRange layerrange = geom_container->get_begin_end();
  for(PHG4CylinderGeomContainer::ConstIterator layeriter = layerrange.first;
      layeriter != layerrange.second;
      ++layeriter) {

    int layer = layeriter->second->get_layer();
    float thickness = (layeriter->second)->get_pixel_thickness();
    float pitch = (layeriter->second)->get_pixel_x();
    float length = (layeriter->second)->get_pixel_z();
   
    float minpath = pitch;
    if (length < minpath) minpath = length;
    if (thickness < minpath) minpath = thickness;
    float mip_e = 0.003876*minpath;  

    if (Verbosity())
    cout << "mip_e = " << mip_e << endl;

    if (_max_adc.find(layer) == _max_adc.end()) {
      _max_adc[layer] = 255;
      _energy_scale[layer] = mip_e / 64;
    }
  }    

  return;
}

void PHG4SvtxDigitizer::DigitizeCylinderCells(PHCompositeNode *topNode) {

  //unsigned int print_layer = 100; // to suppress diagnostic output
  unsigned int print_layer = 47;  // to print diagnostic output for layer 47

  // Digitizes the TPC cells that were created in PHG4CylinderCellTPCReco
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
  // See https://indico.cern.ch/event/489996/timetable/#all.detailed "SAMPA Chip: the New ASIC for the ALICE TPC and MCH Upgrades", M Bregant
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
  // the voltage at the ADC input right then the values of Q_(pad,z) have to be scaled up by 2.4
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
  // increasing time should be (497 -> 249) and (0 -> 248)

  //----------
  // Get Nodes
  //----------

 PHG4CylinderCellGeomContainer* geom_container =
    findNode::getClass<PHG4CylinderCellGeomContainer>(topNode,"CYLINDERCELLGEOM_SVTX");
  if (!geom_container) {
    std::cout << PHWHERE << "ERROR: Can't find node CYLINDERCELLGEOM_SVTX" << std::endl;
  }
 
  PHG4CellContainer* cells = findNode::getClass<PHG4CellContainer>(topNode,"G4CELL_SVTX");
  if (!cells) return; 

  // sort the cells by layer
  // start with an empty vector of vectors of cells for each layer
  layer_sorted_cells.clear();

 PHG4CylinderCellGeomContainer::ConstRange layerrange = geom_container->get_begin_end();
  for(PHG4CylinderCellGeomContainer::ConstIterator layeriter = layerrange.first;
      layeriter != layerrange.second;
      ++layeriter) 
    {
      // add an empty vector of cells for this layer
      layer_sorted_cells.push_back(std::vector<const  PHG4Cell *>());
    }
	
  // now we fill each of the empty vectors with the cells for that layer
  PHG4CellContainer::ConstRange cellrange = cells->getCells();
  for(PHG4CellContainer::ConstIterator celliter = cellrange.first;
      celliter != cellrange.second;
      ++celliter) 
    {    
      PHG4Cell* cell =  celliter->second; 

      if(Verbosity() > 100)
	if( (unsigned int) cell->get_layer() == print_layer)
	  {
	    for (PHG4Cell::EdepConstIterator g4iter = cell->get_g4hits().first;
		 g4iter != cell->get_g4hits().second;
		 ++g4iter) 
	      {
		cout << "Digitizer: input cellid " << cell->get_cellid() << " g4hit ID " << g4iter->first << endl;
	      }
	  }


      if( (unsigned int) cell->get_layer() < TPCMinLayer) 
	{
	  if(Verbosity()>0) std::cout << "Skipping layer " << cell->get_layer() << std::endl;
	  continue;
	}
      //cout << "Digitizer: cell identify:" << endl;
      //cell->identify();
      layer_sorted_cells[cell->get_layer()-TPCMinLayer].push_back(cell);
     }

  // We have the cells sorted by layer, now we loop over the layers and process the hits
  //==========================================================

  for(PHG4CylinderCellGeomContainer::ConstIterator layeriter = layerrange.first;
      layeriter != layerrange.second;
      ++layeriter) 
    {
      // for this layer, make a vector of a vector of cells for each phibin
      phi_sorted_cells.clear();
      
      // start with an empty vector of cells for each phibin
      unsigned int layer = layeriter->second->get_layer();
      int nphibins = layeriter->second->get_phibins();
       for(int iphi = 0;iphi<nphibins;iphi++)
	{
	  phi_sorted_cells.push_back( std::vector<const  PHG4Cell*>() );      
	}

      // Fill the vector of cells for each phibin
      for(unsigned int i = 0; i < layer_sorted_cells[layer-TPCMinLayer].size(); ++i) 
	{
	  unsigned int phibin = PHG4CellDefs::SizeBinning::get_phibin(layer_sorted_cells[layer-TPCMinLayer][i]->get_cellid());	  
	  phi_sorted_cells[phibin].push_back(layer_sorted_cells[layer-TPCMinLayer][i]);
	}
      
      // For this layer we have the cells sorted into vectors for each phi      
      // process these vectors one phi bin at a time
      for(unsigned int iphi=0;iphi<phi_sorted_cells.size();iphi++)
	{	 
	  if( phi_sorted_cells[iphi].size() == 0)
	    continue;
 
	  // Populate a vector of cells ordered by Z for each phibin    
	  int nzbins = layeriter->second->get_zbins();
	  is_populated.clear();
	  is_populated.assign(nzbins,2);  // mark all as noise only for now
	  z_sorted_cells.clear();
	 
	  // add an empty vector for each z bin
	  for(int iz=0;iz<nzbins;iz++)
	    z_sorted_cells.push_back( std::vector<const  PHG4Cell*>() );
 
	  // add a cell for each z bin that has one
	  for(unsigned int iz=0;iz<phi_sorted_cells[iphi].size();iz++)
	    {
	      int zbin = PHG4CellDefs::SizeBinning::get_zbin(phi_sorted_cells[iphi][iz]->get_cellid());
	      is_populated[zbin] = 1;  // this bin is a associated with a cell
	      z_sorted_cells[zbin].push_back(phi_sorted_cells[iphi][iz]);
	    }
	  
	  adc_input.clear(); 
	  adc_cellid.clear(); 
	  // Now for this phibin we process the cells ordered by Z bin into hits with noise
	  //======================================================
	  // For this step we take the edep value and convert it to mV at the ADC input
	  // See comments above for how to do this for signal and noise
	  for(int iz=0;iz<nzbins;iz++)
	    {    
	      if(is_populated[iz]==1)
		{
		  // This zbin has a filled cell, add noise
		  float noise = added_noise();  // in electrons
		  float noise_voltage = (Pedestal + noise) * ADCNoiseConversionGain;  // mV - from definition of noise charge and pedestal charge
		  float adc_input_voltage = z_sorted_cells[iz][0]->get_edep() * ADCSignalConversionGain;  // mV, see comments above

		  adc_input.push_back(adc_input_voltage + noise_voltage);
		  adc_cellid.push_back(z_sorted_cells[iz][0]->get_cellid());
		}
	      else if(is_populated[iz]==2)
		{
		  // This z bin does not have a filled cell, add noise
		  float noise = added_noise();  // returns "electrons"
		  float noise_voltage = (Pedestal + noise) * ADCNoiseConversionGain; // mV - noise "electrons" x conversion gain

		  adc_input.push_back(noise_voltage); // mV
		  adc_cellid.push_back(0);  // there is no cell, just add a placeholder in the vector for now, replace it later
		}
	      else
		{
		  // Cannot happen
		  cout << "Impossible value of is_populated, iz = " << iz << " is_populated = " << is_populated[iz] << endl;
		  exit(-1);
		}
	    }
	  
          // Now we can digitize the stream of z bins 

	  // start with negative z, the first to arrive is bin 0
	  //================================
	  int binpointer = 0;
	  for(int iz=0;iz<nzbins/2;iz++) // 0-247
	    {
	      if(iz < binpointer) continue;
	      
	      if(adc_input[iz] > ADCThreshold*ADCNoiseConversionGain)  // convert threshold in "equivalent electrons" to mV
		{
		  // digitize this bin and the following 4 bins


		  if(Verbosity() > 100)
		    if(layer == print_layer) cout << endl << "  (neg z) Above threshold of " << ADCThreshold*ADCNoiseConversionGain << " for phibin " << iphi 
						  << " iz " << iz << " with adc_input " << adc_input[iz] << " digitize this and 4 following bins: "  << endl;

		  for(int izup=0;izup<5; izup++)
		    {
		      if(iz+izup < nzbins/2 && iz + izup >= 0)
			{			  
			  unsigned int adc_output = (unsigned int) (adc_input[iz+izup] * 1024.0  / 2200.0);  // input voltage x 1024 channels over 2200 mV max range
			  if(adc_input[iz+izup] < 0) adc_output = 0;
			  if (adc_output > 1023) adc_output = 1023;

			  if(Verbosity() > 100)
			    if(layer == print_layer)  cout << "    (neg z) iz+izup " << iz+izup << " adc_cellid " << adc_cellid[iz+izup] 
							   << "  adc_input "  << adc_input[iz+izup] << " ADCThreshold " << ADCThreshold*ADCNoiseConversionGain 
							   << " adc_output " << adc_output << endl;
			  
			  if(is_populated[iz+izup] == 2)
			    {
			      // This is a noise-only hit, so there is no cell
			      // since it is digitized, we will make a hit for it
			      // but first, we have to add a cell for it, so things don't break downstream
			      //cout << " generate key for noise hit" << endl;
			      PHG4CellDefs::keytype akey = PHG4CellDefs::SizeBinning::genkey(layer, iz+izup, iphi);
			      PHG4Cell *cell = new PHG4Cellv1(akey);

			      cell->add_edep(adc_input[iz+izup] / ADCSignalConversionGain);  //convert from voltage back to electrons from GEM
			      adc_cellid[iz+izup]=cell->get_cellid();

			      if(Verbosity() > 100)
				if(layer == print_layer)  cout << " will digitize noise hit for iphi " << iphi << " zbin " << iz+izup 
							       << " created new cell with cellid " << cell->get_cellid() 
							       << " adc_input " << adc_input[iz+izup] << " edep " << cell->get_edep() << endl; 

			      cells->AddCell(cell);
			    }

			  SvtxHit_v1 hit;
			  
			  hit.set_layer(layer);
			  hit.set_cellid(adc_cellid[iz+izup]);
			  
			  float e = adc_output * (2200 / 1024.0) /  ADCSignalConversionGain;          // convert ADU's back to mV, then to input electrons
			  hit.set_adc(adc_output);
			  hit.set_e(e);

			  SvtxHit* ptr = _hitmap->insert(&hit);      
			  if (!ptr->isValid()) 
			    {
			      static bool first = true;
			      if (first) 
				{
				  cout << PHWHERE << "ERROR: Incomplete SvtxHits are being created" << endl;
				  ptr->identify();
				  first = false;
				}
			    }
			  binpointer ++;

			  if(Verbosity() > 100)
			    if(layer == print_layer)
			      { 
				//cout << endl << "Digitizer: Hit identify after insertion:" << endl;
				//ptr->identify();
				PHG4Cell *tmpcell =  cells->findCell(ptr->get_cellid());
				cout << "   Digitizer (neg z): Hit " << ptr->get_id() << " is from cellid " << tmpcell->get_cellid() << " with contributing g4hits: " << endl;
				// list the contrubuting g4 hits for this cell - loop over all the g4hits
				for (PHG4Cell::EdepConstIterator g4iter = tmpcell->get_g4hits().first;
				     g4iter != tmpcell->get_g4hits().second;
				     ++g4iter) 
				  {
				    cout << "       g4hitID " << g4iter->first << " in layer " << tmpcell->get_layer() << " with edep " << g4iter->second << endl; 
				  }
			      }
			  
			} // end nzbins check 
		    } // end izup loop
		  
		}  //  adc threshold if 		  
	      else 
		{
		  // below threshold, move on
		  binpointer++;
		} // end adc threshold if/else 		  
	      
	    } // end iz loop

	  // now positive z, the first to arrive is bin 497
	  //===============================
	  binpointer = nzbins-1;
	  for(int iz=nzbins-1;iz>=nzbins/2;iz--) // 495 - 248
	    {
	      if(iz > binpointer) continue;
	      
	      if(adc_input[iz] > ADCThreshold* ADCNoiseConversionGain)  // convert threshold in electrons to mV
		{
		  // digitize this bin and the following 4 bins
		  
		  if(Verbosity() > 100)
		    if(layer == print_layer) cout << endl << "  (pos z) Above threshold  of " << ADCThreshold*ADCNoiseConversionGain << " for iz " << iz 
						  << " with adc_input " << adc_input[iz] << " digitize this and 4 following bins: " << endl;
		  
		  for(int izup=0;izup<5; izup++)
		    {
		      if(iz-izup < nzbins && iz - izup >= nzbins/2)
			{			  

			  unsigned int adc_output = (unsigned int) (adc_input[iz-izup] * 1024.0 / 2200.0);  // input voltage x 1024 channels over 2200 mV max range
			  if(adc_input[iz-izup] < 0) adc_output = 0;
			  if (adc_output > 1023) adc_output = 1023;

			  if(Verbosity() > 100)
			    if(layer == print_layer)  cout << "    (pos z) iz-izup " << iz-izup << " adc_cellid " << adc_cellid[iz-izup] 
							   << "  adc_input "  << adc_input[iz-izup] << " ADCThreshold " << ADCThreshold*ADCNoiseConversionGain 
							   << " adc_output " << adc_output << endl;
			  
			  if(is_populated[iz-izup] == 2)
			    {
			      // This is a noise-only hit, so there is no cell
			      // since it is digitized, we will make a hit for it
			      // first, we have to add a cell for it so things don't break downstream

			      PHG4CellDefs::keytype akey = PHG4CellDefs::SizeBinning::genkey(layer, iz-izup, iphi);
			      PHG4Cell *cell = new PHG4Cellv1(akey);

			      cell->add_edep(adc_input[iz-izup] / ADCSignalConversionGain);  //convert from voltage back to electrons from GEM stack 
			      adc_cellid[iz-izup]=cell->get_cellid();
			      if(Verbosity() > 100)
				if(layer == print_layer)  cout  << " will digitize noise hit for iphi " << iphi << " zbin " << iz-izup 
								<< " created new cell with cellid " << cell->get_cellid() 
								<< " edep " << cell->get_edep() << endl; 
			      
			      cells->AddCell(cell);
			    }
			  
			  SvtxHit_v1 hit;
			  
			  hit.set_layer(layer);
			  hit.set_cellid(adc_cellid[iz-izup]);
			  
			  float e = adc_output * (2200.0 / 1024.0) /  ADCSignalConversionGain;          // convert ADU's back to mV, then to input electrons from GEM
			  hit.set_adc(adc_output);
			  hit.set_e(e);
 
			  SvtxHit* ptr = _hitmap->insert(&hit);      

			  if (!ptr->isValid()) 
			    {
			      static bool first = true;
			      if (first) 
				{
				  cout << PHWHERE << "ERROR: Incomplete SvtxHits are being created" << endl;
				  ptr->identify();
				  first = false;
				}
			    }
			  binpointer--;

			  if(Verbosity() > 100)
			    if(layer == print_layer)
			      {
				//cout << endl << "Digitizer (pos z): Hit identify after insertion:" << endl;
				//ptr->identify();
				
				PHG4Cell *tmpcell =  cells->findCell(ptr->get_cellid());
				cout << "   Digitizer (pos z): Hit " << ptr->get_id() << " is from cellid " << tmpcell->get_cellid() << " with contributing g4hits: " << endl;			      
				// list the contrubuting g4 hits for this cell -  loop over all the g4hits
				for (PHG4Cell::EdepConstIterator g4iter = tmpcell->get_g4hits().first;
				     g4iter != tmpcell->get_g4hits().second;
				     ++g4iter) 
				  {
				    cout << "       g4hitID " << g4iter->first << " in layer " << tmpcell->get_layer() << " with edep " << g4iter->second << endl; 
				  }
			      }
			  
			} // end nzbins check 
		    } // end izup loop
		  
		}  //  adc threshold if 		  
	      else 
		{
		  // below threshold, move on
		  binpointer--;
		} // end adc threshold if/else 		  
	      
	    } // end iz loop
	  
	  
	} // end phibins loop
      
      
    } // end loop over layers
  
  return;
}

//void PHG4SvtxDigitizer::DigitizeLadderCells(PHCompositeNode *topNode) {
//
//  //----------
//  // Get Nodes
//  //----------
//
//  PHG4CellContainer* cells = findNode::getClass<PHG4CellContainer>(topNode,"G4CELL_SILICON_TRACKER");
//  if (!cells) return;
//
//  //-------------
//  // Digitization
//  //-------------
//
//  vector<PHG4Cell*> cell_list;
//  PHG4CellContainer::ConstRange cellrange = cells->getCells();
//  for(PHG4CellContainer::ConstIterator celliter = cellrange.first;
//      celliter != cellrange.second;
//      ++celliter) {
//
//    PHG4Cell* cell = celliter->second;
//
//    SvtxHit_v1 hit;
//
//    hit.set_layer(cell->get_layer());
//    hit.set_cellid(cell->get_cellid());
//
//    unsigned int adc = cell->get_edep() / _energy_scale[hit.get_layer()];
//    if (adc > _max_adc[hit.get_layer()]) adc = _max_adc[hit.get_layer()];
//    float e = _energy_scale[hit.get_layer()] * adc;
//
//    hit.set_adc(adc);
//    hit.set_e(e);
//
//    SvtxHit* ptr = _hitmap->insert(&hit);
//    if (!ptr->isValid()) {
//      static bool first = true;
//      if (first) {
//	cout << PHWHERE << "ERROR: Incomplete SvtxHits are being created" << endl;
//	ptr->identify();
//	first = false;
//      }
//    }
//  }
//
//  return;
//}

void PHG4SvtxDigitizer::DigitizeMapsLadderCells(PHCompositeNode *topNode) {

  //----------
  // Get Nodes
  //----------
 
  PHG4CellContainer* cells = findNode::getClass<PHG4CellContainer>(topNode,"G4CELL_MAPS");
  if (!cells) return; 
  
  //-------------
  // Digitization
  //-------------

  vector<PHG4Cell*> cell_list;
  PHG4CellContainer::ConstRange cellrange = cells->getCells();
  for(PHG4CellContainer::ConstIterator celliter = cellrange.first;
      celliter != cellrange.second;
      ++celliter) {
    
    PHG4Cell* cell = celliter->second;
    
    SvtxHit_v1 hit;

    hit.set_layer(cell->get_layer());
    hit.set_cellid(cell->get_cellid());

    unsigned int adc = cell->get_edep() / _energy_scale[hit.get_layer()];
    if (adc > _max_adc[hit.get_layer()]) adc = _max_adc[hit.get_layer()]; 
    float e = _energy_scale[hit.get_layer()] * adc;
    
    hit.set_adc(adc);
    hit.set_e(e);
        
    SvtxHit* ptr = _hitmap->insert(&hit);      
    if (!ptr->isValid()) {
      static bool first = true;
      if (first) {
	cout << PHWHERE << "ERROR: Incomplete SvtxHits are being created" << endl;
	ptr->identify();
	first = false;
      }
    }
  }
  
  return;
}

void PHG4SvtxDigitizer::PrintHits(PHCompositeNode *topNode) {

  if (Verbosity() >= 1) {

    SvtxHitMap *hitlist = findNode::getClass<SvtxHitMap>(topNode,"SvtxHitMap");
    if (!hitlist) return;
    
    cout << "================= PHG4SvtxDigitizer::process_event() ====================" << endl;
  

    cout << " Found and recorded the following " << hitlist->size() << " hits: " << endl;

    unsigned int ihit = 0;
    for (SvtxHitMap::Iter iter = hitlist->begin();
	 iter != hitlist->end();
	 ++iter) {

      SvtxHit* hit = iter->second;
      cout << ihit << " of " << hitlist->size() << endl;
      hit->identify();
      ++ihit;
    }
    
    cout << "===========================================================================" << endl;
  }
  
  return;
}

float PHG4SvtxDigitizer::added_noise()
{
  float noise =  gsl_ran_gaussian(RandomGenerator, TPCEnc);

  return noise;
}
