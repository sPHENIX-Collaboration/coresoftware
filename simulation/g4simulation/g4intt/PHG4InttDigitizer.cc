// This is the new trackbase container version

#include "PHG4InttDigitizer.h"

#include "InttDeadMap.h"

#include <g4detectors/PHG4CylinderCellGeom.h>
#include <g4detectors/PHG4CylinderCellGeomContainer.h>
#include <g4detectors/PHG4CylinderGeom.h>
#include <g4detectors/PHG4CylinderGeomContainer.h>

// Move to new storage containers
#include <trackbase/TrkrHitSet.h>
#include <trackbase/TrkrHitSetContainer.h>
#include <trackbase/TrkrHitTruthAssoc.h>
#include <trackbase/TrkrDefs.h>
#include <intt/InttDefs.h>
#include <intt/InttHit.h>

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHRandomSeed.h>
#include <phool/getClass.h>

#include <TSystem.h>

#include <gsl/gsl_randist.h>

#include <cfloat>
#include <cmath>
#include <iostream>
#include <cassert>

using namespace std;

PHG4InttDigitizer::PHG4InttDigitizer(const string &name)
  : SubsysReco(name)
  , PHParameterInterface(name)
  , detector("INTT")
  , mNoiseMean(457.2)
  , mNoiseSigma(166.6)
  , mEnergyPerPair(3.62e-9)  // GeV/e-h
  , m_nCells(0)
  , m_nDeadCells(0)
{
  InitializeParameters();
  unsigned int seed = PHRandomSeed();  // fixed seed is handled in this funtcion
  cout << Name() << " random seed: " << seed << endl;
  RandomGenerator = gsl_rng_alloc(gsl_rng_mt19937);
  gsl_rng_set(RandomGenerator, seed);
}

PHG4InttDigitizer::~PHG4InttDigitizer()
{
  gsl_rng_free(RandomGenerator);
}

int PHG4InttDigitizer::InitRun(PHCompositeNode *topNode)
{
  cout << "PHG4InttDigitizer::InitRun: detector = " << detector << endl;

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

  CalculateLadderCellADCScale(topNode);

  // Create the run and par nodes
  PHCompositeNode *runNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "RUN"));
  PHCompositeNode *parNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "PAR"));

  string paramnodename = "G4CELLPARAM_" + detector;
  string geonodename = "G4CELLGEO_" + detector;

  UpdateParametersWithMacro();
  // save this to the run wise tree to store on DST
  PHNodeIterator runIter(runNode);
  PHCompositeNode *RunDetNode = dynamic_cast<PHCompositeNode *>(runIter.findFirst("PHCompositeNode", detector));
  if (!RunDetNode)
  {
    RunDetNode = new PHCompositeNode(detector);
    runNode->addNode(RunDetNode);
  }
  SaveToNodeTree(RunDetNode, paramnodename);
  // save this to the parNode for use
  PHNodeIterator parIter(parNode);
  PHCompositeNode *ParDetNode = dynamic_cast<PHCompositeNode *>(parIter.findFirst("PHCompositeNode", detector));
  if (!ParDetNode)
  {
    ParDetNode = new PHCompositeNode(detector);
    parNode->addNode(ParDetNode);
  }
  PutOnParNode(ParDetNode, geonodename);
  mNoiseMean = get_double_param("NoiseMean");
  mNoiseSigma = get_double_param("NoiseSigma");
  mEnergyPerPair = get_double_param("EnergyPerPair");

  //----------------
  // Report Settings
  //----------------

  if (Verbosity() > 0)
  {
    cout << "====================== PHG4InttDigitizer::InitRun() =====================" << endl;
    for (std::map<int, unsigned int>::iterator iter1 = _max_adc.begin();
         iter1 != _max_adc.end();
         ++iter1)
    {
      cout << " Max ADC in Layer #" << iter1->first << " = " << iter1->second << endl;
    }
    for (std::map<int, float>::iterator iter2 = _energy_scale.begin();
         iter2 != _energy_scale.end();
         ++iter2)
    {
      cout << " Energy per ADC in Layer #" << iter2->first << " = " << 1.0e6 * iter2->second << " keV" << endl;
    }
    cout << "===========================================================================" << endl;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int PHG4InttDigitizer::process_event(PHCompositeNode *topNode)
{
  DigitizeLadderCells(topNode);

  return Fun4AllReturnCodes::EVENT_OK;
}

void PHG4InttDigitizer::CalculateLadderCellADCScale(PHCompositeNode *topNode)
{
  // FPHX 3-bit ADC, thresholds are set in "set_fphx_adc_scale".

  //PHG4CellContainer *cells = findNode::getClass<PHG4CellContainer>(topNode, "G4CELL_INTT");
  PHG4CylinderGeomContainer *geom_container = findNode::getClass<PHG4CylinderGeomContainer>(topNode, "CYLINDERGEOM_INTT");

  //if (!geom_container || !cells) return;
  if (!geom_container) return;

  PHG4CylinderGeomContainer::ConstRange layerrange = geom_container->get_begin_end();
  for (PHG4CylinderGeomContainer::ConstIterator layeriter = layerrange.first;
       layeriter != layerrange.second;
       ++layeriter)
  {
    int layer = layeriter->second->get_layer();
    if (_max_fphx_adc.find(layer) == _max_fphx_adc.end())
    {
      cout << "Error: _max_fphx_adc is not available." << endl;
      gSystem->Exit(1);
    }
    float thickness = (layeriter->second)->get_thickness();  // cm
    float mip_e = 0.003876 * thickness;                      // GeV
    _energy_scale.insert(std::make_pair(layer, mip_e));
  }

  return;
}

void PHG4InttDigitizer::DigitizeLadderCells(PHCompositeNode *topNode)
{
  //---------------------------
  // Get common Nodes
  //---------------------------
  const InttDeadMap *deadmap = findNode::getClass<InttDeadMap>(topNode, "DEADMAP_INTT");
  if (Verbosity() >= VERBOSITY_MORE)
  {
    if (deadmap)
    {
      cout << "PHG4InttDigitizer::DigitizeLadderCells - Use deadmap ";
      deadmap->identify();
    }
    else
    {
      cout << "PHG4InttDigitizer::DigitizeLadderCells - Can not find deadmap, all channels enabled " << endl;
    }
  }

  // Get the TrkrHitSetContainer node
  TrkrHitSetContainer *trkrhitsetcontainer = findNode::getClass<TrkrHitSetContainer>(topNode, "TRKR_HITSET");
  if(!trkrhitsetcontainer)
    {
      cout << "Could not locate TRKR_HITSET node, quit! " << endl;
      exit(1);
    }

 //-------------
  // Digitization
  //-------------

  // We want all hitsets for the Intt
  TrkrHitSetContainer::ConstRange hitset_range = trkrhitsetcontainer->getHitSets(TrkrDefs::TrkrId::inttId);
  for (TrkrHitSetContainer::ConstIterator hitset_iter = hitset_range.first;
       hitset_iter != hitset_range.second;
       ++hitset_iter)
    {
     // we have an itrator to one TrkrHitSet for the intt from the trkrHitSetContainer
      // get the hitset key so we can find the layer
      TrkrDefs::hitsetkey hitsetkey = hitset_iter->first;
      const int layer = TrkrDefs::getLayer(hitsetkey);
      const int ladder_phi = InttDefs::getLadderPhiId(hitsetkey);
      const int ladder_z = InttDefs::getLadderZId(hitsetkey);

      if(Verbosity() > 1) 
	cout << "PHG4InttDigitizer: found hitset with key: " << hitsetkey << " in layer " << layer << endl;

      // get all of the hits from this hitset      
      TrkrHitSet *hitset = hitset_iter->second;
      TrkrHitSet::ConstRange hit_range = hitset->getHits();
      for(TrkrHitSet::ConstIterator hit_iter = hit_range.first;
	  hit_iter != hit_range.second;
	  ++hit_iter)
	{
	  ++m_nCells;

	  TrkrHit *hit = (InttHit*) hit_iter->second;
	  TrkrDefs::hitkey hitkey = hit_iter->first;
	  int strip_col =  InttDefs::getCol(hitkey);  // strip z index
	  int strip_row =   InttDefs::getRow(hitkey);  // strip phi index

	  // Apply deadmap here if desired
	  if (deadmap)
	    {
	      if (deadmap->isDeadChannelIntt(
					     layer, 
					     ladder_phi,
					     ladder_z,
					     strip_col,
					     strip_row
					     ))
		{
		  ++m_nDeadCells;
		  if (Verbosity() >= VERBOSITY_MORE)
		    {
		      cout << "PHG4InttDigitizer::DigitizeLadderCells - dead strip at layer " << layer << ": ";
		      hit->identify();
		    }
		  continue;
		}
	    }  //    if (deadmap)

	  if (_energy_scale.count(layer) > 1)
	    assert(!"Error: _energy_scale has two or more keys.");

	  const float mip_e = _energy_scale[layer];

	  std::vector<std::pair<double, double> > vadcrange = _max_fphx_adc[layer];

	  int adc = -1;
	  for (unsigned int irange = 0; irange < vadcrange.size(); ++irange)
	    if (hit->getEnergy() >= vadcrange[irange].first * (double) mip_e && hit->getEnergy() < vadcrange[irange].second * (double) mip_e)
	      adc = (int) irange;

	  if(adc == -1)
	    // how do we specify underflow or overflow?
	    adc = 0;
	  
	  hit->setAdc(adc);	      

	  if(Verbosity() > 2)
	    cout << "PHG4InttDigitizer: found hit with layer "  << layer << " ladder_z " << ladder_z << " ladder_phi " << ladder_phi 
		 << " strip_col " << strip_col << " strip_row " << strip_row << " adc " << adc << endl;
 
	} // end loop over hits in this hitset
    } // end loop over hitsets
  
  return;
}

//! end of process
int PHG4InttDigitizer::End(PHCompositeNode *topNode)
{
  if (Verbosity() >= VERBOSITY_SOME)
  {
    cout << "PHG4InttDigitizer::End - processed "
         << m_nCells << " cell with "
         << m_nDeadCells << " dead cells masked"
         << " (" << 100. * m_nDeadCells / m_nCells << "%)" << endl;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

void PHG4InttDigitizer::SetDefaultParameters()
{
  set_default_double_param("NoiseMean", 457.2);
  set_default_double_param("NoiseSigma", 166.6);
  set_default_double_param("EnergyPerPair", 3.62e-9);  // GeV/e-h
  return;
}

float PHG4InttDigitizer::added_noise()
{
//  float noise = gsl_ran_gaussian(RandomGenerator, mNoiseSigma) + mNoiseMean;
//  noise = (noise < 0) ? 0 : noise;

  // Note the noise is bi-polar, i.e. can make ths signal fluctuate up and down.
  // Much of the mNoiseSigma as extracted in https://github.com/sPHENIX-Collaboration/coresoftware/pull/580
  // is statistical fluctuation from the limited calibration data. They does not directly apply here.
  float noise = gsl_ran_gaussian(RandomGenerator, mNoiseMean);

  return noise;
}

void PHG4InttDigitizer::set_adc_scale(const int &layer, const std::vector<double> &userrange)
{
  if (userrange.size() != nadcbins)
  {
    cout << "Error: vector in set_fphx_adc_scale(vector) must have eight elements." << endl;
    gSystem->Exit(1);
  }
  //sort(userrange.begin(), userrange.end()); // TODO, causes GLIBC error

  std::vector<std::pair<double, double> > vadcrange;
  for (unsigned int irange = 0; irange < userrange.size(); ++irange)
  {
    if (irange == userrange.size() - 1)
    {
      vadcrange.push_back(std::make_pair(userrange[irange], FLT_MAX));
    }
    else
    {
      vadcrange.push_back(std::make_pair(userrange[irange], userrange[irange + 1]));
    }
  }
  _max_fphx_adc.insert(std::make_pair(layer, vadcrange));
}
