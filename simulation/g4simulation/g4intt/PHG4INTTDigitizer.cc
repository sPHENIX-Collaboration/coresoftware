#include "PHG4INTTDigitizer.h"

#include "INTTDeadMap.h"

#include <fun4all/Fun4AllReturnCodes.h>
#include <g4detectors/PHG4Cell.h>
#include <g4detectors/PHG4CellContainer.h>
#include <g4detectors/PHG4CylinderCellGeom.h>
#include <g4detectors/PHG4CylinderCellGeomContainer.h>
#include <g4detectors/PHG4CylinderGeom.h>
#include <g4detectors/PHG4CylinderGeomContainer.h>

#include <trackbase_historic/SvtxHit.h>
#include <trackbase_historic/SvtxHitMap.h>
#include <trackbase_historic/SvtxHitMap_v1.h>
#include <trackbase_historic/SvtxHit_v1.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHRandomSeed.h>
#include <phool/getClass.h>

#include <gsl/gsl_randist.h>

#include <cmath>
#include <iostream>

using namespace std;

PHG4INTTDigitizer::PHG4INTTDigitizer(const string &name)
  : SubsysReco(name)
  , mNoiseMean(457.2)
  , mNoiseSigma(166.6)
  , mEnergyPerPair(3.62e-9)  // GeV/e-h
  , _hitmap(NULL)
  , m_nCells(0)
  , m_nDeadCells(0)
{
  unsigned int seed = PHRandomSeed();  // fixed seed is handled in this funtcion
  cout << Name() << " random seed: " << seed << endl;
  RandomGenerator = gsl_rng_alloc(gsl_rng_mt19937);
  gsl_rng_set(RandomGenerator, seed);
}

int PHG4INTTDigitizer::InitRun(PHCompositeNode *topNode)
{
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

  // Create the SVX node if required
  PHCompositeNode *svxNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "SVTX"));
  if (!svxNode)
  {
    svxNode = new PHCompositeNode("SVTX");
    dstNode->addNode(svxNode);
  }

  // Create the Hit node if required
  SvtxHitMap *svxhits = findNode::getClass<SvtxHitMap>(topNode, "SvtxHitMap");
  if (!svxhits)
  {
    svxhits = new SvtxHitMap_v1();
    PHIODataNode<PHObject> *SvtxHitMapNode =
        new PHIODataNode<PHObject>(svxhits, "SvtxHitMap", "PHObject");
    svxNode->addNode(SvtxHitMapNode);
  }

  CalculateLadderCellADCScale(topNode);

  //----------------
  // Report Settings
  //----------------

  if (Verbosity() > 0)
  {
    cout << "====================== PHG4INTTDigitizer::InitRun() =====================" << endl;
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

int PHG4INTTDigitizer::process_event(PHCompositeNode *topNode)
{
  _hitmap = findNode::getClass<SvtxHitMap>(topNode, "SvtxHitMap");
  if (!_hitmap)
  {
    cout << PHWHERE << " ERROR: Can't find SvtxHitMap." << endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  // _hitmap->Reset();

  DigitizeLadderCells(topNode);

  PrintHits(topNode);

  return Fun4AllReturnCodes::EVENT_OK;
}

void PHG4INTTDigitizer::CalculateLadderCellADCScale(PHCompositeNode *topNode)
{
  // FPHX 3-bit ADC, thresholds are set in "set_fphx_adc_scale".

  PHG4CellContainer *cells = findNode::getClass<PHG4CellContainer>(topNode, "G4CELL_INTT");
  PHG4CylinderGeomContainer *geom_container = findNode::getClass<PHG4CylinderGeomContainer>(topNode, "CYLINDERGEOM_INTT");

  if (!geom_container || !cells) return;

  PHG4CylinderGeomContainer::ConstRange layerrange = geom_container->get_begin_end();
  for (PHG4CylinderGeomContainer::ConstIterator layeriter = layerrange.first;
       layeriter != layerrange.second;
       ++layeriter)
  {
    int layer = layeriter->second->get_layer();
    if (_max_fphx_adc.find(layer) == _max_fphx_adc.end())
      assert(!"Error: _max_fphx_adc is not available.");

    float thickness = (layeriter->second)->get_thickness();  // cm
    float mip_e = 0.003876 * thickness;                      // GeV
    _energy_scale.insert(std::make_pair(layer, mip_e));
  }

  return;
}

void PHG4INTTDigitizer::DigitizeLadderCells(PHCompositeNode *topNode)
{
  //----------
  // Get Nodes
  //----------

  PHG4CellContainer *cells = findNode::getClass<PHG4CellContainer>(topNode, "G4CELL_INTT");
  if (!cells) return;

  const INTTDeadMap *deadmap = findNode::getClass<INTTDeadMap>(topNode, "DEADMAP_INTT");
  if (Verbosity() >= VERBOSITY_MORE)
  {
    if (deadmap)
    {
      cout << "PHG4INTTDigitizer::DigitizeLadderCells - Use deadmap ";
      deadmap->identify();
    }
    else
    {
      cout << "PHG4INTTDigitizer::DigitizeLadderCells - Can not find deadmap, all channels enabled " << endl;
    }
  }

  //-------------
  // Digitization
  //-------------

  PHG4CellContainer::ConstRange cellrange = cells->getCells();
  for (PHG4CellContainer::ConstIterator celliter = cellrange.first;
       celliter != cellrange.second;
       ++celliter)
  {
    PHG4Cell *cell = celliter->second;

    ++m_nCells;
    if (deadmap)
    {
      if (deadmap->isDeadChannelINTT(
              cell->get_layer(),             //const int layer,
              cell->get_ladder_phi_index(),  //const int ladder_phi,
              cell->get_ladder_z_index(),    //const int ladder_z,
              cell->get_zbin(),              //const int strip_z,
              cell->get_phibin()             //const int strip_phi
              ))
      {
        ++m_nDeadCells;
        if (Verbosity() >= VERBOSITY_MORE)
        {
          cout << "PHG4INTTDigitizer::DigitizeLadderCells - dead cell at layer " << cell->get_layer() << ": ";
          cell->identify();
        }
        continue;
      }
    }  //    if (deadmap)

    SvtxHit_v1 hit;

    const int layer = cell->get_layer();

    hit.set_layer(layer);
    hit.set_cellid(cell->get_cellid());

    if (_energy_scale.count(layer) > 1)
      assert(!"Error: _energy_scale has two or more keys.");

    // Convert Geant4 true cell energy to # of electrons and add noise electrons
    const int n_true_electron = (int) (cell->get_edep() / mEnergyPerPair);
    const int n_noise_electron = round(added_noise());
    const int n_cell_electron = n_true_electron + n_noise_electron;

    const float mip_e = _energy_scale[layer];

    std::vector<std::pair<double, double> > vadcrange = _max_fphx_adc[layer];

    int adc = -1;
    for (unsigned int irange = 0; irange < vadcrange.size(); ++irange)
    {
      // Convert adc ranges from fraction to the MIP energy to # of electrons.
      // vadcrange uses FLT_MAX and the order of mEnergyPerPair is e-9, so use double.
      const double n_adcrange_electron_first = vadcrange[irange].first * mip_e / mEnergyPerPair;
      const double n_adcrange_electron_second = vadcrange[irange].second * mip_e / mEnergyPerPair;

      if (n_cell_electron >= n_adcrange_electron_first && n_cell_electron < n_adcrange_electron_second)
        adc = (int) irange;
    }
    //
    if (adc >= 0)
    {
      //      adc = 0;

      float e = 0.0;
      if (adc >= 0 && adc < int(vadcrange.size()) - 1)
        e = 0.5 * (vadcrange[adc].second + vadcrange[adc].first) * mip_e;
      else if (adc == int(vadcrange.size()) - 1)  // overflow
        e = vadcrange[adc].first * mip_e;
      else  // underflow
        e = 0.5 * vadcrange[0].first * mip_e;

      hit.set_adc(adc);
      hit.set_e(e);

      SvtxHit *ptr = _hitmap->insert(&hit);
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
    }
  }

  return;
}

//! end of process
int PHG4INTTDigitizer::End(PHCompositeNode *topNode)
{
  if (Verbosity() >= VERBOSITY_SOME)
  {
    cout << "PHG4INTTDigitizer::End - processed "
         << m_nCells << " cell with "
         << m_nDeadCells << " dead cells masked"
         << " (" << 100. * m_nDeadCells / m_nCells << "%)" << endl;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

void PHG4INTTDigitizer::PrintHits(PHCompositeNode *topNode)
{
  if (Verbosity() >= VERBOSITY_EVEN_MORE)
  {
    //if (Verbosity() >= 0) {

    SvtxHitMap *hitlist = findNode::getClass<SvtxHitMap>(topNode, "SvtxHitMap");
    if (!hitlist) return;

    cout << "================= PHG4INTTDigitizer::process_event() ====================" << endl;

    cout << " Found and recorded the following " << hitlist->size() << " hits: " << endl;

    unsigned int ihit = 0;
    for (SvtxHitMap::Iter iter = hitlist->begin();
         iter != hitlist->end();
         ++iter)
    {
      SvtxHit *hit = iter->second;
      cout << ihit << " of " << hitlist->size() << endl;
      hit->identify();
      ++ihit;
    }

    cout << "===========================================================================" << endl;
  }

  return;
}

float PHG4INTTDigitizer::added_noise()
{
  float noise = gsl_ran_gaussian(RandomGenerator, mNoiseSigma) + mNoiseMean;
  noise = (noise < 0) ? 0 : noise;

  return noise;
}
