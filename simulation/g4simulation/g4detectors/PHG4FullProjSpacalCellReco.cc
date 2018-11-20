#include "PHG4FullProjSpacalCellReco.h"

#include "PHG4CylinderGeomContainer.h"
#include "PHG4CylinderGeom_Spacalv3.h"
#include "PHG4CylinderCellGeomContainer.h"
#include "PHG4CylinderCellGeom.h"
#include "PHG4CylinderCellGeom_Spacalv1.h"

#include "PHG4CellContainer.h"
#include "PHG4CellDefs.h"
#include "PHG4Cellv1.h"

#include <g4main/PHG4Hit.h>
#include <g4main/PHG4HitContainer.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/Fun4AllServer.h>

#include <phool/PHNodeIterator.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/getClass.h>

#include <TFile.h>
#include <TH1.h>
#include <TH2.h>

#include <boost/foreach.hpp>

#include <cassert>
#include <cmath>
#include <cstdlib>
#include <exception>
#include <iostream>
#include <limits>       // std::numeric_limits
#include <sstream>


using namespace std;

PHG4FullProjSpacalCellReco::PHG4FullProjSpacalCellReco(const string &name) :
  SubsysReco(name),
  PHParameterInterface(name),
  _timer(PHTimeServer::get()->insert_new(name.c_str())), 
  sum_energy_g4hit(0),
  chkenergyconservation(0),
  tmin(NAN),
  tmax(NAN),
  light_collection_model()
{
  InitializeParameters();
}

int
PHG4FullProjSpacalCellReco::ResetEvent(PHCompositeNode *topNode)
{
  sum_energy_g4hit = 0.;
  return Fun4AllReturnCodes::EVENT_OK;
}

int
PHG4FullProjSpacalCellReco::InitRun(PHCompositeNode *topNode)
{
  PHNodeIterator iter(topNode);

  // Looking for the DST node
  PHCompositeNode *dstNode;
  dstNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
    {
      std::cout << PHWHERE << "DST Node missing, doing nothing." << std::endl;
      exit(1);
    }
  hitnodename = "G4HIT_" + detector;
  PHG4HitContainer *g4hit = findNode::getClass<PHG4HitContainer>(topNode, hitnodename.c_str());
  if (!g4hit)
    {
      cout << "PHG4FullProjSpacalCellReco::InitRun - Could not locate g4 hit node "
           << hitnodename << endl;
      topNode->print();
      exit(1);
    }
  cellnodename = "G4CELL_" + detector;
  PHG4CellContainer *cells = findNode::getClass<PHG4CellContainer>(topNode, cellnodename);
  if (!cells)
    {
      PHNodeIterator dstiter(dstNode);
      PHCompositeNode *DetNode = dynamic_cast<PHCompositeNode*>(dstiter.findFirst("PHCompositeNode", detector));
      if (!DetNode)
        {
          DetNode = new PHCompositeNode(detector);
          dstNode->addNode(DetNode);
        }
      cells = new PHG4CellContainer();
      PHIODataNode<PHObject> *newNode = new PHIODataNode<PHObject>(cells, cellnodename.c_str(), "PHObject");
      DetNode->addNode(newNode);
    }

  geonodename = "CYLINDERGEOM_" + detector;
  PHG4CylinderGeomContainer *geo =
      findNode::getClass<PHG4CylinderGeomContainer>(topNode, geonodename.c_str());
  if (!geo)
    {
      cout << "PHG4FullProjSpacalCellReco::InitRun - Could not locate geometry node "
           << geonodename << endl;
      topNode->print();
      exit(1);
    }
  if (Verbosity() > 0)
    {
      cout << "PHG4FullProjSpacalCellReco::InitRun - incoming geometry:"
           << endl;
      geo->identify();
    }
  seggeonodename = "CYLINDERCELLGEOM_" + detector;
  PHG4CylinderCellGeomContainer *seggeo = findNode::getClass<PHG4CylinderCellGeomContainer>(topNode, seggeonodename.c_str());
  if (!seggeo)
    {
      seggeo = new PHG4CylinderCellGeomContainer();
      PHCompositeNode *runNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "RUN"));
      PHIODataNode<PHObject> *newNode = new PHIODataNode<PHObject>(seggeo, seggeonodename.c_str(), "PHObject");
      runNode->addNode(newNode);
    }

  const PHG4CylinderGeom *layergeom_raw = geo->GetFirstLayerGeom();
  assert(layergeom_raw);

  // a special implimentation of PHG4CylinderGeom is required here.
  const PHG4CylinderGeom_Spacalv3 *layergeom =
    dynamic_cast<const PHG4CylinderGeom_Spacalv3 *>(layergeom_raw);

  if (!layergeom)
    {
      cout << "PHG4FullProjSpacalCellReco::InitRun - Fatal Error -"
	   << " require to work with a version of SPACAL geometry of PHG4CylinderGeom_Spacalv3 or higher. "
	   << "However the incoming geometry has version "
	   << layergeom_raw->ClassName() << endl;
      exit(1);
    }
  if (Verbosity() > 1)
    {
      layergeom->identify();
    }

  layergeom->subtower_consistency_check();

  //      int layer = layergeom->get_layer();

  // create geo object and fill with variables common to all binning methods
  PHG4CylinderCellGeom_Spacalv1 *layerseggeo =
    new PHG4CylinderCellGeom_Spacalv1();
  layerseggeo->set_layer(layergeom->get_layer());
  layerseggeo->set_radius(layergeom->get_radius());
  layerseggeo->set_thickness(layergeom->get_thickness());
  layerseggeo->set_binning(PHG4CellDefs::spacalbinning);

  // construct a map to convert tower_ID into the older eta bins.

  const PHG4CylinderGeom_Spacalv3::tower_map_t & tower_map =
    layergeom->get_sector_tower_map();
  const PHG4CylinderGeom_Spacalv3::sector_map_t & sector_map =
    layergeom->get_sector_map();
  const int nphibin = layergeom->get_azimuthal_n_sec() // sector
    * layergeom->get_max_phi_bin_in_sec() // blocks per sector
    * layergeom->get_n_subtower_phi(); // subtower per block
  const double deltaphi = 2. * M_PI / nphibin;

  typedef map<double, int> map_z_tower_z_ID_t;
  map_z_tower_z_ID_t map_z_tower_z_ID;
  double phi_min = NAN;

  BOOST_FOREACH(const PHG4CylinderGeom_Spacalv3::tower_map_t::value_type& tower_pair, tower_map)
    {

      const int & tower_ID = tower_pair.first;
      const PHG4CylinderGeom_Spacalv3::geom_tower & tower =
	tower_pair.second;

      // inspect index in sector 0
      pair<int, int> tower_z_phi_ID = layergeom->get_tower_z_phi_ID(tower_ID, 0);

      const int & tower_ID_z = tower_z_phi_ID.first;
      const int & tower_ID_phi = tower_z_phi_ID.second;

      if (tower_ID_phi == 0)
	{
	  //assign phi min according phi bin 0
	  phi_min = M_PI_2 - deltaphi *(layergeom->get_max_phi_bin_in_sec()* layergeom->get_n_subtower_phi()/2) // shift of first tower in sector
	    + sector_map.begin()->second;
	}

      if (tower_ID_phi == layergeom->get_max_phi_bin_in_sec() / 2)
	{
	  //assign eta min according phi bin 0
	  map_z_tower_z_ID[tower.centralZ] = tower_ID_z;
	}
      // ...
    } //       BOOST_FOREACH(const PHG4CylinderGeom_Spacalv3::tower_map_t::value_type& tower_pair, tower_map)

  assert(! std::isnan(phi_min));
  layerseggeo->set_phimin(phi_min);
  layerseggeo->set_phistep(deltaphi);
  layerseggeo->set_phibins(nphibin);

  PHG4CylinderCellGeom_Spacalv1::tower_z_ID_eta_bin_map_t tower_z_ID_eta_bin_map;
  int eta_bin = 0;
  BOOST_FOREACH( map_z_tower_z_ID_t::value_type& z_tower_z_ID, map_z_tower_z_ID)
    {
      tower_z_ID_eta_bin_map[z_tower_z_ID.second] = eta_bin;
      eta_bin++;
    }
  layerseggeo->set_tower_z_ID_eta_bin_map(tower_z_ID_eta_bin_map);
  layerseggeo->set_etabins(eta_bin * layergeom->get_n_subtower_eta());
  layerseggeo->set_etamin(NAN);
  layerseggeo->set_etastep(NAN);

  //build eta bin maps
  BOOST_FOREACH(const PHG4CylinderGeom_Spacalv3::tower_map_t::value_type& tower_pair, tower_map)
    {

      const int & tower_ID = tower_pair.first;
      const PHG4CylinderGeom_Spacalv3::geom_tower & tower =
	tower_pair.second;

      // inspect index in sector 0
      std::pair<int, int> tower_z_phi_ID = layergeom->get_tower_z_phi_ID(tower_ID, 0);
      const int & tower_ID_z = tower_z_phi_ID.first;
      const int & tower_ID_phi = tower_z_phi_ID.second;
      const int & etabin = tower_z_ID_eta_bin_map[tower_ID_z];

      if (tower_ID_phi == layergeom->get_max_phi_bin_in_sec() / 2)
	{
	  // half z-range
	  const double dz = fabs(0.5 * (tower.pDy1 + tower.pDy2) / sin(tower.pRotationAngleX));
	  const double tower_radial = layergeom->get_tower_radial_position(tower);

	  auto z_to_eta = [&tower_radial](const double &z){return -log(tan(0.5 * atan2(tower_radial, z)));};

	  const double eta_central = z_to_eta(tower.centralZ);
	  // half eta-range
	  const double deta = (z_to_eta( tower.centralZ + dz) - z_to_eta( tower.centralZ - dz))/2;
	  assert(deta > 0);

	  for (int sub_tower_ID_y = 0; sub_tower_ID_y < tower.NSubtowerY;
	       ++sub_tower_ID_y)
	    {
	      assert(tower.NSubtowerY <=layergeom->get_n_subtower_eta());
	      // do not overlap to the next bin.
	      const int sub_tower_etabin = etabin * layergeom->get_n_subtower_eta() + sub_tower_ID_y;

	      const pair<double, double>etabounds  (eta_central - deta + sub_tower_ID_y * 2 * deta / tower.NSubtowerY,
            eta_central - deta + (sub_tower_ID_y + 1) * 2 * deta / tower.NSubtowerY);

	      const pair<double, double>zbounds  (tower.centralZ - dz + sub_tower_ID_y * 2 * dz / tower.NSubtowerY,
            tower.centralZ - dz + (sub_tower_ID_y + 1) * 2 * dz / tower.NSubtowerY);

	      layerseggeo->set_etabounds(sub_tower_etabin,etabounds);
	      layerseggeo->set_zbounds(sub_tower_etabin,zbounds);

	      if (Verbosity() >= VERBOSITY_SOME)
	        {
	          cout << "PHG4FullProjSpacalCellReco::InitRun::" << Name()
               << "\t tower_ID_z = " << tower_ID_z
               << "\t tower_ID_phi = " << tower_ID_phi
               << "\t sub_tower_ID_y = " << sub_tower_ID_y
               << "\t sub_tower_etabin = " << sub_tower_etabin
               << "\t dz = " << dz
               << "\t tower_radial = " << tower_radial
               << "\t eta_central = " << eta_central
               << "\t deta = " << deta
               << "\t etabounds = [" <<etabounds.first << ", " << etabounds.second<<"]"
               << "\t zbounds = [" <<zbounds.first << ", " << zbounds.second<<"]"
               <<endl;
	        }

	    }

	}
      // ...
    } //       BOOST_FOREACH(const PHG4CylinderGeom_Spacalv3::tower_map_t::value_type& tower_pair, tower_map)

      // add geo object filled by different binning methods
  seggeo->AddLayerCellGeom(layerseggeo);
  if (Verbosity() >= VERBOSITY_SOME)
    {
      cout << "PHG4FullProjSpacalCellReco::InitRun::" << Name()
	   << " - Done layer" << (layergeom->get_layer())
	   << ". Print out the cell geometry:" << endl;
      layerseggeo->identify();
    }
  UpdateParametersWithMacro();
  // save this to the run wise tree to store on DST
  PHCompositeNode *runNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "RUN" ));
  PHCompositeNode *parNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "PAR" ));
  PHNodeIterator runIter(runNode);
  PHCompositeNode *RunDetNode =  dynamic_cast<PHCompositeNode*>(runIter.findFirst("PHCompositeNode",detector));
  if (! RunDetNode)
    {
      RunDetNode = new PHCompositeNode(detector);
      runNode->addNode(RunDetNode);
    }
  string paramnodename = "G4CELLPARAM_" + detector;
  SaveToNodeTree(RunDetNode,paramnodename);
  // save this to the parNode for use
  PHNodeIterator parIter(parNode);
  PHCompositeNode *ParDetNode =  dynamic_cast<PHCompositeNode*>(parIter.findFirst("PHCompositeNode",detector));
  if (! ParDetNode)
    {
      ParDetNode = new PHCompositeNode(detector);
      parNode->addNode(ParDetNode);
    }
  string geonodename = "G4CELLGEO_" + detector;
  PutOnParNode(ParDetNode,geonodename);
  tmin = get_double_param("tmin");
  tmax = get_double_param("tmax");

  return Fun4AllReturnCodes::EVENT_OK;
}

int
PHG4FullProjSpacalCellReco::process_event(PHCompositeNode *topNode)
{
  _timer.get()->restart();
  PHG4HitContainer *g4hit = findNode::getClass<PHG4HitContainer>(topNode, hitnodename.c_str());
  if (!g4hit)
    {
      cout << "PHG4FullProjSpacalCellReco::process_event - Fatal Error - Could not locate g4 hit node "
           << hitnodename << endl;
      exit(1);
    }
  PHG4CellContainer *cells = findNode::getClass<PHG4CellContainer>(topNode, cellnodename);
  if (!cells)
    {
      cout << "PHG4FullProjSpacalCellReco::process_event - Fatal Error - could not locate cell node "
           << cellnodename << endl;
      exit(1);
    }

  PHG4CylinderCellGeomContainer *seggeo = findNode::getClass<PHG4CylinderCellGeomContainer>(topNode, seggeonodename.c_str());
  if (!seggeo)
    {
      cout << "PHG4FullProjSpacalCellReco::process_event - Fatal Error - could not locate cell geometry node "
           << seggeonodename << endl;
      exit(1);
    }

  PHG4CylinderGeomContainer *layergeo = findNode::getClass<PHG4CylinderGeomContainer>(topNode, geonodename.c_str());
  if (!layergeo)
    {
      cout << "PHG4FullProjSpacalCellReco::process_event - Fatal Error - Could not locate sim geometry node "
           << geonodename << endl;
      exit(1);
    }

    
  PHG4HitContainer::ConstIterator hiter;
  PHG4HitContainer::ConstRange hit_begin_end = g4hit->getHits();

  PHG4CylinderCellGeom *geo_raw = seggeo->GetFirstLayerCellGeom();
  PHG4CylinderCellGeom_Spacalv1 * geo = dynamic_cast<PHG4CylinderCellGeom_Spacalv1 *>(geo_raw);
  assert(geo);
  const PHG4CylinderGeom *layergeom_raw = layergeo->GetFirstLayerGeom();
  assert(layergeom_raw);
  // a special implimentation of PHG4CylinderGeom is required here.
  const PHG4CylinderGeom_Spacalv3 *layergeom =
    dynamic_cast<const PHG4CylinderGeom_Spacalv3 *>(layergeom_raw);
  assert(layergeom);

  for (hiter = hit_begin_end.first; hiter != hit_begin_end.second; ++hiter)
    {
      // checking ADC timing integration window cut
      if (hiter->second->get_t(0)>tmax) continue;
      if (hiter->second->get_t(1)<tmin) continue;

      sum_energy_g4hit += hiter->second->get_edep();

      // hit loop
      int scint_id = hiter->second->get_scint_id();

      // decode scint_id
      PHG4CylinderGeom_Spacalv3::scint_id_coder decoder(scint_id);

      // convert to z_ID, phi_ID
      std::pair<int, int> tower_z_phi_ID = layergeom->get_tower_z_phi_ID(decoder.tower_ID, decoder.sector_ID);
      const int & tower_ID_z = tower_z_phi_ID.first;
      const int & tower_ID_phi = tower_z_phi_ID.second;

      PHG4CylinderGeom_Spacalv3::tower_map_t::const_iterator it_tower =
	layergeom->get_sector_tower_map().find(decoder.tower_ID);
      assert(it_tower != layergeom->get_sector_tower_map().end());

      unsigned int key = static_cast<unsigned int>(scint_id);
      PHG4Cell *cell = nullptr;
      map<unsigned int, PHG4Cell *>::iterator it = celllist.find(key);
      if (it != celllist.end())
	{
	  cell = it->second;
	}
      else
	{

	  // convert tower_ID_z to to eta bin number
	  int etabin = -1;
	  try
	    {
	      etabin = geo->get_etabin_block(tower_ID_z); // block eta bin
	    }
	  catch (exception & e)
	    {
	      cout << "Print cell geometry:" << endl;
	      geo->identify();
	      cout << "Print scint_id_coder:" << endl;
	      decoder.identify();
	      cout << "Print the hit:" << endl;
	      hiter->second->print();
	      cout << "PHG4FullProjSpacalCellReco::process_event::"
		   << Name() << " - Fatal Error - " << e.what() << endl;
	      exit(1);
	    }

	  const int sub_tower_ID_x = it_tower->second.get_sub_tower_ID_x(decoder.fiber_ID);
	  const int sub_tower_ID_y = it_tower->second.get_sub_tower_ID_y(decoder.fiber_ID);
	  unsigned short fiber_ID = decoder.fiber_ID;
	  unsigned short etabinshort  =  etabin * layergeom->get_n_subtower_eta() + sub_tower_ID_y;
	  unsigned short phibin = tower_ID_phi * layergeom->get_n_subtower_phi() + sub_tower_ID_x;
	  PHG4CellDefs::keytype cellkey = PHG4CellDefs::SpacalBinning::genkey(etabinshort,phibin,fiber_ID);
	  cell = new PHG4Cellv1(cellkey);
	  celllist[key] = cell;
	}

      double light_yield = hiter->second->get_light_yield();

      // light yield correction from fiber attenuation:
      if (light_collection_model.use_fiber_model())
	{
	  const double z = 0.5
	    * (hiter->second->get_local_z(0)
	       + hiter->second->get_local_z(1));
	  assert(not std::isnan(z));

	  light_yield *= light_collection_model.get_fiber_transmission(z);
	}

      // light yield correction from light guide collection efficiency:
      if (light_collection_model.use_fiber_model())
	{
	  const double x = it_tower->second.get_position_fraction_x_in_sub_tower(decoder.fiber_ID);
	  const double y = it_tower->second.get_position_fraction_y_in_sub_tower(decoder.fiber_ID);

	  light_yield *= light_collection_model.get_light_guide_efficiency(x, y);
	}

      cell->add_edep(hiter->first, hiter->second->get_edep());
      cell->add_edep(hiter->second->get_edep());
      cell->add_light_yield(light_yield);
      cell->add_shower_edep(hiter->second->get_shower_id(), hiter->second->get_edep());

    } // end loop over g4hits
  int numcells = 0;
  for (map<unsigned int, PHG4Cell *>::const_iterator mapiter =
	 celllist.begin(); mapiter != celllist.end(); ++mapiter)
    {
      cells->AddCell(mapiter->second);
      numcells++;
      if (Verbosity() > 1)
	{
	  cout << "PHG4FullProjSpacalCellReco::process_event::" << Name()
	       << " - " << "Adding cell in bin eta "
	       << PHG4CellDefs::SpacalBinning::get_etabin(mapiter->second->get_cellid()) 
	       << " phi "
	       << PHG4CellDefs::SpacalBinning::get_phibin(mapiter->second->get_cellid()) 
	       << " fiber "
	       << PHG4CellDefs::SpacalBinning::get_fiberid(mapiter->second->get_cellid()) 
	       << ", energy dep: "
	       << mapiter->second->get_edep() << ", light yield: "
	       << mapiter->second->get_light_yield() << endl;
	}
    }
  celllist.clear();
  if (Verbosity() > 0)
    {
      cout << "PHG4FullProjSpacalCellReco::process_event::" << Name()
	   << " - " << " found " << numcells
	   << " fibers with energy deposition" << endl;
    }
    

  if (chkenergyconservation || Verbosity() > 4)
    {
      CheckEnergy(topNode);
    }
  _timer.get()->stop();
  return Fun4AllReturnCodes::EVENT_OK;
}

int
PHG4FullProjSpacalCellReco::CheckEnergy(PHCompositeNode *topNode)
{
  PHG4CellContainer *cells = findNode::getClass<PHG4CellContainer>(topNode, cellnodename);
  double sum_energy_cells = 0.;
  PHG4CellContainer::ConstRange cell_begin_end = cells->getCells();
  PHG4CellContainer::ConstIterator citer;
  for (citer = cell_begin_end.first; citer != cell_begin_end.second; ++citer)
    {
      sum_energy_cells += citer->second->get_edep();

    }
  // the fractional eloss for particles traversing eta bins leads to minute rounding errors
  if (fabs(sum_energy_cells - sum_energy_g4hit) / sum_energy_g4hit > 1e-6)
    {
      cout
          << "PHG4FullProjSpacalCellReco::CheckEnergy - energy mismatch between cells: "
          << sum_energy_cells << " and hits: " << sum_energy_g4hit
          << " diff sum(cells) - sum(hits): "
          << sum_energy_cells - sum_energy_g4hit << endl;
      return -1;
    }
  else
    {
      if (Verbosity() > 0)
        {
          cout << "PHG4FullProjSpacalCellReco::CheckEnergy::" << Name()
              << " - total energy for this event: " << sum_energy_g4hit
              << " GeV. Passed CheckEnergy" << endl;
        }
    }
  return 0;
}

PHG4FullProjSpacalCellReco::LightCollectionModel::LightCollectionModel() :
    data_grid_light_guide_efficiency(nullptr), 
    data_grid_fiber_trans(nullptr)
{

  data_grid_light_guide_efficiency_verify = new TH2F("data_grid_light_guide_efficiency_verify",
      "light collection efficiency as used in PHG4FullProjSpacalCellReco::LightCollectionModel;x positio fraction;y position fraction", //
      100, 0., 1., 100, 0., 1.);

  data_grid_fiber_trans_verify = new TH1F("data_grid_fiber_trans",
      "SCSF-78 Fiber Transmission as used in PHG4FullProjSpacalCellReco::LightCollectionModel;position in fiber (cm);Effective transmission",
      100, -15, 15);

  // register histograms
  Fun4AllServer *se = Fun4AllServer::instance();

  se->registerHisto(data_grid_light_guide_efficiency_verify);
  se->registerHisto(data_grid_fiber_trans_verify);

}

PHG4FullProjSpacalCellReco::LightCollectionModel::~LightCollectionModel()
{
  delete data_grid_light_guide_efficiency;
  delete data_grid_fiber_trans;
}

void
PHG4FullProjSpacalCellReco::LightCollectionModel::load_data_file(
    const std::string & input_file,
    const std::string & histogram_light_guide_model,
    const std::string & histogram_fiber_model)
{
  TFile * fin = TFile::Open(input_file.c_str());

  assert(fin);
  assert(fin->IsOpen());

  if (data_grid_light_guide_efficiency) delete data_grid_light_guide_efficiency;
  data_grid_light_guide_efficiency = dynamic_cast<TH2 *>(fin->Get(histogram_light_guide_model.c_str()));
  assert(data_grid_light_guide_efficiency);
  data_grid_light_guide_efficiency->SetDirectory(nullptr);

  if (data_grid_fiber_trans) delete data_grid_fiber_trans;
  data_grid_fiber_trans = dynamic_cast<TH1 *>(fin->Get(histogram_fiber_model.c_str()));
  assert(data_grid_fiber_trans);
  data_grid_fiber_trans->SetDirectory(nullptr);

  delete fin;
}

double
PHG4FullProjSpacalCellReco::LightCollectionModel::get_light_guide_efficiency(const double x_fraction, const double y_fraction)
{
  assert(data_grid_light_guide_efficiency);
  assert(x_fraction >= 0);
  assert(x_fraction <= 1);
  assert(y_fraction >= 0);
  assert(y_fraction <= 1);

  const double eff = data_grid_light_guide_efficiency->Interpolate(x_fraction,
      y_fraction);

  data_grid_light_guide_efficiency_verify->SetBinContent( //
      data_grid_light_guide_efficiency_verify->GetXaxis()->FindBin(x_fraction), //
      data_grid_light_guide_efficiency_verify->GetYaxis()->FindBin(y_fraction), //
      eff //
      );

  return eff;
}

double
PHG4FullProjSpacalCellReco::LightCollectionModel::get_fiber_transmission(const double z_distance)
{
  assert(data_grid_fiber_trans);

  const double eff = data_grid_fiber_trans->Interpolate(z_distance);

  data_grid_fiber_trans_verify->SetBinContent( //
      data_grid_fiber_trans_verify->GetXaxis()->FindBin(z_distance), //
      eff //
      );

  return eff;
}

void
PHG4FullProjSpacalCellReco::SetDefaultParameters()
{
  set_default_double_param("tmax",60.0);
  set_default_double_param("tmin",-20.0); // collision has a timing spread around the triggered event. Accepting negative time too.
  return;
}

void
PHG4FullProjSpacalCellReco::set_timing_window(const double tmi, const double tma)
{
  set_double_param("tmin",tmi);
  set_double_param("tmax",tma);
}
