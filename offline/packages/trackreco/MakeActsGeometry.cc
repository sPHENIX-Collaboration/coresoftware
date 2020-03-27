/*!
 *  \file		MakeActsGeometry.C
 *  \brief		Refit SvtxTracks with PHActs.
 *  \details	Refit SvtxTracks with PHActs.
 *  \author	        Tony Frawley <afrawley@fsu.edu>
 */

/*
//modularizing:

MakeActsGeometry  
  // Create acts geometry
  // Make (hitsetkey,surface) map and place on node tree

Create acts geometry (makes (hitsetkey,surface) maps)
  BuildSiliconLayers();  
      MakeActsGeometry  // makes (hitsetkey,surface) map for silicon
  BuildTpcSurfaceMap  // makes (hitsetkey,surface) map for silicon
  MakeTGeoNodeMap(topNode);  // makes cluster-node map - not used?
    isActive  // not used
*/

#include "MakeActsGeometry.h"

#include <trackbase/TrkrDefs.h>

#include <intt/CylinderGeomIntt.h>
#include <intt/InttDefs.h>

#include <mvtx/CylinderGeom_Mvtx.h>
#include <mvtx/MvtxDefs.h>

#include <tpc/TpcDefs.h>

#include <g4detectors/PHG4CylinderGeom.h>           // for PHG4CylinderGeom
#include <g4detectors/PHG4CylinderGeomContainer.h>
#include <g4detectors/PHG4CylinderCellGeom.h>
#include <g4detectors/PHG4CylinderCellGeomContainer.h>

#include <phgeom/PHGeomUtility.h>

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/getClass.h>
#include <phool/phool.h>

#include <Acts/Geometry/GeometryContext.hpp>
#include <Acts/Geometry/GeometryContext.hpp>
#include <Acts/Geometry/TrackingGeometry.hpp>
#include <Acts/Geometry/TrackingVolume.hpp>
#include <Acts/Surfaces/Surface.hpp>
#include <Acts/Surfaces/PlaneSurface.hpp>
#include <Acts/Surfaces/PerigeeSurface.hpp>
#include <Acts/EventData/TrackParameters.hpp>
#include <Acts/Utilities/CalibrationContext.hpp>
#include <Acts/MagneticField/MagneticFieldContext.hpp>

#include <ACTFW/Detector/IBaseDetector.hpp>
#include <ACTFW/EventData/Track.hpp>
#include <ACTFW/Framework/AlgorithmContext.hpp>
#include <ACTFW/Framework/IContextDecorator.hpp>
#include <ACTFW/Framework/WhiteBoard.hpp>
#include <ACTFW/Plugins/BField/BFieldOptions.hpp>
#include <ACTFW/Geometry/CommonGeometry.hpp>
#include <ACTFW/Options/CommonOptions.hpp>
#include <ACTFW/Plugins/Obj/ObjWriterOptions.hpp>

#include <ACTFW/Utilities/Options.hpp>


#include <TVector3.h>
#include <TMatrixT.h>                               // for TMatrixT, operator*
#include <TObject.h>
#include <TGeoManager.h>
#include <TSystem.h>
#include <cmath>                              // for sqrt, NAN
#include <cstddef>                                              // for size_t
#include <cstdlib>                                              // for atoi
#include <iostream>
#include <map>
#include <memory>
#include <utility>
#include <vector>

using namespace std;

/*
 * Constructor
 */
MakeActsGeometry::MakeActsGeometry(const string& name)
  : _geom_container_mvtx(nullptr)
  , _geom_container_intt(nullptr)
  , _geom_container_tpc(nullptr)
  , _geomanager(nullptr)
  , MinSurfZ(0.0)
  , MaxSurfZ(110.0)
  , NSurfZ(11)
  , NSurfPhi(10)
  ,_verbosity(0)
{

  // These are arbitrary subdivisions, and may change
  SurfStepZ = (MaxSurfZ - MinSurfZ) / (double) NSurfZ;
  ModuleStepPhi = 2.0 * M_PI / 12.0;
  ModulePhiStart = - M_PI;
  SurfStepPhi = 2.0 * M_PI / (double) (NSurfPhi * NTpcModulesPerLayer);

}

int MakeActsGeometry::BuildAllGeometry(PHCompositeNode *topNode)
{
  GetNodes(topNode);    // for geometry nodes

  CreateNodes(topNode);   // for writing surface map

  // run Acts layer builder
  BuildSiliconLayers();

  // create a map of sensor TGeoNode pointers using the TrkrDefs:: hitsetkey as the key  
  MakeTGeoNodeMap(topNode);

  // TPC continuous readout geometry does not exist within ACTS, so we build our own surfaces
  BuildTpcSurfaceMap();

 
  return Fun4AllReturnCodes::EVENT_OK;
}

void MakeActsGeometry::BuildTpcSurfaceMap()
{
  // Make a map of surfaces corresponding to each TPC sectorid and side.
  // There are 12 sectors and 2 sides of the membrane
  // We additionally subdivide these surfaces into approximately plane surfaces segmented in Z and phi
  // Subdivisions for creating TPC Acts::Surfaces are defined in constructor and header file
  for(unsigned int iz = 0; iz < NSurfZ; ++iz)
    {
      for(unsigned int side = 0; side < NTpcSides; ++side)
	{
	  double z_center =  SurfStepZ / 2.0 + (double) iz * SurfStepZ;
	  if(side == 0)  z_center = - z_center;

	  for(unsigned int ilayer = 0; ilayer < NTpcLayers; ++ilayer)
	    {
	      int tpc_layer = ilayer + 7;  // 3 MVTX and 4 INTT layers ahead of TPC

	      PHG4CylinderCellGeom *layergeom = _geom_container_tpc->GetLayerCellGeom(tpc_layer);
	      if(!layergeom)
		{
		  std::cout << PHWHERE << "Did not get layergeom for layer " <<tpc_layer  << std::endl;
		  gSystem->Exit(1);
		}
	      
	      double radius = layergeom->get_radius();

	      for(unsigned int imod = 0; imod < NTpcModulesPerLayer; ++imod)
		{
		  for(unsigned int iphi = 0; iphi < NSurfPhi; ++iphi)
		    {
		      double min_phi = ModulePhiStart + (double) imod * ModuleStepPhi + (double) iphi * SurfStepPhi;

		      double phi_center = min_phi +  SurfStepPhi / 2.0;		      
		      double x_center = radius * cos(phi_center);
		      double y_center = radius * sin(phi_center);

		      // There is a surface constructor that is a Vector3D for the center and a Vector3D for the normal. 
		      Acts::Vector3D center_vec(x_center, y_center, z_center);
		      
		      // The normal vector is a unit vector pointing to the center of the beam line at z = z_center
		      double x_norm = - x_center / sqrt(x_center*x_center + y_center*y_center);
		      double y_norm = -y_center / sqrt(x_center*x_center + y_center*y_center);
		      double z_norm = 0;
		      
		      Acts::Vector3D center_norm(x_norm, y_norm, z_norm);
		      
		      auto surf = Acts::Surface::makeShared<Acts::PlaneSurface>(center_vec, center_norm);

		      // make up a fake cluskey with the "cluster id" being composed of iphi and iz
		      unsigned int i_phi_z = iphi + 100*iz;
		      TrkrDefs::cluskey surfkey = TpcDefs::genClusKey(tpc_layer, imod, side, i_phi_z);

		      // Add this surface to the map
		      std::pair<TrkrDefs::cluskey, std::shared_ptr<const Acts::Surface>> tmp = make_pair(surfkey, surf);
		      _cluster_surface_map_tpc.insert(tmp);
		    }
		}		  
	    }
	}
    }
 }  


/**
 * Builds silicon layers in the ACTS surface world
 */
void MakeActsGeometry::BuildSiliconLayers()
{
  // define int argc and char* argv to provide options to processGeometry
  
  // Can hard code geometry options since the TGeo options are fixed by our detector design
  const int argc = 28;
  char *arg[argc];
  const std::string argstr[argc]{"-n1", "-l0", "--geo-tgeo-filename=none", "--geo-tgeo-worldvolume=\"World\"", "--geo-subdetectors", "MVTX", "Silicon", "--geo-tgeo-nlayers=0", "0", "--geo-tgeo-clayers=1",  "1", "--geo-tgeo-players=0", "0", "--geo-tgeo-clayernames", "MVTX", "siactive", "--geo-tgeo-cmodulenames", "MVTXSensor",  "siactive",  "--geo-tgeo-cmoduleaxes", "xzy", "yzx",  "--output-obj", "true","--bf-values",  "0", "0", "1.4"};

  // Set vector of chars to arguments needed
  for(int i = 0;i < argc; ++i)
    {
      // need a copy, since .c_str() returns a const char * and process geometry will not take a const
      arg[i] = strdup(argstr[i].c_str());
    }

  // We replicate the relevant functionality of  acts-framework/Examples/Common/src/GeometryExampleBase::ProcessGeometry() in MakeActsGeometry()
  // so we get access to the results. The layer builder magically gets the TGeoManager
  
  MakeSiliconGeometry(argc, arg, detector);

}

int MakeActsGeometry::MakeSiliconGeometry(int argc, char* argv[], FW::IBaseDetector& detector)
{
  // setup and parse options
  auto desc = FW::Options::makeDefaultOptions();
  FW::Options::addSequencerOptions(desc);
  FW::Options::addGeometryOptions(desc);
  FW::Options::addMaterialOptions(desc);
  FW::Options::addObjWriterOptions(desc);
  //  FW::Options::addCsvWriterOptions(desc);
  FW::Options::addOutputOptions(desc);
  FW::Options::addBFieldOptions(desc);

  // Add specific options for this geometry
  detector.addOptions(desc);
  auto vm = FW::Options::parse(desc, argc, argv);
  if (vm.empty()) { return EXIT_FAILURE; }

  // Now read the standard options
  logLevel = FW::Options::readLogLevel(vm);
  //auto nEvents  = FW::Options::readSequencerConfig(vm).events;

  // The geometry, material and decoration
  auto geometry          = FW::Geometry::build(vm, detector);
  // geometry is a pair of (tgeoTrackingGeometry, tgeoContextDecorators)
  auto tGeometry         = geometry.first;
  //std::vector<std::shared_ptr<FW::IContextDecorator> > contextDecorators = geometry.second; 
  contextDecorators = geometry.second; 

  auto magneticField = FW::Options::readBField(vm);

  // The detectors
  // "MVTX" and "Silicon"
  read_strings subDetectors = vm["geo-subdetectors"].as<read_strings>();

  size_t ievt = 0;
  size_t ialg = 0;

  // Setup the event and algorithm context
  FW::WhiteBoard eventStore(Acts::getDefaultLogger("EventStore#" + std::to_string(ievt), 
						   logLevel));
  
  // The geometry context
  FW::AlgorithmContext context(ialg, ievt, eventStore);

  // Make a fit configuration 
  fitCfg.fit = FW::TrkrClusterFittingAlgorithm::makeFitterFunction(tGeometry, 
  						    magneticField,
  						    logLevel);
  m_fitCfgOptions = new FitCfgOptions(fitCfg, 
				      context.calibContext,
				      context.geoContext,
				      context.magFieldContext);
  

  // this is not executed because contextDecorators has size 0
  /// Decorate the context
  for (auto& cdr : contextDecorators) {
    if (cdr->decorate(context) != FW::ProcessCode::SUCCESS)
      throw std::runtime_error("Failed to decorate event context");
  }

  geo_ctxt = context.geoContext;

  // tGeometry is a TrackingGeometry pointer, acquired in the first 20 lines/ of this file
  auto vol = tGeometry->highestTrackingVolume();
  // vol is a TrackingVolume pointer

  // Get the confined volumes in the highest tracking volume
  // confinedVolumes is a shared_ptr<TrackingVolumeArray>
  auto confinedVolumes = vol->confinedVolumes();

  // volumeVector is a std::vector<TrackingVolumePtrs>
  auto volumeVector = confinedVolumes->arrayObjects();

  // The first entry is the MVTX
  //=====================
  auto mvtxVolumes = volumeVector.at(0)->confinedVolumes();
  // mvtxVolumes is a shared_ptr<TrackingVolumeArray>
  // Now get the individual TrackingVolumePtrs corresponding to each MVTX volume
  auto mvtxBarrel = mvtxVolumes->arrayObjects().at(1);

  // Now get the LayerArrays corresponding to each volume
  auto mvtxBarrelLayerArray = mvtxBarrel->confinedLayers();  // the barrel is all we care about

  // This is a BinnedArray<LayerPtr>, but we care only about the sensitive surfaces
  auto mvtxLayerVector1 = mvtxBarrelLayerArray->arrayObjects();  // the barrel is all we care about

  // mvtxLayerVector has size 3, but only index 1 has any surfaces since the clayersplits option was removed
  // the layer number has to be deduced from the radius
  for(unsigned int i=0;i<mvtxLayerVector1.size(); ++i)
    {

      // Get the Acts::SurfaceArray from each MVTX LayerPtr being iterated over
      auto surfaceArray = mvtxLayerVector1.at(i)->surfaceArray();
      if(surfaceArray == NULL)
	continue;

      double ref_rad[3] = {2.556, 3.359, 4.134};
      
      // surfaceVector is an Acts::SurfaceVector, vector of surfaces
      //std::vector<const Surface*>
      auto surfaceVector = surfaceArray->surfaces();
      for(unsigned int j=0; j<surfaceVector.size(); j++){

	auto surf = surfaceVector.at(j)->getSharedPtr();

	auto vec3d =  surf->center(context.geoContext);
	std::vector<double> world_center = { vec3d(0)/10.0, vec3d(1)/10.0, vec3d(2)/10.0 };  // convert from mm to cm
	double layer_rad = sqrt(pow(world_center[0],2) + pow(world_center[1],2));
	unsigned int layer = 0;
	for(unsigned int i = 0;i<3;++i)
	  {
	    if( fabs(layer_rad - ref_rad[i]) < 0.1)
	      layer = i;
	  }

	TrkrDefs::hitsetkey hitsetkey = GetMvtxHitSetKeyFromCoords(layer, world_center);

	// Add this surface to the map
	std::pair<TrkrDefs::hitsetkey, std::shared_ptr<const Acts::Surface>> tmp = make_pair(hitsetkey, surf);

	_cluster_surface_map_silicon.insert(tmp);


	if(_verbosity > 0)
	  {
	    unsigned int stave = MvtxDefs::getStaveId(hitsetkey);
	    unsigned int chip = MvtxDefs::getChipId(hitsetkey);

	    // check it is in there
	    std::cout << "Layer radius " << layer_rad << " Layer " << layer << " stave " << stave << " chip " << chip 
		      << " recover surface from _cluster_surface_map_silicon " << std::endl;
	    std::map<TrkrDefs::hitsetkey, std::shared_ptr<const Acts::Surface>>::iterator surf_iter = _cluster_surface_map_silicon.find(hitsetkey);
	    std::cout << " surface type " << surf_iter->second->type() << std::endl;
	    surf_iter->second->toStream(geo_ctxt,std::cout);	
	    auto assoc_layer = surf->associatedLayer();
	    std::cout << std::endl << " Layer type "  << assoc_layer->layerType() << std::endl;

	    auto assoc_det_element = surf->associatedDetectorElement();
	    if(assoc_det_element != nullptr)
	      {
		std::cout << " Associated detElement has non-null pointer " << assoc_det_element << std::endl;
		std::cout << std::endl << " Associated detElement found, thickness = "  << assoc_det_element->thickness() << std::endl;
	      }
	    else
	      std::cout << std::endl << " Associated detElement is nullptr " << std::endl;
	    
	  }

      }
    }

  // end of MVTX
  //=================================

  // INTT only has one volume, so there is not an added iterator like for the  MVTX
  //========================================================
  auto inttVolume = volumeVector.at(1);
  auto inttLayerArray = inttVolume->confinedLayers();

  auto inttLayerVector = inttLayerArray->arrayObjects();
  // inttLayerVector is a std::vector<LayerPtr>
  for(unsigned int i=0; i<inttLayerVector.size(); i++){
  
    //auto vol = inttLayerVector.at(i)->trackingVolume();

    // Get the ACTS::SurfaceArray from each INTT LayerPtr being iterated over
    auto surfaceArray = inttLayerVector.at(i)->surfaceArray();
    if(surfaceArray == NULL)
      continue;
    
    // surfaceVector is an Acts::SurfaceVector, vector of surfaces
    auto surfaceVector = surfaceArray->surfaces();
    for(unsigned int j=0; j<surfaceVector.size(); j++){
      
      auto surf = surfaceVector.at(j)->getSharedPtr();
   
      auto vec3d =  surf->center(context.geoContext);

      double ref_rad[4] = {8.987, 9.545, 10.835, 11.361};
      
      std::vector<double> world_center = { vec3d(0)/10.0, vec3d(1)/10.0, vec3d(2)/10.0 };  // convert from mm to cm
      // The Acts geometry builder combines layers 4 and 5 together, and layers 6 and 7 together. We need to use the radius to figure 
      // out which layer to use to get the layergeom
      double layer_rad = sqrt(pow(world_center[0],2) + pow(world_center[1],2));
      
	unsigned int layer = 0;
	for(unsigned int i = 0;i<4;++i)
	  {
	    if( fabs(layer_rad - ref_rad[i]) < 0.1)
	      layer = i + 3;
	  }

	TrkrDefs::hitsetkey hitsetkey = GetInttHitSetKeyFromCoords(layer, world_center);

	// Add this surface to the map
	std::pair<TrkrDefs::hitsetkey, std::shared_ptr<const Acts::Surface>> tmp = make_pair(hitsetkey, surf);
	_cluster_surface_map_silicon.insert(tmp);

	if(_verbosity > 0)
	  {
	    // check it is in there
	    unsigned int ladderPhi = InttDefs::getLadderPhiId(hitsetkey);
	    unsigned int ladderZ = InttDefs::getLadderZId(hitsetkey);
 
	    std::cout <<"Layer radius " << layer_rad << " layer " << layer << " ladderPhi " << ladderPhi << " ladderZ " << ladderZ 
		      << " recover surface from _cluster_surface_map_silicon " << std::endl;
	    std::map<TrkrDefs::hitsetkey, std::shared_ptr<const Acts::Surface>>::iterator surf_iter = _cluster_surface_map_silicon.find(hitsetkey);
	    std::cout << " surface type " << surf->type() << std::endl;
	    surf_iter->second->toStream(geo_ctxt,std::cout);
	    auto assoc_layer = surf->associatedLayer();
	    std::cout << std::endl << " Layer type "  << assoc_layer->layerType() << std::endl;

	    auto assoc_det_element = surf->associatedDetectorElement();
	    if(assoc_det_element != nullptr)
	      {
		std::cout << " Associated detElement has non-null pointer " << assoc_det_element << std::endl;
		std::cout << std::endl << " Associated detElement found, thickness = "  << assoc_det_element->thickness() << std::endl;
	      }
	    else
	      std::cout << std::endl << " Associated detElement is nullptr " << std::endl;
	  }
    }
  }

  // end of INTT
  //=================================
  
  return 0;
}

TrkrDefs::hitsetkey MakeActsGeometry::GetMvtxHitSetKeyFromCoords(unsigned int layer, std::vector<double> &world)
{
  // Look up the MVTX sensor index values from the world position of the surface center

  CylinderGeom_Mvtx *layergeom = dynamic_cast<CylinderGeom_Mvtx *>(_geom_container_mvtx->GetLayerGeom(layer));
  if(!layergeom)
    {
      std::cout << PHWHERE << "Did not get layergeom for layer " << layer  << std::endl;
    }

  unsigned int stave = 0;
  unsigned int chip = 0;
  layergeom->get_sensor_indices_from_world_coords(world, stave, chip);

  double check_pos[3] = {0,0,0};
  layergeom->find_sensor_center(stave, 0, 0, chip, check_pos);

  TrkrDefs::hitsetkey mvtx_hitsetkey = MvtxDefs::genHitSetKey(layer, stave, chip);
  
  return mvtx_hitsetkey;
}

TrkrDefs::hitsetkey MakeActsGeometry::GetInttHitSetKeyFromCoords(unsigned int layer, std::vector<double> &world)
{
  // Look up the INTT sensor index values from the world position of the surface center

  CylinderGeomIntt *layergeom = dynamic_cast<CylinderGeomIntt *>(_geom_container_intt->GetLayerGeom(layer));
  if(!layergeom)
    {
      std::cout << PHWHERE << "Did not get layergeom for layer " << layer  << std::endl;
    }

  double location[3] = {world[0], world[1], world[2]};
  int segment_z_bin = 0;
  int segment_phi_bin = 0;
  layergeom->find_indices_from_segment_center(segment_z_bin, segment_phi_bin, location);

  double check_pos[3] = {0,0,0};
  layergeom->find_segment_center(segment_z_bin, segment_phi_bin, check_pos);

  TrkrDefs::hitsetkey intt_hitsetkey = InttDefs::genHitSetKey(layer, segment_z_bin, segment_phi_bin);
  
  return intt_hitsetkey;
}


void MakeActsGeometry::MakeTGeoNodeMap(PHCompositeNode *topNode)
{
  _geomanager = PHGeomUtility::GetTGeoManager(topNode);
  if(!_geomanager )
    {
      cout << PHWHERE << " Did not find TGeoManager, quit! " << endl;
      return;
    }
  TGeoVolume *topVol = _geomanager->GetTopVolume();
  TObjArray *nodeArray = topVol->GetNodes();

  TIter iObj(nodeArray); 
  while(TObject *obj = iObj())
    {
      TGeoNode *node = dynamic_cast<TGeoNode*>(obj);
      std::string node_str = node->GetName();

      std::string mvtx("av_1");
      std::string intt("ladder");
      std::string intt_ext("ladderext");

      if ( node_str.compare(0, mvtx.length(), mvtx) == 0 )       // is it in the MVTX?
	{
	  if(_verbosity > 100)  cout << " node " << node->GetName() << " is in the MVTX" << endl;
	  getMvtxKeyFromNode(node);
	}
      else if ( node_str.compare(0, intt.length(), intt) == 0 ) 	      // is it in the INTT?
	{
	  // We do not want the "ladderext" nodes
	  if ( node_str.compare(0, intt_ext.length(), intt_ext) == 0 ) 
	    continue;
	  
	  if(_verbosity > 100) cout << " node " << node->GetName() << " is in the INTT" << endl;	  
	  getInttKeyFromNode(node);
	}
      else
	continue;

      bool print_sensor_paths = false;  // normally false
      if(print_sensor_paths)
	{
	  // Descends the node tree to find the active silicon nodes - used for information only
	  cout<< " Top Node is " << node->GetName() << " volume name is " << node->GetVolume()->GetName()  << endl;
	  cout << " Top Node mother volume name is " << node->GetMotherVolume()->GetName() << endl;
	  isActive(node);
	}
    }
}

void  MakeActsGeometry::getInttKeyFromNode(TGeoNode *gnode)
{
  int layer = -1;           // sPHENIX layer number
  int itype = -1;           // specifies inner (0) or outer (1) sensor
  int ladder_phi = -1;  // copy number of ladder in phi
  int zposneg = -1;                // specifies positive (1) or negative (0) z
  int ladder_z = -1;      // 0-3, from most negative z to most positive
  
  std::string s = gnode->GetName();
  std::string delimiter = "_";
  std::string posz("posz");
  std::string negz("negz");
  
  size_t pos = 0;
  std::string token;

  int counter = 0;
  while ((pos = s.find(delimiter)) != std::string::npos) {
    token = s.substr(0, pos);
    //std::cout << token << std::endl;
    s.erase(0, pos + delimiter.length());
    if(counter == 1) 
      layer = std::atoi(token.c_str()) + 3;
    if(counter == 2)
      itype = std::atoi(token.c_str());
    if(counter == 3)
      {
	ladder_phi = std::atoi(token.c_str());
	if( s.compare(0, negz.length(), negz) ==0 ) zposneg = 0; 
	if( s.compare(0, posz.length(), posz) ==0 ) zposneg = 1; 
      }	
    counter ++;
  }

  ladder_z = itype  + zposneg*2;  

  // The active sensor is a daughter of gnode
  int ndaught = gnode->GetNdaughters();
  if(ndaught == 0)
    {
      cout << PHWHERE << "OOPS: Did not find INTT sensor! Quit." << endl;
      exit(1);
    }

  std::string intt_refactive("siactive");  
  TGeoNode *sensor_node = 0;
  for(int i=0; i<ndaught; ++i)
    {
      std::string node_str = gnode->GetDaughter(i)->GetName();

      if (node_str.compare(0, intt_refactive.length(), intt_refactive) == 0)
	{
	  sensor_node = gnode->GetDaughter(i);      
	  break;      
	}
    } 
 
  // unique key identifying this sensor
  TrkrDefs::hitsetkey node_key = InttDefs::genHitSetKey(layer, ladder_z, ladder_phi);

  std::pair<TrkrDefs::hitsetkey, TGeoNode*> tmp = make_pair(node_key, sensor_node);
  _cluster_node_map.insert(tmp);

  if(_verbosity > 1)    
    std::cout << " INTT layer " << layer << " ladder_phi " << ladder_phi << " itype " << itype << " zposneg " << zposneg << " ladder_z " << ladder_z << " name " << sensor_node->GetName() << std::endl;
  
  return;
}

void MakeActsGeometry::getMvtxKeyFromNode(TGeoNode *gnode)
{
  int counter = 0;
  int impr = -1;   // stave number, 1-48 in TGeo
  int layer = -1;
  int stave = -1;  // derived from impr
  int chip = -1;   // 9 chips per stave

  std::string s = gnode->GetName();
  std::string delimiter = "_";
  
  size_t pos = 0;
  std::string token;

  while ((pos = s.find(delimiter)) != std::string::npos) {
    token = s.substr(0, pos);
    //std::cout << token << std::endl;
    s.erase(0, pos + delimiter.length());
    if(counter == 3) 
      impr = std::atoi(token.c_str());
 
    counter ++;
  }

  // extract layer and stave info from impr
  // int staves_in_layer[3] = {12, 16, 20}; 
  // note - impr stave count starts from 1, not 0, but TrkrCluster counting starts from 0, so we reduce it by 1 here
  impr -= 1;
 
 if(impr < 12)
    {
      layer = 0;
      stave = impr;
    }
  else if(impr > 11 && impr < 28)
    {
      layer = 1;
      stave = impr - 12;
    }
  else
    {
      layer = 2;
      stave = impr - 28;
    }

  // Now descend node tree to find chip ID's - there are multiple chips per stave
  TGeoNode *module_node = gnode->GetDaughter(0);    
  int mnd = module_node->GetNdaughters();
  std::string mvtx_chip("MVTXChip");
  for(int i=0; i<mnd; ++i)
    {
      std::string dstr = module_node->GetDaughter(i)->GetName();
      if (dstr.compare(0, mvtx_chip.length(), mvtx_chip) == 0)
	{
	  if(_verbosity > 1) 
	    cout << "Found MVTX layer " << layer << " stave " << stave << " chip  " << i << " with node name " <<  module_node->GetDaughter(i)->GetName() << endl;

	  // Make key for this chip
	  TrkrDefs::hitsetkey node_key = MvtxDefs::genHitSetKey(layer, stave, i);

	  // add sensor node to map
	  TGeoNode *sensor_node = module_node->GetDaughter(i)->GetDaughter(0);
	  std::pair<TrkrDefs::hitsetkey, TGeoNode*> tmp = make_pair(node_key, sensor_node);
	  _cluster_node_map.insert(tmp);
	  
	  if(_verbosity > 1)    
	    std::cout << " MVTX layer " << layer << " stave " << stave << " chip " << chip << " name " << sensor_node->GetName() << std::endl;
	}
    }
  
  return;
}

void MakeActsGeometry::isActive(TGeoNode *gnode)
{
  // Not used in analysis: diagnostic, for looking at the node tree only.
  // Recursively searches gnode for silicon sensors, prints out heirarchy

  std::string node_str = gnode->GetName();

  std::string intt_refactive("siactive");
  std::string mvtx_refactive("MVTXSensor");

  if (node_str.compare(0, intt_refactive.length(), intt_refactive) == 0)
    {
      cout << "          ******* Found INTT active volume,  node is " << gnode->GetName() 
	   << " volume name is "   << gnode->GetVolume()->GetName() << endl;

      //const TGeoMatrix* tgMatrix = gnode->GetMatrix();
      //tgMatrix->Print();

      return;
    }
  else if (node_str.compare(0, mvtx_refactive.length(), mvtx_refactive) == 0)
    {
      cout << "          ******* Found MVTX active volume,  node is " << gnode->GetName() 
	   << " volume name is " << gnode->GetVolume()->GetName() << endl;

      //const TGeoMatrix* tgMatrix = gnode->GetMatrix();
      //tgMatrix->Print();

      return;
    }

  int ndaught = gnode->GetNdaughters();
  if(ndaught == 0)
    {
      cout << "     No further daughters" << endl;
    }

  for(int i=0; i<ndaught; ++i)
    {
      cout << "     " << gnode->GetVolume()->GetName() << "  daughter " << i 
	   << " has name " << gnode->GetDaughter(i)->GetVolume()->GetName() << endl;
      isActive(gnode->GetDaughter(i));      
    }
}

MakeActsGeometry::~MakeActsGeometry()
{

}


int MakeActsGeometry::CreateNodes(PHCompositeNode* topNode)
{
  
  return Fun4AllReturnCodes::EVENT_OK;
}

/*
 * GetNodes():
 *  Get all the all the required nodes off the node tree
 */

int MakeActsGeometry::GetNodes(PHCompositeNode* topNode)
{
  _geomanager = PHGeomUtility::GetTGeoManager(topNode);
  if(!_geomanager )
    {
      cout << PHWHERE << " Did not find TGeoManager, quit! " << endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }
  
  _geom_container_mvtx = findNode::getClass<
    PHG4CylinderGeomContainer>(topNode, "CYLINDERGEOM_MVTX");
  if (!_geom_container_mvtx)
  {
    cout << PHWHERE << " CYLINDERGEOM_MVTX  node not found on node tree"
         << endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  _geom_container_tpc =
    findNode::getClass<PHG4CylinderCellGeomContainer>(topNode, "CYLINDERCELLGEOM_SVTX");
  if (!_geom_container_tpc)
    {
      std::cout << PHWHERE << "ERROR: Can't find node CYLINDERCELLGEOM_SVTX" << std::endl;
      return Fun4AllReturnCodes::ABORTRUN;
    }


  _geom_container_intt = findNode::getClass<
    PHG4CylinderGeomContainer>(topNode, "CYLINDERGEOM_INTT");
  if (!_geom_container_intt)
    {
      cout << PHWHERE << " CYLINDERGEOM_INTT  node not found on node tree"
	   << endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }
  
  return Fun4AllReturnCodes::EVENT_OK;
}



