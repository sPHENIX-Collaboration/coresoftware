/*!
 *  \file		PHActsTrkFitter.C
 *  \brief		Refit SvtxTracks with PHActs.
 *  \details	Refit SvtxTracks with PHActs.
 *  \author	        Tony Frawley <afrawley@fsu.edu>
 */

#include "PHActsTrkFitter.h"

#include <trackbase/TrkrCluster.h>                  // for TrkrCluster
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrDefs.h>

#include <intt/CylinderGeomIntt.h>
#include <intt/InttDefs.h>

#include <mvtx/CylinderGeom_Mvtx.h>
#include <mvtx/MvtxDefs.h>

#include <tpc/TpcDefs.h>

#include <trackbase_historic/SvtxTrack.h>
#include <trackbase_historic/SvtxTrackMap.h>

#include <g4detectors/PHG4CylinderGeom.h>           // for PHG4CylinderGeom
#include <g4detectors/PHG4CylinderGeomContainer.h>
#include <g4detectors/PHG4CylinderCellGeom.h>
#include <g4detectors/PHG4CylinderCellGeomContainer.h>

#include <phgeom/PHGeomUtility.h>


#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/getClass.h>
#include <phool/phool.h>


// Acts classes
#include <Acts/Geometry/GeometryContext.hpp>
#include <Acts/Geometry/TrackingGeometry.hpp>
#include <Acts/Geometry/TrackingVolume.hpp>

#include <Acts/Plugins/Digitization/DigitizationCell.hpp>
#include <Acts/Plugins/Digitization/PlanarModuleCluster.hpp>

#include <Acts/Surfaces/Surface.hpp>
#include <Acts/Surfaces/PlaneSurface.hpp>

#include <Acts/Utilities/Definitions.hpp>
#include <Acts/Utilities/BinnedArray.hpp>                       // for Binne...
#include <Acts/Utilities/Logger.hpp>                            // for getDe...

// This needs to stay after these other Acts includes
#include "TrkrClusterSourceLink.hpp"


#include <ACTFW/Detector/IBaseDetector.hpp>

#include <ACTFW/Framework/AlgorithmContext.hpp>
#include <ACTFW/Framework/IContextDecorator.hpp>
#include <ACTFW/Framework/WhiteBoard.hpp>

#include <ACTFW/Geometry/CommonGeometry.hpp>

#include <ACTFW/Options/CommonOptions.hpp>

#include <ACTFW/Plugins/Obj/ObjWriterOptions.hpp>

#include <ACTFW/TGeoDetector/TGeoDetector.hpp>

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
PHActsTrkFitter::PHActsTrkFitter(const string& name)
  : PHTrackFitting(name)
  , _geom_container_mvtx(nullptr)
  , _geom_container_intt(nullptr)
  , _geom_container_tpc(nullptr)
  , _trackmap(nullptr)
  , _clustermap(nullptr)
  , _geomanager(nullptr)
  , MinSurfZ(0.0)
  , MaxSurfZ(100.0)
  , NSurfZ(10)
  , NSurfPhi(10)
{
  Verbosity(0);

  // These are arbitrary subdivisions, and may change
  SurfStepZ = (MaxSurfZ - MinSurfZ) / (double) NSurfZ;
  ModuleStepPhi = 2.0 * M_PI / 12.0;
  ModulePhiStart = - M_PI;
  SurfStepPhi = 2.0 * M_PI / (double) (NSurfPhi * NTpcModulesPerLayer);

  _event = 0;
}

int PHActsTrkFitter::Setup(PHCompositeNode *topNode)
{
  GetNodes(topNode);

  // run Acts layer builder
  BuildSiliconLayers();

  // create a map of sensor TGeoNode pointers using the TrkrDefs:: hitsetkey as the key  
  MakeTGeoNodeMap(topNode);

  // TPC continuous readout geometry does not exist within ACTS, so we build our own surfaces
  BuildTpcSurfaceMap();
  
  return Fun4AllReturnCodes::EVENT_OK;
}

void PHActsTrkFitter::BuildTpcSurfaceMap()
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
		      
		      // The normsal vector is a unit vector pointing to the center of the beam line at z = z_center
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
void PHActsTrkFitter::BuildSiliconLayers()
{
  // define int argc and char* argv to provide options to processGeometry

  /*
  int argc = 27;
    char *arg[27];
    std::string argstr[27]{"-n1", "-l0", "--geo-tgeo-filename=none", "--geo-tgeo-worldvolume=\"World\"", "--geo-subdetectors", "MVTX", "Silicon", "--geo-tgeo-nlayers=0", "0", "--geo-tgeo-clayers=1",  "1", "--geo-tgeo-players=0", "0", "--geo-tgeo-clayernames", "MVTX", "siactive", "--geo-tgeo-cmodulenames", "MVTXSensor",  "siactive",  "--geo-tgeo-cmoduleaxes", "xzy", "yzx", "--geo-tgeo-clayersplits", "5.",  "10.",  "--output-obj", "true"};
  */

  // Can hard code geometry options since the TGeo options are fixed by our detector design
  int argc = 24;
  char *arg[24];
  std::string argstr[24]{"-n1", "-l0", "--geo-tgeo-filename=none", "--geo-tgeo-worldvolume=\"World\"", "--geo-subdetectors", "MVTX", "Silicon", "--geo-tgeo-nlayers=0", "0", "--geo-tgeo-clayers=1",  "1", "--geo-tgeo-players=0", "0", "--geo-tgeo-clayernames", "MVTX", "siactive", "--geo-tgeo-cmodulenames", "MVTXSensor",  "siactive",  "--geo-tgeo-cmoduleaxes", "xzy", "yzx",  "--output-obj", "true"};


  for(int i = 0;i < argc; ++i)
    {
      // need a copy, since .c_str() returns a const char * and process geometry will not take a const
      arg[i] = strdup(argstr[i].c_str());
    }

  /*
    acts-framework:
    (compiled binary = ACTFWTGeoGeometryExample)
    (basically just feeds options to GeometryExampleBase::ProcessGeometry)
        (- calls GeometryExampleBase = acts-framework/Examples/Common/src/GeometryExampleBase(.cpp)::ProcessGeometry() with arguments and the instance of TGeoDetector)
           (sets up options, adds options to detector)
    **** All of the above replaced by the method here called MakeActsGeometry()
           - calls Geometry::build()  =  in acts-framework/Examples/Common/include/ACTFW/Geometry/CommonGeometry.hpp (= namespace Geometry) 
              material decoration
              returns TGeoDetector::finalize =  std::pair<FW::IBaseDetector::TrackingGeometryPtr, ContextDecorators>
              TrackingGeometry = acts-core/Core/include/Acts/Geometry/TrackingGeometry.hpp 
              - calls TGeoDetecrtor::finalize() = cts-framework/Detectors/TGeoDetector/include/ACTFW/TGeoDetector/TGeoDetector(.hpp)::finalize (inherits from FW::IBaseDetector)
                    creates DetectorStore
                    returns   std::pair<FW::IBaseDetector::TrackingGeometryPtr, ContextDecorators>     (TrackingGeometryPtr = std::shared_ptr<const Acts::TrackingGeometry>)
                  - calls TGeo::buildTGeoDetector()  = acts-framework/Detectors/TGeoDetector/include/ACTFW/TGeoDetector/BuildTGeoDetector.hpp ( = namespace TGeo)
                      initializes TGeoManager from file
   acts-core from here on:
                      - calls TGeoLayerBuilder::buildLayers = acts-core/Plugins/TGeo/include/Acts/Plugins/TGeo/TGeoLayerBuilder(.hpp)::buildLayers  (inherits from Acts::ILayerBuilder)
  */

  // We replicate the relevant functionality of  acts-framework/Examples/Common/src/GeometryExampleBase::ProcessGeometry() in MakeActsGeometry()
  // so we get access to the results. The layer builder magically gets the TGeoManager
  TGeoDetector detector;
  MakeActsGeometry(argc, arg, detector);

}

int PHActsTrkFitter::MakeActsGeometry(int argc, char* argv[], FW::IBaseDetector& detector)
{
  // setup and parse options
  auto desc = FW::Options::makeDefaultOptions();
  FW::Options::addSequencerOptions(desc);
  FW::Options::addGeometryOptions(desc);
  FW::Options::addMaterialOptions(desc);
  FW::Options::addObjWriterOptions(desc);
  //  FW::Options::addCsvWriterOptions(desc);
  FW::Options::addOutputOptions(desc);

  // Add specific options for this geometry
  detector.addOptions(desc);
  auto vm = FW::Options::parse(desc, argc, argv);
  if (vm.empty()) { return EXIT_FAILURE; }

  // Now read the standard options
  auto logLevel = FW::Options::readLogLevel(vm);
  //auto nEvents  = FW::Options::readSequencerConfig(vm).events;

  // The geometry, material and decoration
  auto geometry          = FW::Geometry::build(vm, detector);
  // geometry is a pair of (tgeoTrackingGeometry, tgeoContextDecorators)
  auto tGeometry         = geometry.first;
  //std::vector<std::shared_ptr<FW::IContextDecorator> > contextDecorators = geometry.second; 
  contextDecorators = geometry.second; 

  // The detectors
  // "MVTX" and "Silicon"
  read_strings subDetectors = vm["geo-subdetectors"].as<read_strings>();

  size_t ievt = 0;
  size_t ialg = 0;

  // Setup the event and algorithm context
  FW::WhiteBoard eventStore(
			    Acts::getDefaultLogger("EventStore#" + std::to_string(ievt), logLevel));
  
  // The geometry context
  FW::AlgorithmContext context(ialg, ievt, eventStore);

  // this is not executed because contextDecorators has size 0
  /// Decorate the context
  for (auto& cdr : contextDecorators) {
    if (cdr->decorate(context) != FW::ProcessCode::SUCCESS)
      throw std::runtime_error("Failed to decorate event context");
  }

  geo_ctxt = context.geoContext;

  /*
  // This is not executed because size is 0
  std::string geoContextStr = "";
  if (contextDecorators.size() > 0) {
    // We need indeed a context object
    if (nEvents > 1) { geoContextStr = "_geoContext" + std::to_string(ievt); }
  }
  */

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

  // Thiis is a BinnedArray<LayerPtr>, but we care only about the sensitive surfaces
  auto mvtxLayerVector1 = mvtxBarrelLayerArray->arrayObjects();  // the barrel is all we care about

  // mvtxLayerVector has size 3, but only index 1 has any surfaces since the clayersplits option was removed
  // the layer number has to be deduced from the radius
  for(unsigned int i=0;i<mvtxLayerVector1.size(); ++i)
    {

      // Get the ACTS::SurfaceArray from each MVTX LayerPtr being iterated over
      auto surfaceArray = mvtxLayerVector1.at(i)->surfaceArray();
      if(surfaceArray == NULL)
	continue;

      double ref_rad[3] = {2.556, 3.359, 4.134};
      
      // surfaceVector is an ACTS::SurfaceVector, vector of surfaces
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

	_cluster_surface_map.insert(tmp);

	if(Verbosity() > 0)
	  {
	    unsigned int stave = MvtxDefs::getStaveId(hitsetkey);
	    unsigned int chip = MvtxDefs::getChipId(hitsetkey);

	    // check it is in there
	    std::cout << "Layer radius " << layer_rad << " Layer " << layer << " stave " << stave << " chip " << chip 
		      << " recover surface from _cluster_surface_map " << std::endl;
	    std::map<TrkrDefs::hitsetkey, std::shared_ptr<const Acts::Surface>>::iterator surf_iter = _cluster_surface_map.find(hitsetkey);
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
    
    // surfaceVector is an ACTS::SurfaceVector, vector of surfaces
    auto surfaceVector = surfaceArray->surfaces();
    for(unsigned int j=0; j<surfaceVector.size(); j++){
      
      auto surf =  surfaceVector.at(j)->getSharedPtr();
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
	_cluster_surface_map.insert(tmp);

	if(Verbosity() > 0)
	  {
	    // check it is in there
	    unsigned int ladderPhi = InttDefs::getLadderPhiId(hitsetkey);
	    unsigned int ladderZ = InttDefs::getLadderZId(hitsetkey);
 
	    std::cout <<"Layer radius " << layer_rad << " layer " << layer << " ladderPhi " << ladderPhi << " ladderZ " << ladderZ 
		      << " recover surface from _cluster_surface_map " << std::endl;
	    std::map<TrkrDefs::hitsetkey, std::shared_ptr<const Acts::Surface>>::iterator surf_iter = _cluster_surface_map.find(hitsetkey);
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

TrkrDefs::hitsetkey PHActsTrkFitter::GetMvtxHitSetKeyFromCoords(unsigned int layer, std::vector<double> &world)
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

TrkrDefs::hitsetkey PHActsTrkFitter::GetInttHitSetKeyFromCoords(unsigned int layer, std::vector<double> &world)
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


void PHActsTrkFitter::MakeTGeoNodeMap(PHCompositeNode *topNode)
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
	  if(Verbosity() > 100)  cout << " node " << node->GetName() << " is in the MVTX" << endl;
	  getMvtxKeyFromNode(node);
	}
      else if ( node_str.compare(0, intt.length(), intt) == 0 ) 	      // is it in the INTT?
	{
	  // We do not want the "ladderext" nodes
	  if ( node_str.compare(0, intt_ext.length(), intt_ext) == 0 ) 
	    continue;
	  
	  if(Verbosity() > 100) cout << " node " << node->GetName() << " is in the INTT" << endl;	  
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

void  PHActsTrkFitter::getInttKeyFromNode(TGeoNode *gnode)
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

  if(Verbosity() > 1)    
    std::cout << " INTT layer " << layer << " ladder_phi " << ladder_phi << " itype " << itype << " zposneg " << zposneg << " ladder_z " << ladder_z << " name " << sensor_node->GetName() << std::endl;
  
  return;
}

void PHActsTrkFitter::getMvtxKeyFromNode(TGeoNode *gnode)
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
	  if(Verbosity() > 1) 
	    cout << "Found MVTX layer " << layer << " stave " << stave << " chip  " << i << " with node name " <<  module_node->GetDaughter(i)->GetName() << endl;

	  // Make key for this chip
	  TrkrDefs::hitsetkey node_key = MvtxDefs::genHitSetKey(layer, stave, i);

	  // add sensor node to map
	  TGeoNode *sensor_node = module_node->GetDaughter(i)->GetDaughter(0);
	  std::pair<TrkrDefs::hitsetkey, TGeoNode*> tmp = make_pair(node_key, sensor_node);
	  _cluster_node_map.insert(tmp);
	  
	  if(Verbosity() > 1)    
	    std::cout << " MVTX layer " << layer << " stave " << stave << " chip " << chip << " name " << sensor_node->GetName() << std::endl;
	}
    }
  
  return;
}

void PHActsTrkFitter::isActive(TGeoNode *gnode)
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

int PHActsTrkFitter::Process()
{
  _event++;

  if (Verbosity() > 1)
    {
      std::cout << PHWHERE << "Events processed: " << _event << std::endl;
            cout << "Start PHActsTrkfitter::process_event" << endl;
    }

  TrkrDefs::hitsetkey hsetkey;

  TrkrClusterContainer::ConstRange clusrange = _clustermap->getClusters();
  for(TrkrClusterContainer::ConstIterator clusiter = clusrange.first; clusiter != clusrange.second; ++clusiter)
    {
      TrkrCluster *cluster = clusiter->second;
      TrkrDefs::cluskey cluskey = clusiter->first;

      // get the cluster parameters in global coordinates
      float x = cluster->getPosition(0);
      float y = cluster->getPosition(1);
      float z = cluster->getPosition(2);
      double radius = sqrt(x*x+y*y);

      // In local coords the covariances are in the  r*phi vs z frame
      // They have been rotated into global coordinates in TrkrCluster
      TMatrixD world_err(3,3);
      for(int i=0; i < 3; ++i)
	for(int j =0; j<3; j++)
	  {
	    world_err[i][j] = cluster->getError(i,j);
	  }

      // extract detector element identifiers from cluskey so we can access the Surface and TGeoNode
      unsigned int layer = TrkrDefs::getLayer(cluskey);
      if(Verbosity() > 0) cout << std::endl << " layer " << layer << endl;

      TVector3 world(x,y,z);
      TVector3 local(0,0,0);
      TMatrixD local_err(3, 3);

      double local_2D[2] = {0};

      unsigned int trkrid = TrkrDefs::getTrkrId(cluskey);  // 0 for MVTX, 1 for INTT, 2 for TPC
      TGeoNode *sensor_node;
      std::shared_ptr<const Acts::Surface> surf;

      // Getting the hitsetkey from the cluskey is detector specific
      if(trkrid == TrkrDefs::mvtxId)   	  // MVTX
	{
	  unsigned int staveid = MvtxDefs::getStaveId(cluskey);
	  unsigned int chipid = MvtxDefs::getChipId(cluskey);
	  if(Verbosity() > 0)   cout << "   MVTX cluster with staveid " << staveid << " chipid " << chipid << endl; 

	  // make key for this sensor
	  hsetkey = MvtxDefs::genHitSetKey(layer, staveid, chipid);

	  // get the TGeoNode for it
	  std::map<TrkrDefs::hitsetkey, TGeoNode*>::iterator it;
	  it = _cluster_node_map.find(hsetkey);
	  if(it != _cluster_node_map.end())
	    {
	      sensor_node = it->second;
	      if(Verbosity() > 0) 
		cout << "    Found in _cluster_node_map: layer " << layer << " staveid " << staveid << " chipid " << chipid 
		     <<  " node " << sensor_node->GetName() << endl;
	    }
	  else
	    {
	      cout << PHWHERE << "No entry in TGeo map for cluster: layer " << layer << " staveid " << staveid << " chipid " << chipid  << " - should be impossible!" << endl;
	      return Fun4AllReturnCodes::ABORTEVENT;
	    }
	 
	  // Find Acts surface corresponding to it
	  std::map<TrkrDefs::hitsetkey, std::shared_ptr<const Acts::Surface>>::iterator surf_iter;
	  surf_iter = _cluster_surface_map.find(hsetkey);  
	  if(surf_iter != _cluster_surface_map.end())
	    {	      
	      //TrkrDefs::hitsetkey found_key = surf_iter->first;
	      surf = surf_iter->second->getSharedPtr();
	      //cout<< "Got surface pair " << surf->name() << " surface type " << surf->type() << std::endl;
	      //surf->toStream(geo_ctxt, std::cout);  // this causes a segfault for mvtx and intt
	    }
	  else
	    {
	      cout << PHWHERE << "Failed to find associated surface element - should be impossible!" << endl;
	      return Fun4AllReturnCodes::ABORTEVENT;
	    }

	  CylinderGeom_Mvtx *layergeom = dynamic_cast<CylinderGeom_Mvtx *>(_geom_container_mvtx->GetLayerGeom(layer));
	  local = layergeom->get_local_from_world_coords(staveid, 0, 0, chipid, world);
	  local_2D[0] = local[0];
	  local_2D[1] = local[2];

	  if(Verbosity() > 10)
	    {
	      double segcent[3];
	      layergeom->find_sensor_center(staveid, 0, 0, chipid, segcent);
	      cout << "   segment center: " << segcent[0] << " " << segcent[1] << " " << segcent[2] << endl;
	      cout << "   world; " << world[0] << " " << world[1] << " " << world[2] << endl;
	      cout << "   local; " << local[0] << " " << local[1] << " " << local[2] << endl;
	    }

	  // transform covariance matrix back to local coords on chip
	  local_err = GetMvtxCovarLocal(layer, staveid, chipid, world_err);

	}
      else if (trkrid == TrkrDefs::inttId)  	  // INTT
	{
	  unsigned int ladderzid = InttDefs::getLadderZId(cluskey);
	  unsigned int ladderphiid = InttDefs::getLadderPhiId(cluskey);
	  if(Verbosity() > 0) 
	    cout << "    Intt cluster with ladderzid " << ladderzid << " ladderphid " << ladderphiid << endl; 

	  // make identifier for this sensor
	  hsetkey = InttDefs::genHitSetKey(layer, ladderzid, ladderphiid);
	  // get the TGeoNode for it
	  std::map<TrkrDefs::hitsetkey, TGeoNode*>::iterator it;
	  it = _cluster_node_map.find(hsetkey);
	  if(it == _cluster_node_map.end())
	    {
	      cout << PHWHERE << " Did not find entry in TGeo map for this cluster. That should be impossible!" << endl;
	      return Fun4AllReturnCodes::ABORTEVENT; 	      
	    }

	  sensor_node = it->second;
	  if(Verbosity() > 0)
	    cout << "      Found in _cluster_node_map:  layer " << layer << " ladderzid " << ladderzid << " ladderphiid " << ladderphiid 
		 <<  " node " << sensor_node->GetName() << endl;;

	  // Find Acts surface corresponding to this cluster
	  
	  std::map<TrkrDefs::hitsetkey, std::shared_ptr<const Acts::Surface>>::iterator surf_iter;
	  surf_iter = _cluster_surface_map.find(hsetkey);  
	  if(surf_iter == _cluster_surface_map.end())
	    {
	      cout << "Failed to find associated surface element - should be impossible " << endl;
	      return Fun4AllReturnCodes::ABORTEVENT;
	    }

	  //TrkrDefs::hitsetkey found_key = surf_iter->first;
	  surf = surf_iter->second->getSharedPtr();
	  // cout<< "Got surface pair " << surf->name() << " surface type " << surf->type() << std::endl;
	  // surf->toStream(geo_ctxt, std::cout);
	  
	  // transform position back to local coords on sensor
	  CylinderGeomIntt *layergeom = dynamic_cast<CylinderGeomIntt *>(_geom_container_intt->GetLayerGeom(layer));
	  local = layergeom->get_local_from_world_coords(ladderzid, ladderphiid, world);
	  local_2D[0] = local[1];     // r*phi
	  local_2D[1] = local[2];    // z

	  if(Verbosity() > 10)
	    {
	      double segcent[3];
	      layergeom->find_segment_center(ladderzid, ladderphiid,segcent);
	      cout << "   segment center: " << segcent[0] << " " << segcent[1] << " " << segcent[2] << endl;
	      cout << "   world; " << world[0] << " " << world[1] << " " << world[2] << endl;
	      cout << "   local; " << local[0] << " " << local[1] << " " << local[2] << endl;
	    }

	  local_err = GetInttCovarLocal(layer, ladderzid, ladderphiid, world_err);

	}
      else        // TPC
	{
	  double clusphi = atan2(world[1], world[0]);
	  double r_clusphi = radius*clusphi;
	  double ztpc = world[2];
	  
	  // figure out which layer, module, side this cluster is in
	  unsigned int tpc_layer = layer;
	  unsigned int sectorid = TpcDefs::getSectorId(cluskey);
	  unsigned int side = TpcDefs::getSide(cluskey);

	  // use the cluster coords to assign the subsurface indices
	  double module_phi_low = ModulePhiStart + (double) sectorid * ModuleStepPhi;
	  unsigned int iphi = (clusphi - module_phi_low) / SurfStepPhi;
	  unsigned int iz = fabs(ztpc) / SurfStepZ;
	  unsigned int i_phi_z = iphi + 100*iz;  // for making map key

	  if(Verbosity() > 0)
	    {
	      double check_surf_rphi_center = radius * (module_phi_low + (double) iphi * SurfStepPhi + SurfStepPhi / 2.0);
	      double check_surf_z_center = (double) iz * SurfStepZ + SurfStepZ / 2.0;
	      if(side == 0)  check_surf_z_center = - check_surf_z_center;
	      
	      std::cout << " xworld " << world[0] << " yworld " << world[1] << " clusphi " << clusphi << std::endl;
	      std::cout << " sectorid " << sectorid << " side " << side << " iphi " << iphi << " iz " << iz << " i_phi_z " << i_phi_z << std::endl;
	      cout <<     "    r_clusphi " << r_clusphi << " check_surf_rphi center " << check_surf_rphi_center 
		   << " ztpc " <<  ztpc  << " check_surf_z_center " << check_surf_z_center 
		   << std::endl;
	    }

	  // get surface
	  TrkrDefs::cluskey surfkey = TpcDefs::genClusKey(tpc_layer, sectorid, side, i_phi_z);
	  std::map<TrkrDefs::cluskey, std::shared_ptr<const Acts::Surface>>::iterator surf_iter;
	  surf_iter = _cluster_surface_map_tpc.find(surfkey);  
	  if(surf_iter == _cluster_surface_map_tpc.end())
	    {
	      std::cout << PHWHERE << "Failed to find surface, should be impossible!" << std::endl;
	      return Fun4AllReturnCodes::ABORTEVENT;
	    }

	  //TrkrDefs::cluskey found_key = surf_iter->first;
	  surf = surf_iter->second->getSharedPtr();

	  // transformation of cluster to local surface coords
	  // Coords are r*phi relative to surface r-phi center, and z relative to surface z center
	  Acts::Vector3D center = surf->center(geo_ctxt);
	  double surf_rphi_center = atan2(center[1], center[0]) * radius;
	  double surf_z_center = center[2];
	  if(Verbosity() > 0)
	    {
	      std::cout << "    surface center readback:   x " << center[0] << " y " << center[1] << " phi " << atan2(center[1], center[0]) << " z " << center[2] 
			<< " surf_rphi_center " << surf_rphi_center << " surf_z_center " << surf_z_center 
			<< std::endl;
	    }

	  local_2D[0] = r_clusphi - surf_rphi_center;
	  local_2D[1] = ztpc - surf_z_center;

	  if(Verbosity() > 0)
	    {
	      cout <<     "    r_clusphi " << r_clusphi << " surf_rphi center " << surf_rphi_center 
		   << " ztpc " <<  ztpc  << " surf_z_center " << surf_z_center 
		   << " local rphi " << local_2D[0] << " local z " << local_2D[1]
		   << std::endl;
	    }

	  local_err = TransformCovarToLocal(clusphi, world_err);

	}
      
      //====================================================
      // Finished with detector specific cluster stuff 
      // We have the data needed to construct an Acts  measurement for this cluster
      //====================================================

      // is local_err the correct covariance?
      Acts::ActsSymMatrixD<3> cov;   	  // ActsSymMatrix = Eigen::Matrix<double, rows, rows> 
      cov << local_err[0][0],  local_err[0][1],  local_err[0][2],
	local_err[1][0], local_err[1][1], local_err[1][2],
	local_err[2][0], local_err[2][1], local_err[2][2];
      
      // local and local_err contain the position and covariance matrix in local coords
      if(Verbosity() > 0)
	{
	  std::cout << "    layer " << layer << std::endl;
	  for(int i=0;i<2;++i)
	    {
	      cout << "    i " << i << "   local 2D position " << local_2D[i]  << endl;
	    }
	  
	  std::cout << "    local covariance matrix:" << std::endl;
	  std::cout << cov << std::endl;
	}
 
      // Cluster positions on GeoObject/Surface
      Acts::BoundVector loc = Acts::BoundVector::Zero();     
      loc[Acts::eLOC_0]     = local_2D[0];  
      loc[Acts::eLOC_1]     = local_2D[1];

      if(Verbosity() > 0)
	{      
	  std::cout << "    Layer " << layer << " create measurement for trkrid " << trkrid 
		    << " surface " << surf->name() << " surface type " << surf->type() 
		    << " local x " << loc[Acts::eLOC_0] << " local y " << loc[Acts::eLOC_1]
		    << endl;
	}

      if(trkrid == TrkrDefs::mvtxId || trkrid == TrkrDefs::inttId || trkrid == TrkrDefs::tpcId)
	{
	  // TrkrClusterSourceLink takes care of creating an Acts::FittableMeasurement
	  TrkrClusterSourceLink sourceLink(cluskey,surf,loc,cov);
	  // Can obtain the measurement like this
	  auto measurement = *sourceLink;
	  // Have to store this in a container 
	  // DigitizationAlgorithm.cpp makes a map
	  // what does tracking need?
	  // KF tracking needs SourceLinks, track seeds (I think), and options for how
	  // to run KF tracking
	  // Does TGeoNode come into it?
	  
	}
    }
  
  // _trackmap is SvtxTrackMap from the node tree
  // We need to convert to Acts tracks
  for (SvtxTrackMap::Iter iter = _trackmap->begin(); iter != _trackmap->end();
       ++iter)
    {
      SvtxTrack* svtx_track = iter->second;
      if(Verbosity() > 0)
	{
	  cout << "   found SVTXTrack " << iter->first << endl;
	  svtx_track->identify();
	}
      if (!svtx_track)
	continue;
      
    }
  return 0;
}

// methods for converting TrkrCluster data to what Acts needs

TMatrixD PHActsTrkFitter::GetMvtxCovarLocal(const unsigned int layer, const unsigned int staveid, const unsigned int chipid, TMatrixD world_err)
{
  TMatrixD local_err(3,3);

  // rotate errors back to local coords
  double ladder_location[3] = {0.0, 0.0, 0.0};

  // returns the center of the sensor in world coordinates - used to get the ladder phi location
  CylinderGeom_Mvtx *layergeom = dynamic_cast<CylinderGeom_Mvtx *>(_geom_container_mvtx->GetLayerGeom(layer));
  layergeom->find_sensor_center(staveid, 0, 0, chipid, ladder_location);
  double ladderphi = atan2(ladder_location[1], ladder_location[0]);
  ladderphi += layergeom->get_stave_phi_tilt();

  local_err = TransformCovarToLocal(ladderphi, world_err);
  
  if(Verbosity() > 10)
    {
      for(int i=0;i<3;++i)
	for(int j = 0; j<3; ++j)
	  {
	    cout << "  " << i << "    " << j << " local_err " << local_err[i][j] << endl;
	  }
    }
    
  return local_err;
}

TMatrixD PHActsTrkFitter::GetInttCovarLocal(const unsigned int layer, const unsigned int ladderzid, const unsigned int ladderphiid, TMatrixD world_err)
{
  TMatrixD local_err(3,3);

  // rotate errors back to local coords
  double ladder_location[3] = {0.0, 0.0, 0.0};

  // rotate errors back to local coords 
  CylinderGeomIntt *layergeom = dynamic_cast<CylinderGeomIntt *>(_geom_container_intt->GetLayerGeom(layer));
  layergeom->find_segment_center(ladderzid, ladderphiid, ladder_location);
  double ladderphi = atan2(ladder_location[1], ladder_location[0]);

  local_err = TransformCovarToLocal(ladderphi, world_err);
  
  if(Verbosity() > 10)
    {
      for(int i=0;i<3;++i)
	for(int j = 0; j<3; ++j)
	  {
	    cout << "  INTT: " << i << "    " << j << " local_err " << local_err[i][j] << endl;
	  }
    }
    
  return local_err;
}

TMatrixD PHActsTrkFitter::TransformCovarToLocal(const double ladderphi, TMatrixD world_err)
{
  TMatrixD local_err(3,3);
  
  // this is the matrix that was used to rotate from local to global coords 
  TMatrixD ROT(3, 3);
  ROT[0][0] = cos(ladderphi);
  ROT[0][1] = -1.0 * sin(ladderphi);
  ROT[0][2] = 0.0;
  ROT[1][0] = sin(ladderphi);
  ROT[1][1] = cos(ladderphi);
  ROT[1][2] = 0.0;
  ROT[2][0] = 0.0; 
  ROT[2][1] = 0.0;
  ROT[2][2] = 1.0;
  // we want the inverse rotation
  ROT.Invert();
  
  TMatrixD ROT_T(3, 3);
  ROT_T.Transpose(ROT);
  
  local_err = ROT * world_err * ROT_T;

  return local_err;

}
 
int PHActsTrkFitter::End(PHCompositeNode* topNode)
{

  return Fun4AllReturnCodes::EVENT_OK;
}

PHActsTrkFitter::~PHActsTrkFitter()
{

}


int PHActsTrkFitter::CreateNodes(PHCompositeNode* topNode)
{
  
  return Fun4AllReturnCodes::EVENT_OK;
}

/*
 * GetNodes():
 *  Get all the all the required nodes off the node tree
 */

int PHActsTrkFitter::GetNodes(PHCompositeNode* topNode)
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
  

  // Input Trkr Clusters
  _clustermap = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");
  if (!_clustermap)
  {
    cout << PHWHERE << " TRKR_CLUSTER node not found on node tree"
         << endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  // Input Svtx Tracks
  _trackmap = findNode::getClass<SvtxTrackMap>(topNode, "SvtxTrackMap");
  if (!_trackmap)
  {
    cout << PHWHERE << " SvtxTrackMap node not found on node tree"
         << endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }


  /*
  // Input Svtx Vertices
  _vertexmap = findNode::getClass<SvtxVertexMap>(topNode, "SvtxVertexMap");
  if (!_vertexmap && _event < 2)
  {
    cout << PHWHERE << " SvtxVertexrMap node not found on node tree"
         << endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }
  */

  /*(
  // Output Svtx Tracks
  if (!(_over_write_svtxtrackmap) || _output_mode == DebugMode)
  {
    _trackmap_refit = findNode::getClass<SvtxTrackMap>(topNode,
                                                       "SvtxTrackMapRefit");
    if (!_trackmap_refit && _event < 2)
    {
      cout << PHWHERE << " SvtxTrackMapRefit node not found on node tree"
           << endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }
  }
  */

  return Fun4AllReturnCodes::EVENT_OK;
}



