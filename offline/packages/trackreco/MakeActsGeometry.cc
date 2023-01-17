/*!
 *  \file		MakeActsGeometry.cc
 *  \brief		Refit SvtxTracks with PHActs.
 *  \details	Refit SvtxTracks with PHActs.
 *  \author	        Tony Frawley <afrawley@fsu.edu>
 */

#include "MakeActsGeometry.h"

#include <trackbase/TrkrDefs.h>
#include <trackbase/InttDefs.h>
#include <trackbase/MvtxDefs.h>
#include <trackbase/TpcDefs.h>
#include <trackbase/sPHENIXActsDetectorElement.h>
#include <trackbase/AlignmentTransformation.h>
#include <trackbase/alignmentTransformationContainer.h>

#include <intt/CylinderGeomIntt.h>

#include <mvtx/CylinderGeom_Mvtx.h>

#include <micromegas/CylinderGeomMicromegas.h>

#include <micromegas/MicromegasDefs.h>

#include <g4detectors/PHG4TpcCylinderGeom.h>
#include <g4detectors/PHG4TpcCylinderGeomContainer.h>
#include <g4detectors/PHG4CylinderGeom.h>  // for PHG4CylinderGeom
#include <g4detectors/PHG4CylinderGeomContainer.h>

#include <phgeom/PHGeomUtility.h>
#include <phgeom/PHGeomIOTGeo.h>
#include <phgeom/PHGeomTGeo.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/Fun4AllServer.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHDataNode.h>
#include <phool/PHNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>
#include <phool/getClass.h>
#include <phool/phool.h>

#include <Acts/Geometry/GeometryContext.hpp>
#include <Acts/Geometry/TrackingVolume.hpp>
#include <Acts/MagneticField/MagneticFieldContext.hpp>
#include <Acts/Surfaces/PerigeeSurface.hpp>
#include <Acts/Surfaces/PlaneSurface.hpp>
#include <Acts/Surfaces/Surface.hpp>
#include <Acts/Utilities/CalibrationContext.hpp>

#include <ActsExamples/Framework/IContextDecorator.hpp>
#include <ActsExamples/Geometry/CommonGeometry.hpp>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#include <ActsExamples/Options/CommonOptions.hpp>
#pragma GCC diagnostic pop

#include <ActsExamples/Utilities/Options.hpp>
#include <ActsExamples/Options/MagneticFieldOptions.hpp>

#include <ActsExamples/TGeoDetector/JsonTGeoDetectorConfig.hpp>

#include <ActsExamples/Geometry/MaterialWiper.hpp>
#include <Acts/Material/IMaterialDecorator.hpp>
#include <Acts/Plugins/Json/JsonMaterialDecorator.hpp>
#include <Acts/Plugins/Json/MaterialMapJsonConverter.hpp>

#include <TGeoManager.h>
#include <TMatrixT.h>
#include <TObject.h>
#include <TSystem.h>
#include <TVector3.h>

#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <map>
#include <memory>
#include <utility>
#include <vector>
#include <filesystem>

namespace
{
  /// navigate Acts volumes to find one matching a given name (recursive)
  TrackingVolumePtr find_volume_by_name( const Acts::TrackingVolume* master, const std::string& name )
  {
    // skip if name is not composite
    if( master->volumeName().empty() || master->volumeName()[0] != '{' ) return nullptr;

    // loop over children
    for( const auto& child:master->confinedVolumes()->arrayObjects() )
    {
      if( child->volumeName() == name ) return child;
      else if( auto found = find_volume_by_name( child.get(), name ) ) return found;
    }

    // not found
    return nullptr;
  }
}

MakeActsGeometry::MakeActsGeometry(const std::string &name)
: SubsysReco(name)
{
  for ( const auto id : { TrkrDefs::mvtxId, TrkrDefs::inttId, TrkrDefs::tpcId, TrkrDefs::micromegasId })
    {
      m_misalignmentFactor.insert(std::make_pair(id, 1.));
    }
}

int MakeActsGeometry::Init(PHCompositeNode */*topNode*/)
{  

  return Fun4AllReturnCodes::EVENT_OK;
}

int MakeActsGeometry::InitRun(PHCompositeNode *topNode)
{
    m_geomContainerTpc =
      findNode::getClass<PHG4TpcCylinderGeomContainer>(topNode, "CYLINDERCELLGEOM_SVTX");

  // Alignment Transformation declaration of instance - must be here to set initial alignment flag
  AlignmentTransformation alignment_transformation;
  alignment_transformation.createAlignmentTransformContainer(topNode);

  //set parameter for sampling probability distribution
  if(mvtxParam){alignment_transformation.setMVTXParams(m_mvtxDevs);}
  if(inttParam){alignment_transformation.setINTTParams(m_inttDevs);}
  if(tpcParam){alignment_transformation.setTPCParams(m_tpcDevs);}
  if(mmParam){alignment_transformation.setMMParams(m_mmDevs);}

  if(buildAllGeometry(topNode) != Fun4AllReturnCodes::EVENT_OK)
    return Fun4AllReturnCodes::ABORTEVENT;

  /// Set the actsGeometry struct to be put on the node tree
  ActsTrackingGeometry trackingGeometry;
  trackingGeometry.tGeometry = m_tGeometry;
  trackingGeometry.magField = m_magneticField;
  trackingGeometry.getGeoContext() = m_geoCtxt;  // set reference as plain geocontext
  trackingGeometry.tpcSurfStepPhi = m_surfStepPhi;
  trackingGeometry.tpcSurfStepZ = m_surfStepZ;

  // fill ActsSurfaceMap content
  ActsSurfaceMaps surfMaps;
  surfMaps.m_siliconSurfaceMap = m_clusterSurfaceMapSilicon;
  surfMaps.m_tpcSurfaceMap = m_clusterSurfaceMapTpcEdit;
  surfMaps.m_mmSurfaceMap = m_clusterSurfaceMapMmEdit ;
  surfMaps.m_tGeoNodeMap = m_clusterNodeMap;

  // fill TPC volume ids
  for( const auto& [hitsetid, surfaceVector]:m_clusterSurfaceMapTpcEdit )
    for( const auto& surface:surfaceVector )
      { surfMaps.m_tpcVolumeIds.insert( surface->geometryId().volume() ); }
  
  // fill Micromegas volume ids
  for( const auto& [hitsetid, surface]:m_clusterSurfaceMapMmEdit )
    { surfMaps.m_micromegasVolumeIds.insert( surface->geometryId().volume() ); } 
 
  m_actsGeometry->setGeometry(trackingGeometry);
  m_actsGeometry->setSurfMaps(surfMaps);
  m_actsGeometry->set_drift_velocity(m_drift_velocity);

  alignment_transformation.createMap(topNode);
  for(auto& [id, factor] : m_misalignmentFactor)
    {
      alignment_transformation.misalignmentFactor(id, factor);
    }

 // print
  if( Verbosity() )
    {
      for( const auto& id:surfMaps.m_tpcVolumeIds )
	{ std::cout << "MakeActsGeometry::InitRun - TPC volume id: " << id << std::endl; }
      
      for( const auto& id:surfMaps.m_micromegasVolumeIds )
	{ std::cout << "MakeActsGeometry::InitRun - Micromegas volume id: " << id << std::endl; }
    }
  
  return Fun4AllReturnCodes::EVENT_OK;
}

int MakeActsGeometry::buildAllGeometry(PHCompositeNode *topNode)
{

  /// Add the TPC surfaces to the copy of the TGeoManager. 
  // this also adds the micromegas surfaces
  // Do this before anything else, so that the geometry is finalized
  
  if(getNodes(topNode) != Fun4AllReturnCodes::EVENT_OK)
    return Fun4AllReturnCodes::ABORTEVENT;  

  setPlanarSurfaceDivisions();//eshulga
  // This should be done only on the first tracking pass, to avoid adding surfaces twice. 
  // There is a check for existing acts fake surfaces in editTPCGeometry
  editTPCGeometry(topNode); 
 
  /// Export the new geometry to a root file for examination
  if(Verbosity() > 3)
    {
      PHGeomUtility::ExportGeomtry(topNode, "sPHENIXActsGeom.root"); 
      PHGeomUtility::ExportGeomtry(topNode, "sPHENIXActsGeom.gdml");
    }

  // need to get nodes first, in order to be able to build the proper micromegas geometry
  //if(getNodes(topNode) != Fun4AllReturnCodes::EVENT_OK)
  //  return Fun4AllReturnCodes::ABORTEVENT;  

  // In case we did not call EditTpcGeometry, we still want to make the MMs surface map
  if( m_buildMMs && !m_geomContainerMicromegas )
  {
    std::cout << "MakeActsGeometry::buildAllGeometry - could not find node CYLINDERGEOM_MICROMEGAS_FULL. Disabling building MM surfaces" << std::endl;
    m_buildMMs = false;
  }

  if(createNodes(topNode) != Fun4AllReturnCodes::EVENT_OK)
    return Fun4AllReturnCodes::ABORTEVENT;

  /// Run Acts layer builder
  buildActsSurfaces();

  /// Create a map of sensor TGeoNode pointers using the TrkrDefs:: hitsetkey as the key
  //makeTGeoNodeMap(topNode);



  return Fun4AllReturnCodes::EVENT_OK;
}

void MakeActsGeometry::editTPCGeometry(PHCompositeNode *topNode)
{
  // Because we reset and rebuild the geomNode, we do edits of the TPC geometry in the same module

  PHGeomTGeo *geomNode = PHGeomUtility::GetGeomTGeoNode(topNode, true);
  assert(geomNode);

  /// Reset the geometry node, which we will recreate with the TPC edits
  if(geomNode->isValid())
    {
      geomNode->Reset();
    }

  PHGeomIOTGeo *dstGeomIO = PHGeomUtility::GetGeomIOTGeoNode(topNode, false);
  assert(dstGeomIO);
  assert(dstGeomIO->isValid());

  TGeoManager *geoManager = dstGeomIO->ConstructTGeoManager();
  geomNode->SetGeometry(geoManager);
  assert(geoManager);
 
  if(TGeoManager::GetDefaultUnits() != TGeoManager::EDefaultUnits::kRootUnits)
    {
      std::cerr << "There is a potential unit mismatch in the ROOT geometry."
		<< " It is dangerous to continue and may lead to inconsistencies in the Acts geometry. Exiting." 
		<< std::endl;
      gSystem->Exit(1);
    }

  TGeoVolume *World_vol = geoManager->GetTopVolume();
  
  // TPC geometry edits
  //===============

  TGeoNode *tpc_envelope_node = nullptr;
  TGeoNode *tpc_gas_north_node = nullptr;

  // find tpc north gas volume at path of World*/tpc_envelope*
  if (Verbosity()> 3)
  {
    std::cout << "EditTPCGeometry - searching under volume: ";
    World_vol->Print();
  }
  for (int i = 0; i < World_vol->GetNdaughters(); i++)
  {
    TString node_name = World_vol->GetNode(i)->GetName();
    if (node_name.BeginsWith("tpc_envelope"))
    {
      if (Verbosity())
        std::cout << "EditTPCGeometry - found " << node_name << std::endl;

      tpc_envelope_node = World_vol->GetNode(i);
      break;
    }
  }

  assert(tpc_envelope_node);

  // find tpc north gas volume at path of World*/tpc_envelope*/tpc_gas_north*
  TGeoVolume *tpc_envelope_vol = tpc_envelope_node->GetVolume();
  assert(tpc_envelope_vol);
  if (Verbosity() > 3)
    {
      std::cout << "EditTPCGeometry - searching under volume: ";
      tpc_envelope_vol->Print();
    }
  
  for (int i = 0; i < tpc_envelope_vol->GetNdaughters(); i++)
    {
      TString node_name = tpc_envelope_vol->GetNode(i)->GetName();
      
      if (node_name.BeginsWith("tpc_gas_north"))
	{
	  if (Verbosity())
	    std::cout << "EditTPCGeometry - found " << node_name << std::endl;
	  
	  tpc_gas_north_node = tpc_envelope_vol->GetNode(i);
	  break;
	}
    }
  
  assert(tpc_gas_north_node);
  TGeoVolume *tpc_gas_north_vol = tpc_gas_north_node->GetVolume();
  assert(tpc_gas_north_vol);

  int nfakesurfaces = 0;
  for(int i=0; i<tpc_gas_north_vol->GetNdaughters(); i++)
    {
      TString node_name = tpc_gas_north_vol->GetNode(i)->GetName();
      if(node_name.BeginsWith("tpc_gas_measurement_"))
	{ nfakesurfaces++; }
    }

  /// Make a check for the fake surfaces. If we have more than 0
  /// then we've built the fake surfaces and we should not do it again
  if(nfakesurfaces > 0)
    { return; }

  if (Verbosity() > 3)
  {
    std::cout << "EditTPCGeometry - gas volume: ";
    tpc_gas_north_vol->Print();
  }

  // adds surfaces to the underlying volume, so both north and south placements get them
  addActsTpcSurfaces(tpc_gas_north_vol, geoManager);

  // done
  geoManager->CloseGeometry();

  // save the edited geometry to DST persistent IO node for downstream DST files
  PHGeomUtility::UpdateIONode(topNode);
 
}

void MakeActsGeometry::addActsTpcSurfaces(TGeoVolume *tpc_gas_vol,
					  TGeoManager *geoManager)
{
  TGeoMedium *tpc_gas_medium = tpc_gas_vol->GetMedium();
  assert(tpc_gas_medium);

  // printout
  std::cout << "MakeActsGeometry::addActsTpcSurfaces - m_nSurfPhi: " << m_nSurfPhi << std::endl;
  
  TGeoVolume *tpc_gas_measurement_vol[m_nTpcLayers];
  double tan_half_phi = tan(m_surfStepPhi / 2.0);
  int copy = 0;
  for(unsigned int ilayer = 0; ilayer < m_nTpcLayers; ++ilayer)
    {
      // make a box for this layer
      char bname[500];
      sprintf(bname,"tpc_gas_measurement_%u",ilayer);

      // Because we use a box, not a section of a cylinder, we need this to prevent overlaps
      // set the nominal r*phi dimension of the box so they just touch at the inner edge when placed 
      double box_r_phi = 2.0 * tan_half_phi * 
	(m_layerRadius[ilayer] - m_layerThickness[ilayer] / 2.0);

  
      tpc_gas_measurement_vol[ilayer] = geoManager->MakeBox(bname, tpc_gas_medium, 
							    m_layerThickness[ilayer]*half_width_clearance_thick, 
							    box_r_phi*half_width_clearance_phi, 
							    m_surfStepZ*half_width_clearance_z);

      tpc_gas_measurement_vol[ilayer]->SetLineColor(kBlack);
      tpc_gas_measurement_vol[ilayer]->SetFillColor(kYellow);
      tpc_gas_measurement_vol[ilayer]->SetVisibility(kTRUE);

      if(Verbosity() > 3)
	{
	  std::cout << " Made box for layer " << ilayer 
		    << " with dx " << m_layerThickness[ilayer] << " dy " 
		    << box_r_phi << " ref arc " 
		    << m_surfStepPhi * m_layerRadius[ilayer] << " dz " 
		    << m_surfStepZ << std::endl;
	  tpc_gas_measurement_vol[ilayer]->Print();
	  tpc_gas_measurement_vol[ilayer]->CheckOverlaps();
	}
      	      
      for (unsigned int iz = 0; iz < m_nSurfZ; ++iz)
	{
	  // The (half) tpc gas volume is 105.5 cm long and is symmetric around (x,y,z) = (0,0,0) in its frame
	  double z_center = 0.0;
	  
	  for (unsigned int imod = 0; imod < m_nTpcModulesPerLayer; ++imod)
	    {
	      for (unsigned int iphi = 0; iphi < m_nSurfPhi; ++iphi)
		{
		  
		  double min_phi = m_modulePhiStart + 
		    (double) imod * m_moduleStepPhi + 
		    (double) iphi * m_surfStepPhi;
		  double phi_center = min_phi + m_surfStepPhi / 2.0;
		  double phi_center_degrees = phi_center * 180.0 / M_PI;
        
		  // place copies of the gas volume to fill up the layer
		  
		  double x_center = m_layerRadius[ilayer] * cos(phi_center);
		  double y_center = m_layerRadius[ilayer] * sin(phi_center);
		  
		  char rot_name[500];
		  sprintf(rot_name,"tpc_gas_rotation_%i", copy);
		  TGeoCombiTrans *tpc_gas_measurement_location 
		    = new TGeoCombiTrans(x_center, y_center, z_center,
					 new TGeoRotation(rot_name,
							  phi_center_degrees, 
							  0, 0));
		 
		  tpc_gas_vol->AddNode(tpc_gas_measurement_vol[ilayer], 
				       copy, tpc_gas_measurement_location);
        
		  copy++;
		  if(Verbosity() > 3 ) 
		    {
        
		      std::cout << "Box center : ("<<x_center<<", " <<y_center
				<< ", " << z_center << ")" << " and in rphiz "
				<< sqrt(x_center*x_center+y_center*y_center) 
				<< ", " << atan2(y_center,x_center) << ", " 
				<< z_center << std::endl;
		      std::cout << "Box dimensions " <<m_layerThickness[ilayer] *half_width_clearance_thick
				<<" , " << box_r_phi/(m_layerRadius[ilayer]-m_layerThickness[ilayer]/2.) * half_width_clearance_phi 
				<< ", " << m_surfStepZ*half_width_clearance_z << " and in xyz " 
				<< m_layerThickness[ilayer]*half_width_clearance_thick*cos(box_r_phi/(m_layerRadius[ilayer]-m_layerThickness[ilayer]/2.)*half_width_clearance_phi) << ", " 
				<< m_layerThickness[ilayer]*half_width_clearance_thick*sin(box_r_phi/(m_layerRadius[ilayer]-m_layerThickness[ilayer]/2.)*half_width_clearance_phi) << ", " 
				<<m_surfStepZ*half_width_clearance_z<<std::endl;
		    }
		}
	    }
	}
    }
	}

/**
 * Builds silicon layers and TPC geometry in the ACTS surface world
 */
void MakeActsGeometry::buildActsSurfaces()
{

  // define int argc and char* argv to provide options to processGeometry
  const int argc = 20;
  char* arg[argc];
  m_magFieldRescale = 1;
  // if(Verbosity() > 0)
    std::cout << PHWHERE << "Magnetic field " << m_magField 
	      << " with rescale " << m_magFieldRescale << std::endl;

  std::string responseFile, materialFile;
  setMaterialResponseFile(responseFile, materialFile);

  // Response file contains arguments necessary for geometry building
  std::string argstr[argc]{
    "-n1",
    "--geo-tgeo-jsonconfig", responseFile,
      "--mat-input-type","file",
      "--mat-input-file", materialFile,
      "--bf-constant-tesla","0:0:1.4",
      "--bf-bscalor"};
  
  argstr[9] = std::to_string(m_magFieldRescale);
     

  /// Alter args if using field map
  if(m_magField.find(".root") != std::string::npos)
    {
      if(m_magField.find("2d") != std::string::npos)
	{        
	  m_magFieldRescale = 1;
	}
      char *calibrationsroot = getenv("CALIBRATIONROOT");
      m_magField = "sphenix3dtrackingmapxyz.root";
      if (calibrationsroot != nullptr)
      {
	m_magField = std::string(calibrationsroot) + std::string("/Field/Map/") + m_magField;
      }
      //m_magField = std::string("/phenix/u/bogui/data/Field/sphenix3dtrackingmapxyz.root");
      //m_magField = std::string("/phenix/u/bogui/data/Field/sphenix3dbigmapxyz.root");
      
      argstr[7] = "--bf-map-file";
      argstr[8] = m_magField;
      argstr[9]= "--bf-map-tree";
      argstr[10] = "fieldmap";
      argstr[11] = "--bf-map-lengthscale-mm";
      argstr[12] = "10";
      argstr[13] = "--bf-map-fieldscale-tesla";
      argstr[14] = std::to_string(m_magFieldRescale);  
      
    }

  //  if(Verbosity() > 0)
    std::cout << "Mag field now " << m_magField << " with rescale "
	      << m_magFieldRescale << std::endl;

  // Set vector of chars to arguments needed
  for (int i = 0; i < argc; ++i)
    {
      if(Verbosity() > 1)
	std::cout << argstr[i] << ", ";
      // need a copy, since .c_str() returns a const char * and process geometry will not take a const
      arg[i] = strdup(argstr[i].c_str());
    }

  // We replicate the relevant functionality of  
  //acts/Examples/Run/Common/src/GeometryExampleBase::ProcessGeometry() in MakeActsGeometry()
  // so we get access to the results. The layer builder magically gets the TGeoManager
  
  makeGeometry(argc, arg, m_detector);

  for(int i=0; i<argc; i++)
    free(arg[i]);

}
void MakeActsGeometry::setMaterialResponseFile(std::string& responseFile,
					       std::string& materialFile)
{
 
  responseFile = m_buildMMs ? "tgeo-sphenix-mms.json":"tgeo-sphenix.json";
  materialFile = m_buildMMs ? "sphenix-mm-material.json":"sphenix-material.json";
  /// Check to see if files exist locally - if not, use defaults
  std::ifstream file;

  file.open(responseFile);
  if(!file.is_open())
    {
      std::cout << responseFile
		<< " not found locally, use repo version"
		<< std::endl;
      char *offline_main = getenv("OFFLINE_MAIN");
      assert(offline_main);
      responseFile = std::string(offline_main) +
	(m_buildMMs ? "/share/tgeo-sphenix-mms.json":"/share/tgeo-sphenix.json");
    }

  file.open(materialFile);
  if(!file.is_open())
    {
      std::cout << materialFile 
		<< " not found locally, use repo version" 
		<< std::endl;
      materialFile = std::string(getenv("CALIBRATIONROOT")) +
	(m_buildMMs ? "/ACTS/sphenix-mm-material.json":"/ACTS/sphenix-material.json");
    }
  
  if(Verbosity() > -1)
    {
      std::cout << "using Acts material file : " << materialFile 
		<< std::endl;
      std::cout << "Using Acts TGeoResponse file : " << responseFile
		<< std::endl;
    }
  
  return;

}
void MakeActsGeometry::makeGeometry(int argc, char* argv[], 
				    ActsExamples::TGeoDetector &detector)
{
  
  /// setup and parse options
  boost::program_options::options_description desc;
  ActsExamples::Options::addGeometryOptions(desc);
  ActsExamples::Options::addMaterialOptions(desc);
  ActsExamples::Options::addMagneticFieldOptions(desc);

  /// Add specific options for this geometry
  detector.addOptions(desc);
  auto vm = ActsExamples::Options::parse(desc, argc, argv);
 
  /// The geometry, material and decoration
  auto geometry = build(vm, detector);
  /// Geometry is a pair of (tgeoTrackingGeometry, tgeoContextDecorators)

  m_tGeometry = geometry.first;
  if(m_useField)
    { m_magneticField = ActsExamples::Options::readMagneticField(vm); }
  else
    { m_magneticField = nullptr; }

  m_geoCtxt = Acts::GeometryContext();
 
  unpackVolumes();
  
  return;
}


std::pair<std::shared_ptr<const Acts::TrackingGeometry>,
          std::vector<std::shared_ptr<ActsExamples::IContextDecorator>>>
MakeActsGeometry::build(const boost::program_options::variables_map& vm,
			ActsExamples::TGeoDetector& detector) {
  // Material decoration
  std::shared_ptr<const Acts::IMaterialDecorator> matDeco = nullptr;
 
  // Retrieve the filename
  auto fileName = vm["mat-input-file"].template as<std::string>();
  // json or root based decorator
  if (fileName.find(".json") != std::string::npos ||
      fileName.find(".cbor") != std::string::npos) 
    {
      // Set up the converter first
      Acts::MaterialMapJsonConverter::Config jsonGeoConvConfig;
      // Set up the json-based decorator
      matDeco = std::make_shared<const Acts::JsonMaterialDecorator>(
	   jsonGeoConvConfig, fileName, Acts::Logging::INFO);
    } 
  else
    {
      matDeco = std::make_shared<const Acts::MaterialWiper>();
    }

  ActsExamples::TGeoDetector::Config config;

  config.elementFactory = sPHENIXElementFactory;

  config.fileName = vm["geo-tgeo-filename"].as<std::string>();

  const auto path = vm["geo-tgeo-jsonconfig"].template as<std::string>();

  readTGeoLayerBuilderConfigsFile(path, config);

  /// Return the geometry and context decorators
  return detector.finalize(config, matDeco);
 }
  
void MakeActsGeometry::readTGeoLayerBuilderConfigsFile(const std::string& path,
						       ActsExamples::TGeoDetector::Config& config) {
  if (path.empty()) {
    std::cout << "There is no acts geometry response file loaded. Cannot build, exiting"
	      << std::endl;
    exit(1);
  }

  nlohmann::json djson;
  std::ifstream infile(path, std::ifstream::in | std::ifstream::binary);
  infile >> djson;

  config.unitScalor = djson["geo-tgeo-unit-scalor"];

  config.buildBeamPipe = djson["geo-tgeo-build-beampipe"];
  if (config.buildBeamPipe) {
    const auto beamPipeParameters =
        djson["geo-tgeo-beampipe-parameters"].get<std::array<double, 3>>();
    config.beamPipeRadius = beamPipeParameters[0];
    config.beamPipeHalflengthZ = beamPipeParameters[1];
    config.beamPipeLayerThickness = beamPipeParameters[2];
  }

  // Fill nested volume configs
  for (const auto& volume : djson["Volumes"]) {
    auto& vol = config.volumes.emplace_back();
    vol = volume;
  }
}

void MakeActsGeometry::unpackVolumes()
{
  /// m_tGeometry is a TrackingGeometry pointer
  /// vol is a TrackingVolume pointer  
  auto vol = m_tGeometry->highestTrackingVolume();

  if(Verbosity() > 3 )
  { 
    std::cout << "MakeActsGeometry::unpackVolumes - top volume: " << vol->volumeName() << std::endl; 
    std::cout << "Before: Mvtx: m_clusterSurfaceMapSilicon size    " << m_clusterSurfaceMapSilicon.size() << std::endl;
    std::cout << "Before: m_clusterSurfaceMapTpc size    " << m_clusterSurfaceMapTpcEdit.size() << std::endl;
  }

  if(m_buildMMs)
  {
    // micromegas
    auto mmBarrel = find_volume_by_name( vol, "MICROMEGAS::Barrel" );
    assert( mmBarrel );
    makeMmMapPairs(mmBarrel);
  }
  else
    {
      std::cout << "WARNING: You are not building the micromegas in your macro! If you intended to, make sure you set Enable::MICROMEGAS=true; otherwise, your macro will seg fault" << std::endl;
    }
  
  {
    // MVTX
    auto mvtxBarrel = find_volume_by_name( vol, "MVTX::Barrel" );
    assert( mvtxBarrel );
    makeMvtxMapPairs(mvtxBarrel);
    if(Verbosity() > 3)
      { std::cout << "After: Mvtx: m_clusterSurfaceMapSilicon size    " << m_clusterSurfaceMapSilicon.size() << std::endl; }
  }

  {
    // INTT
    if(Verbosity() > 3)
      { std::cout << "Before: INTT: m_clusterSurfaceMapSilicon size    " << m_clusterSurfaceMapSilicon.size() << std::endl; }
    auto inttBarrel = find_volume_by_name( vol, "Silicon::Barrel" );
    assert( inttBarrel );
    makeInttMapPairs(inttBarrel);
    if(Verbosity() > 3)
      { std::cout << "After: INTT: m_clusterSurfaceMapSilicon size    " << m_clusterSurfaceMapSilicon.size() << std::endl; }
  }

  {
    // TPC
    auto tpcBarrel = find_volume_by_name( vol, "TPC::Barrel" );
    assert( tpcBarrel );
    makeTpcMapPairs(tpcBarrel);
    if(Verbosity() > 3)
      { std::cout << "After: m_clusterSurfaceMapTpc size    " << m_clusterSurfaceMapTpcEdit.size() << std::endl; }
  }

  return;
}

void MakeActsGeometry::makeTpcMapPairs(TrackingVolumePtr &tpcVolume)
{
  if(Verbosity() > 10)
  { std::cout << "MakeActsGeometry::makeTpcMapPairs - tpcVolume: " << tpcVolume->volumeName() << std::endl; }
   
  auto tpcLayerArray = tpcVolume->confinedLayers();
  auto tpcLayerVector = tpcLayerArray->arrayObjects();

  /// Need to unfold each layer that Acts builds
  for(unsigned int i = 0; i < tpcLayerVector.size(); i++)
    {
      auto surfaceArray = tpcLayerVector.at(i)->surfaceArray();
      if(surfaceArray == NULL){
	continue;
      }
      /// surfaceVector is a vector of surfaces corresponding to the tpc layer
      /// that acts builds
      auto surfaceVector = surfaceArray->surfaces();
      for( unsigned int j = 0; j < surfaceVector.size(); j++)
	{
	  auto surf = surfaceVector.at(j)->getSharedPtr();
	  auto vec3d = surf->center(m_geoCtxt);

	  /// convert to cm
	  std::vector<double> world_center = {vec3d(0) / 10.0, 
					      vec3d(1) / 10.0,
					      vec3d(2) / 10.0};
	
	  TrkrDefs::hitsetkey hitsetkey = getTpcHitSetKeyFromCoords(world_center);
	  unsigned int layer = TrkrDefs::getLayer(hitsetkey);

	  /// If there is already an entry for this hitsetkey, add the surface
	  /// to its corresponding vector
	  //std::map<TrkrDefs::hitsetkey, std::vector<Surface>>::iterator mapIter;
	  std::map<unsigned int, std::vector<Surface>>::iterator mapIter;
	  //mapIter = m_clusterSurfaceMapTpcEdit.find(hitsetkey);
	  mapIter = m_clusterSurfaceMapTpcEdit.find(layer);
	  
	  if(mapIter != m_clusterSurfaceMapTpcEdit.end())
	    {
	      mapIter->second.push_back(surf);
	    }
	  else
	    {
	      /// Otherwise make a new map entry
	      std::vector<Surface> dumvec;
	      dumvec.push_back(surf);
	      std::pair<unsigned int, std::vector<Surface>> tmp = 
		std::make_pair(layer, dumvec);
	      m_clusterSurfaceMapTpcEdit.insert(tmp);
	    }
	  
	}
    }

}

//____________________________________________________________________________________________
void MakeActsGeometry::makeMmMapPairs(TrackingVolumePtr &mmVolume)
{
  if(Verbosity() > 10)
  { std::cout << "MakeActsGeometry::makeMmMapPairs - mmVolume: " << mmVolume->volumeName() << std::endl; }
  const auto mmLayerArray = mmVolume->confinedLayers();
  const auto mmLayerVector = mmLayerArray->arrayObjects();
  
  /// Need to unfold each layer that Acts builds
  for(unsigned int i = 0; i < mmLayerVector.size(); i++)
  {
    auto surfaceArray = mmLayerVector.at(i)->surfaceArray();
    if(!surfaceArray) continue;
    
    /// surfaceVector is a vector of surfaces corresponding to the micromegas layer
    /// that acts builds
    const auto surfaceVector = surfaceArray->surfaces();

    for( unsigned int j = 0; j < surfaceVector.size(); j++)
    {
      auto surface = surfaceVector.at(j)->getSharedPtr();
      auto vec3d = surface->center(m_geoCtxt);
      
      /// convert to cm
      TVector3 world_center( 
        vec3d(0)/Acts::UnitConstants::cm, 
        vec3d(1)/Acts::UnitConstants::cm,
        vec3d(2)/Acts::UnitConstants::cm
      );
      
      // get relevant layer
      int layer = -1;
      CylinderGeomMicromegas* layergeom = nullptr;
      const auto range = m_geomContainerMicromegas->get_begin_end();
      for( auto iter = range.first; iter != range.second; ++iter )
      {
        auto this_layergeom =  static_cast<CylinderGeomMicromegas*>( iter->second );
        if(this_layergeom->check_radius(world_center))
        { 
          layer = iter->first;
          layergeom = this_layergeom;
          break;
        }
      }
      
      if( !layergeom ) 
      {
        std::cout << "MakeActsGeometry::makeMmMapPairs - could not file CylinderGeomMicromegas matching ACTS surface" << std::endl;
        continue;
      }

      // get matching tile
      int tileid = layergeom->find_tile_cylindrical( world_center );
      if( tileid < 0 ) 
      {
        std::cout << "MakeActsGeometry::makeMmMapPairs - could not file Micromegas tile matching ACTS surface" << std::endl;
        continue;
      } 
            
      // get segmentation type
      const auto segmentation_type = layergeom->get_segmentation_type();
      
      // create hitset key and insert surface in map
      const auto hitsetkey = MicromegasDefs::genHitSetKey(layer, segmentation_type, tileid);
      const auto [iter, inserted] = m_clusterSurfaceMapMmEdit.insert( std::make_pair( hitsetkey, surface ) );
      assert( inserted );
    }
  }
}

void MakeActsGeometry::makeInttMapPairs(TrackingVolumePtr &inttVolume)
{
  
  if(Verbosity() > 10)
  { std::cout << "MakeActsGeometry::makeInttMapPairs - inttVolume: " << inttVolume->volumeName() << std::endl; }

  auto inttLayerArray = inttVolume->confinedLayers();

  auto inttLayerVector = inttLayerArray->arrayObjects();
  // inttLayerVector is a std::vector<LayerPtr>
  for (unsigned int i = 0; i < inttLayerVector.size(); i++)
  {
    // Get the ACTS::SurfaceArray from each INTT LayerPtr being iterated over
    auto surfaceArray = inttLayerVector.at(i)->surfaceArray();
    if (surfaceArray == NULL)
      continue;

    // surfaceVector is an Acts::SurfaceVector, vector of surfaces
    auto surfaceVector = surfaceArray->surfaces();

    for (unsigned int j = 0; j < surfaceVector.size(); j++)
    {
      auto surf = surfaceVector.at(j)->getSharedPtr();
      auto vec3d = surf->center(m_geoCtxt);

      double ref_rad[4] = {7.188, 7.732, 9.680, 10.262};

      std::vector<double> world_center = {vec3d(0) / 10.0, vec3d(1) / 10.0, vec3d(2) / 10.0};  // convert from mm to cm
      
      /// The Acts geometry builder combines layers 4 and 5 together, 
      /// and layers 6 and 7 together. We need to use the radius to figure
      /// out which layer to use to get the layergeom
      double layer_rad = sqrt(pow(world_center[0], 2) + pow(world_center[1], 2));

      unsigned int layer = 0;
      for (unsigned int i = 0; i < 4; ++i)
      {
        if (fabs(layer_rad - ref_rad[i]) < 0.1)
          layer = i + 3;
      }

      TrkrDefs::hitsetkey hitsetkey = getInttHitSetKeyFromCoords(layer, world_center);

      // Add this surface to the map
      std::pair<TrkrDefs::hitsetkey, Surface> tmp = make_pair(hitsetkey, surf);
      m_clusterSurfaceMapSilicon.insert(tmp);

      if (Verbosity() > 10)
      {
        // check it is in there
        unsigned int ladderPhi = InttDefs::getLadderPhiId(hitsetkey);
        unsigned int ladderZ = InttDefs::getLadderZId(hitsetkey);

        std::cout << "Layer radius " << layer_rad << " layer " << layer << " ladderPhi " << ladderPhi << " ladderZ " << ladderZ
                  << " recover surface from m_clusterSurfaceMapSilicon " << std::endl;
        std::cout << " surface type " << surf->type() << std::endl;
        auto assoc_layer = surf->associatedLayer();
        std::cout << std::endl
                  << " Layer type " << assoc_layer->layerType() << std::endl;

        auto assoc_det_element = surf->associatedDetectorElement();
        if (assoc_det_element != nullptr)
        {
          std::cout << " Associated detElement has non-null pointer " 
		    << assoc_det_element << std::endl;
          std::cout << std::endl
                    << " Associated detElement found, thickness = " 
		    << assoc_det_element->thickness() << std::endl;
        }
        else
          std::cout << std::endl
                    << " Associated detElement is nullptr " << std::endl;
      }
    }
  }

}

void MakeActsGeometry::makeMvtxMapPairs(TrackingVolumePtr &mvtxVolume)
{
  
  if(Verbosity() > 10)
  { std::cout << "MakeActsGeometry::makeMvtxMapPairs - mvtxVolume: " << mvtxVolume->volumeName() << std::endl; }

  // Now get the LayerArrays corresponding to each volume
  auto mvtxBarrelLayerArray = mvtxVolume->confinedLayers();  // the barrel is all we care about

  // This is a BinnedArray<LayerPtr>, but we care only about the sensitive surfaces
  auto mvtxLayerVector1 = mvtxBarrelLayerArray->arrayObjects();  // the barrel is all we care about

  // mvtxLayerVector has size 3, but only index 1 has any surfaces since the clayersplits option was removed
  // the layer number has to be deduced from the radius
  for (unsigned int i = 0; i < mvtxLayerVector1.size(); ++i)
  {
    // Get the Acts::SurfaceArray from each MVTX LayerPtr being iterated over
    auto surfaceArray = mvtxLayerVector1.at(i)->surfaceArray();
    if (surfaceArray == NULL)
      continue;

    double ref_rad[3] = {2.556, 3.359, 4.134};

    // surfaceVector is an Acts::SurfaceVector, vector of surfaces
    //std::vector<const Surface*>
    auto surfaceVector = surfaceArray->surfaces();
    for (unsigned int j = 0; j < surfaceVector.size(); j++)
    {
      auto surf = surfaceVector.at(j)->getSharedPtr();
      auto vec3d = surf->center(m_geoCtxt);
      std::vector<double> world_center = {vec3d(0) / 10.0, vec3d(1) / 10.0, vec3d(2) / 10.0};  // convert from mm to cm
      double layer_rad = sqrt(pow(world_center[0], 2) + pow(world_center[1], 2));
      unsigned int layer = 0;
      for (unsigned int i = 0; i < 3; ++i)
      {
        if (fabs(layer_rad - ref_rad[i]) < 0.1)
          layer = i;
      }

      TrkrDefs::hitsetkey hitsetkey = getMvtxHitSetKeyFromCoords(layer, world_center);

      // Add this surface to the map
      std::pair<TrkrDefs::hitsetkey, Surface> tmp = make_pair(hitsetkey, surf);

      m_clusterSurfaceMapSilicon.insert(tmp);

      if (Verbosity() > 10)
      {
        unsigned int stave = MvtxDefs::getStaveId(hitsetkey);
        unsigned int chip = MvtxDefs::getChipId(hitsetkey);

        // check it is in there
        std::cout << "Layer radius " << layer_rad << " Layer " 
		  << layer << " stave " << stave << " chip " << chip
                  << " recover surface from m_clusterSurfaceMapSilicon " 
		  << std::endl;
	
        std::map<TrkrDefs::hitsetkey, Surface>::iterator surf_iter = m_clusterSurfaceMapSilicon.find(hitsetkey);
        std::cout << " surface type " << surf_iter->second->type() 
		  << std::endl;
        auto assoc_layer = surf->associatedLayer();
        std::cout << std::endl << " Layer type " 
		  << assoc_layer->layerType() << std::endl;

        auto assoc_det_element = surf->associatedDetectorElement();
        if (assoc_det_element != nullptr)
        {
          std::cout << " Associated detElement has non-null pointer " 
		    << assoc_det_element << std::endl;
          std::cout << std::endl
                    << " Associated detElement found, thickness = " 
		    << assoc_det_element->thickness() << std::endl;
        }
        else
          std::cout << std::endl
                    << " Associated detElement is nullptr " << std::endl;
      }
    }
  }
}

TrkrDefs::hitsetkey MakeActsGeometry::getTpcHitSetKeyFromCoords(std::vector<double> &world)
{
  // Look up TPC surface index values from world position of surface center
  // layer
  unsigned int layer = 999;
  double layer_rad = sqrt(pow(world[0],2) + pow(world[1],2));
  for(unsigned int ilayer=0;ilayer<m_nTpcLayers;++ilayer)
    {
      double tpc_ref_radius_low = 
	m_layerRadius[ilayer] - m_layerThickness[ilayer] / 2.0;
      double tpc_ref_radius_high = 
	m_layerRadius[ilayer] + m_layerThickness[ilayer] / 2.0;
      if(layer_rad >= tpc_ref_radius_low && layer_rad < tpc_ref_radius_high)
	{
	  layer = ilayer;
	  break;
	}
    }

  if(layer >= m_nTpcLayers) 
    {
      std::cout << PHWHERE 
		<< "Error: undefined layer, do nothing world =  " 
		<< world[0] << "  " << world[1] << "  " << world[2] 
		<< " layer " << layer << std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }
  
  layer += 7;

  // side -  from world z sign 
  unsigned int side;
  if(world[2] < 0)
    side = 0;
  else
    side = 1;

  // readout module 
  unsigned int readout_mod = 999;
  double phi_world = atan2(world[1], world[0]);
  for(unsigned int imod=0; imod<m_nTpcModulesPerLayer; ++imod)
    {
      double min_phi = m_modulePhiStart + (double) imod * m_moduleStepPhi;
      double max_phi = m_modulePhiStart + (double) (imod+1) * m_moduleStepPhi;
      if(phi_world >=min_phi && phi_world < max_phi)
	{
	  readout_mod = imod;
	  break;
	}
    }
  if(readout_mod >= m_nTpcModulesPerLayer)
    {
      std::cout << PHWHERE 
		<< "Error: readout_mod is undefined, do nothing  phi_world = " 
		<< phi_world << std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }

  TrkrDefs::hitsetkey hitset_key = TpcDefs::genHitSetKey(layer, readout_mod, side);
  if(Verbosity() > 3)
    if(layer == 30)
      std::cout << "   world = " << world[0] << "  " << world[1] 
		<< "  " << world[2] << " phi_world " 
		<< phi_world*180/M_PI << " layer " << layer
		<< " readout_mod " << readout_mod << " side " << side 
		<< " hitsetkey " << hitset_key<< std::endl;
  
  return hitset_key;
}

TrkrDefs::hitsetkey MakeActsGeometry::getMvtxHitSetKeyFromCoords(unsigned int layer, std::vector<double> &world)
{
  // Look up the MVTX sensor index values from the world position of the surface center

  CylinderGeom_Mvtx *layergeom = dynamic_cast<CylinderGeom_Mvtx *>(m_geomContainerMvtx->GetLayerGeom(layer));
  if (!layergeom)
  {
    std::cout << PHWHERE << "Did not get layergeom for layer " 
	      << layer << std::endl;
    return 0;
  }

  unsigned int stave = 0;
  unsigned int chip = 0;
  layergeom->get_sensor_indices_from_world_coords(world, stave, chip);

  unsigned int strobe = 0;
  TrkrDefs::hitsetkey mvtx_hitsetkey = MvtxDefs::genHitSetKey(layer, stave, chip, strobe);

  return mvtx_hitsetkey;
}

TrkrDefs::hitsetkey MakeActsGeometry::getInttHitSetKeyFromCoords(unsigned int layer, std::vector<double> &world)
{
  // Look up the INTT sensor index values from the world position of the surface center

  CylinderGeomIntt *layergeom = dynamic_cast<CylinderGeomIntt *>(m_geomContainerIntt->GetLayerGeom(layer));
  if (!layergeom)
  {
    std::cout << PHWHERE << "Did not get layergeom for layer " 
	      << layer << std::endl;
    return 0;
  }

  double location[3] = {world[0], world[1], world[2]};
  int segment_z_bin = 0;
  int segment_phi_bin = 0;
  layergeom->find_indices_from_segment_center(segment_z_bin, 
					      segment_phi_bin, location);

  int crossing = 0;
  TrkrDefs::hitsetkey intt_hitsetkey = InttDefs::genHitSetKey(layer, segment_z_bin, segment_phi_bin, crossing);

  return intt_hitsetkey;
}


// void MakeActsGeometry::makeTGeoNodeMap(PHCompositeNode * /*topNode*/)
// {
//   // This is just a diagnostic method
//   // it lets you list all of the nodes by setting print_sensors = true
// 
//   if (!m_geoManager)
//   {
//     std::cout << PHWHERE << " Did not find TGeoManager, quit! " << std::endl;
//     return;
//   }
//   TGeoVolume *topVol = m_geoManager->GetTopVolume();
//   TObjArray *nodeArray = topVol->GetNodes();
// 
//   TIter iObj(nodeArray);
//   while (TObject *obj = iObj())
//   {
//     TGeoNode *node = dynamic_cast<TGeoNode *>(obj);
//     std::string node_str = node->GetName();
// 
//     std::string mvtx("MVTX_Wrapper");
//     std::string intt("ladder");
//     std::string intt_ext("ladderext");
//     std::string tpc("tpc_envelope");
//     std::string micromegas("MICROMEGAS_55");
// 
//     if (node_str.compare(0, mvtx.length(), mvtx) == 0)  // is it in the MVTX?
//     {
//       if (Verbosity() > 2)
// 	std::cout << " node " << node->GetName() << " is the MVTX wrapper"
// 		  << std::endl;
// 
//       /// The Mvtx has an additional wrapper that needs to be unpacked
//       TObjArray *mvtxArray = node->GetNodes();
//       TIter mvtxObj(mvtxArray);
//       while(TObject *mvtx = mvtxObj())
// 	{
// 	  TGeoNode *mvtxNode = dynamic_cast<TGeoNode *>(mvtx);
// 	  if(Verbosity() > 2)
// 	    std::cout << "mvtx node name is " << mvtxNode->GetName()
// 		      << std::endl;
// 	  std::string mvtxav1("av_1");
// 	  std::string mvtxNodeName = mvtxNode->GetName();
// 
// 	  /// We only want the av_1 nodes
// 	  if(mvtxNodeName.compare(0, mvtxav1.length(), mvtxav1) == 0)
// 	    getMvtxKeyFromNode(mvtxNode);
// 	}
//     }
//     else if (node_str.compare(0, intt.length(), intt) == 0)  // is it in the INTT?
//     {
//       // We do not want the "ladderext" nodes
//       if (node_str.compare(0, intt_ext.length(), intt_ext) == 0)
//         continue;
// 
//       if (Verbosity() > 2)
// 	std::cout << " node " << node->GetName() << " is in the INTT"
// 		  << std::endl;
//       getInttKeyFromNode(node);
//     }
//     /// Put placeholders for the TPC and MMs. Because we modify the geometry
//     /// within TGeoVolume, we don't need a mapping to the TGeoNode
//     else if (node_str.compare(0, tpc.length(), tpc) == 0)  // is it in the TPC?
//       {
// 	if(Verbosity() > 2)
// 	  std::cout << " node " << node->GetName()
// 		    << " is in the TPC " << std::endl;
//       }
//     else if (node_str.compare(0, micromegas.length(), micromegas) == 0)  // is it in the Micromegas?
//       {
// 	if(Verbosity() > 2)
// 	  std::cout << " node " << node->GetName()
// 		    << " is in the MMs " << std::endl;
//       }
//     else
//       continue;
// 
//     bool print_sensor_paths = false;  // normally false
//     if (print_sensor_paths)
//     {
//       // Descends the node tree to find the active silicon nodes - used for information only
//       std::cout << " Top Node is " << node->GetName() << " volume name is " << node->GetVolume()->GetName() << std::endl;
// 
//       int nmax_print = 20;
//       isActive(node, nmax_print);
//     }
//   }
// }

void MakeActsGeometry::getInttKeyFromNode(TGeoNode *gnode)
{
  int layer = -1;       // sPHENIX layer number
  int itype = -1;       // specifies inner (0) or outer (1) sensor
  int ladder_phi = -1;  // copy number of ladder in phi
  int zposneg = -1;     // specifies positive (1) or negative (0) z
  int ladder_z = -1;    // 0-3, from most negative z to most positive

  std::string s = gnode->GetName();
  std::string delimiter = "_";
  std::string posz("posz");
  std::string negz("negz");

  size_t pos = 0;
  std::string token;

  int counter = 0;
  while ((pos = s.find(delimiter)) != std::string::npos)
  {
    token = s.substr(0, pos);
 
    s.erase(0, pos + delimiter.length());
    if (counter == 1)
      layer = std::atoi(token.c_str()) + 3;
    if (counter == 2)
      itype = std::atoi(token.c_str());
    if (counter == 3)
    {
      ladder_phi = std::atoi(token.c_str());
      if (s.compare(0, negz.length(), negz) == 0) zposneg = 0;
      if (s.compare(0, posz.length(), posz) == 0) zposneg = 1;
    }
    counter++;
  }

  ladder_z = itype + zposneg * 2;

  // The active sensor is a daughter of gnode
  int ndaught = gnode->GetNdaughters();
  if (ndaught == 0)
  {
    std::cout << PHWHERE << "OOPS: Did not find INTT sensor! Quit." 
	      << std::endl;
    exit(1);
  }

  std::string intt_refactive("siactive");
  TGeoNode *sensor_node = 0;
  for (int i = 0; i < ndaught; ++i)
  {
    std::string node_str = gnode->GetDaughter(i)->GetName();

    if (node_str.compare(0, intt_refactive.length(), intt_refactive) == 0)
    {
      sensor_node = gnode->GetDaughter(i);
      break;
    }
  }

  // unique key identifying this sensor
  int crossing = 0;
  TrkrDefs::hitsetkey node_key = InttDefs::genHitSetKey(layer, ladder_z, ladder_phi, crossing);

  std::pair<TrkrDefs::hitsetkey, TGeoNode *> tmp = std::make_pair(
					         node_key, sensor_node);
  m_clusterNodeMap.insert(tmp);

  if (Verbosity() > 1)
    std::cout << " INTT layer " << layer << " ladder_phi " << ladder_phi 
	      << " itype " << itype << " zposneg " << zposneg 
	      << " ladder_z " << ladder_z << " name " 
	      << sensor_node->GetName() << std::endl;

  return;
}
void MakeActsGeometry::getTpcKeyFromNode(TGeoNode * /*gnode*/)
{

}
void MakeActsGeometry::getMvtxKeyFromNode(TGeoNode *gnode)
{
  int counter = 0;
  int impr = -1;  // stave number, 1-48 in TGeo
  int layer = -1;
  int stave = -1;  // derived from impr
  int chip = -1;   // 9 chips per stave

  std::string s = gnode->GetName();
  std::string delimiter = "_";

  size_t pos = 0;
  std::string token;

  while ((pos = s.find(delimiter)) != std::string::npos)
  {
    token = s.substr(0, pos);
    //std::cout << token << std::endl;
    s.erase(0, pos + delimiter.length());
    if (counter == 3)
      impr = std::atoi(token.c_str());

    counter++;
  }

  // extract layer and stave info from impr
  // int staves_in_layer[3] = {12, 16, 20};
  // note - impr stave count starts from 1, not 0, but TrkrCluster counting starts from 0, so we reduce it by 1 here
  impr -= 1;

  if (impr < 12)
  {
    layer = 0;
    stave = impr;
  }
  else if (impr > 11 && impr < 28)
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
  for (int i = 0; i < mnd; ++i)
  {
    std::string dstr = module_node->GetDaughter(i)->GetName();
    if (dstr.compare(0, mvtx_chip.length(), mvtx_chip) == 0)
    {
      if (Verbosity() > 3)
        std::cout << "Found MVTX layer " << layer << " stave " << stave 
		  << " chip  " << i << " with node name " 
		  << module_node->GetDaughter(i)->GetName() << std::endl;

      // Make key for this chip
      unsigned int strobe = 0;
      TrkrDefs::hitsetkey node_key = MvtxDefs::genHitSetKey(layer, 
							    stave, i, strobe);

      // add sensor node to map
      TGeoNode *sensor_node = module_node->GetDaughter(i)->GetDaughter(0);
      std::pair<TrkrDefs::hitsetkey, TGeoNode *> tmp = std::make_pair(
						       node_key, sensor_node);
      m_clusterNodeMap.insert(tmp);

      if (Verbosity() > 3)
        std::cout << " MVTX layer " << layer << " stave " << stave 
		  << " chip " << chip << " name " 
		  << sensor_node->GetName() << std::endl;
    }
  }

  return;
}

// void MakeActsGeometry::isActive(TGeoNode *gnode, int nmax_print)
// {
//   // Not used in analysis: diagnostic, for looking at the node tree only.
//   // Recursively searches gnode for silicon sensors, prints out heirarchy
// 
//   std::string node_str = gnode->GetName();
// 
//   std::string intt_refactive("siactive");
//   std::string mvtx_refactive("MVTXSensor");
//   std::string tpc_refactive("tpc_gas_measurement");
//   std::string micromegas_refactive("MICROMEGAS_55");
// 
//   if (node_str.compare(0, intt_refactive.length(), intt_refactive) == 0)
//   {
//     std::cout << "          ******* Found INTT active volume,  node is "
// 	      << gnode->GetName()
// 	      << " volume name is " << gnode->GetVolume()->GetName()
// 	      << std::endl;
// 
//     return;
//   }
//   else if (node_str.compare(0, mvtx_refactive.length(), mvtx_refactive) == 0)
//   {
//     std::cout << "        ******* Found MVTX active volume,  node is "
// 	      << gnode->GetName()
// 	      << " volume name is " << gnode->GetVolume()->GetName()
// 	      << std::endl;
// 
//      return;
//   }
//   else if (node_str.compare(0, tpc_refactive.length(), tpc_refactive) == 0)
//   {
//     if(nprint_tpc < nmax_print)
//       {
// 	std::cout << "     ******* Found TPC  active volume,  node is "
// 		  << gnode->GetName()
// 		  << " volume name is " << gnode->GetVolume()->GetName()
// 		  << std::endl;
//       }
//     nprint_tpc++;
// 
//     return;
//   }
//   else if (node_str.compare(0, micromegas_refactive.length(), micromegas_refactive) == 0)
//   {
//     std::cout << "     ******* Found Micromegas  active volume,  node is "
// 	      << gnode->GetName()
// 	      << " volume name is " << gnode->GetVolume()->GetName()
// 	      << std::endl;
// 
//     return;
//   }
//   else
//     {
//       if(nprint_tpc < nmax_print)
// 	{
// 	  std::cout << "          ******* Found  node "
// 		    << gnode->GetName()
// 		    << " volume name is " << gnode->GetVolume()->GetName()
// 		    << std::endl;
// 	}
//       nprint_tpc++;
// 
//       return;
//     }
// 
// 
//   int ndaught = gnode->GetNdaughters();
// 
//   for (int i = 0; i < ndaught; ++i)
//   {
//     std::cout << "     " << gnode->GetVolume()->GetName()
// 	      << "  daughter " << i << " has name "
// 	      << gnode->GetDaughter(i)->GetVolume()->GetName() << std::endl;
// 
//     isActive(gnode->GetDaughter(i), nmax_print);
//   }
// }


void MakeActsGeometry::setPlanarSurfaceDivisions()
{
  /// These are arbitrary tpc subdivisions, and may change
  /// Setup how TPC boxes will be built for Acts::Surfaces
  m_surfStepZ = (m_maxSurfZ - m_minSurfZ) / (double) m_nSurfZ;
  m_moduleStepPhi = 2.0 * M_PI / 12.0;
  m_modulePhiStart = -M_PI;
  m_surfStepPhi = 2.0 * M_PI / (double) (m_nSurfPhi * m_nTpcModulesPerLayer);

  int layer=0;
  PHG4TpcCylinderGeomContainer::ConstRange layerrange = m_geomContainerTpc->get_begin_end();
  for (PHG4TpcCylinderGeomContainer::ConstIterator layeriter = layerrange.first;
       layeriter != layerrange.second;
       ++layeriter)
  {
    m_layerRadius[layer] = layeriter->second->get_radius();
    m_layerThickness[layer] = layeriter->second->get_thickness();
    layer++;
  }

}

int MakeActsGeometry::createNodes(PHCompositeNode *topNode)
{

  PHNodeIterator iter(topNode);
  /// Get the DST Node
  PHCompositeNode *parNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "PAR"));
  
  /// Check that it is there
  if (!parNode)
    {
      std::cout << "PAR Node missing, creating it" << std::endl;
      parNode = new PHCompositeNode("PAR");
      topNode->addNode(parNode);
    }
  
  PHNodeIterator pariter(parNode);
  /// Get the tracking subnode
  PHCompositeNode *svtxNode = dynamic_cast<PHCompositeNode *>(pariter.findFirst("PHCompositeNode", "SVTX"));
  
  /// Check that it is there
  if (!svtxNode)
    {
      svtxNode = new PHCompositeNode("SVTX");
      parNode->addNode(svtxNode);
    }

  m_actsGeometry = findNode::getClass<ActsGeometry>(topNode, "ActsGeometry");
  if(!m_actsGeometry)
    {
      m_actsGeometry = new ActsGeometry();
      PHDataNode<ActsGeometry> *tGeoNode 
	= new PHDataNode<ActsGeometry>(m_actsGeometry, "ActsGeometry");
      svtxNode->addNode(tGeoNode);
    }
  
  return Fun4AllReturnCodes::EVENT_OK;
}

/*
 * GetNodes():
 *  Get all the all the required nodes off the node tree
 */
int MakeActsGeometry::getNodes(PHCompositeNode *topNode)
{
  m_geoManager = PHGeomUtility::GetTGeoManager(topNode);
  if (!m_geoManager)
  {
    std::cout << PHWHERE << " Did not find TGeoManager, quit! " 
	      << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  m_geomContainerMvtx = findNode::getClass<
      PHG4CylinderGeomContainer>(topNode, "CYLINDERGEOM_MVTX");
  if (!m_geomContainerMvtx)
  {
    std::cout << PHWHERE 
	      << " CYLINDERGEOM_MVTX  node not found on node tree"
	      << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  m_geomContainerTpc =
      findNode::getClass<PHG4TpcCylinderGeomContainer>(topNode, "CYLINDERCELLGEOM_SVTX");
  if (!m_geomContainerTpc)
  {
    std::cout << PHWHERE 
	      << "ERROR: Can't find node CYLINDERCELLGEOM_SVTX" 
	      << std::endl;
    topNode->print();
    auto se = Fun4AllServer::instance();
    se->Print();
    return Fun4AllReturnCodes::ABORTRUN;
  }

  m_geomContainerIntt = findNode::getClass<
      PHG4CylinderGeomContainer>(topNode, "CYLINDERGEOM_INTT");
  if (!m_geomContainerIntt)
  {
    std::cout << PHWHERE 
	      << " CYLINDERGEOM_INTT  node not found on node tree"
	      << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  // load micromegas geometry
  // do not abort if not found
  m_geomContainerMicromegas = findNode::getClass<PHG4CylinderGeomContainer>(topNode, "CYLINDERGEOM_MICROMEGAS_FULL");
  if (!m_geomContainerMicromegas)
  {
    std::cout << PHWHERE 
	      << " CYLINDERGEOM_MICROMEGAS_FULL  node not found on node tree"
	      << std::endl;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

