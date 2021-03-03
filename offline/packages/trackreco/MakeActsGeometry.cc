/*!
 *  \file		MakeActsGeometry.C
 *  \brief		Refit SvtxTracks with PHActs.
 *  \details	Refit SvtxTracks with PHActs.
 *  \author	        Tony Frawley <afrawley@fsu.edu>
 */

#include "MakeActsGeometry.h"

#include <trackbase/TrkrDefs.h>

#include <intt/CylinderGeomIntt.h>
#include <intt/InttDefs.h>

#include <mvtx/CylinderGeom_Mvtx.h>
#include <mvtx/MvtxDefs.h>

#include <tpc/TpcDefs.h>

#include <micromegas/MicromegasDefs.h>

#include <g4detectors/PHG4CylinderCellGeom.h>
#include <g4detectors/PHG4CylinderCellGeomContainer.h>
#include <g4detectors/PHG4CylinderGeom.h>  // for PHG4CylinderGeom
#include <g4detectors/PHG4CylinderGeomContainer.h>

#include <phgeom/PHGeomUtility.h>
#include <phgeom/PHGeomIOTGeo.h>
#include <phgeom/PHGeomTGeo.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHDataNode.h>
#include <phool/PHNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>
#include <phool/getClass.h>
#include <phool/phool.h>

#include <Acts/EventData/TrackParameters.hpp>
#include <Acts/Geometry/GeometryContext.hpp>
#include <Acts/Geometry/TrackingVolume.hpp>
#include <Acts/MagneticField/MagneticFieldContext.hpp>
#include <Acts/Surfaces/PerigeeSurface.hpp>
#include <Acts/Surfaces/PlaneSurface.hpp>
#include <Acts/Surfaces/Surface.hpp>
#include <Acts/Utilities/CalibrationContext.hpp>

#include <ActsExamples/Detector/IBaseDetector.hpp>
#include <ActsExamples/EventData/Track.hpp>
#include <ActsExamples/Framework/AlgorithmContext.hpp>
#include <ActsExamples/Framework/IContextDecorator.hpp>
#include <ActsExamples/Framework/WhiteBoard.hpp>
#include <ActsExamples/Geometry/CommonGeometry.hpp>
#include <ActsExamples/Options/CommonOptions.hpp>
#include <ActsExamples/Plugins/Obj/ObjWriterOptions.hpp>
#include <ActsExamples/Utilities/Options.hpp>

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

MakeActsGeometry::MakeActsGeometry(const std::string &name)
: SubsysReco(name)
{
  setPlanarSurfaceDivisions();
  nprint_tpc = 0;
}

MakeActsGeometry::~MakeActsGeometry()
{}

int MakeActsGeometry::Init(PHCompositeNode *topNode)
{  
  return Fun4AllReturnCodes::EVENT_OK;
}
int MakeActsGeometry::InitRun(PHCompositeNode *topNode)
{

  if(buildAllGeometry(topNode) != Fun4AllReturnCodes::EVENT_OK)
    return Fun4AllReturnCodes::ABORTEVENT;

  /// Set the actsGeometry struct to be put on the node tree
  m_actsGeometry->tGeometry = m_tGeometry;
  m_actsGeometry->magField = m_magneticField;
  m_actsGeometry->calibContext = m_calibContext;
  m_actsGeometry->magFieldContext = m_magFieldContext;
  m_actsGeometry->geoContext = m_geoCtxt;

  /// Same for the surface maps
  m_surfMaps->siliconSurfaceMap = m_clusterSurfaceMapSilicon;
  m_surfMaps->tpcSurfaceMap = m_clusterSurfaceMapTpcEdit;
  m_surfMaps->tGeoNodeMap = m_clusterNodeMap;
  m_surfMaps->mmSurfaceMap = m_clusterSurfaceMapMmEdit;

  return Fun4AllReturnCodes::EVENT_OK;

}

int MakeActsGeometry::process_event(PHCompositeNode *topNode)
{
  return Fun4AllReturnCodes::EVENT_OK;
}
int MakeActsGeometry::End(PHCompositeNode *topNode)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

int MakeActsGeometry::buildAllGeometry(PHCompositeNode *topNode)
{

  /// Add the TPC surfaces to the copy of the TGeoManager. 
  // this also adds the micromegas surfaces
  // Do this before anything else, so that the geometry is finalized
  editTPCGeometry(topNode);

  if(getNodes(topNode) != Fun4AllReturnCodes::EVENT_OK)
    return Fun4AllReturnCodes::ABORTEVENT;  
  
  if(createNodes(topNode) != Fun4AllReturnCodes::EVENT_OK)
    return Fun4AllReturnCodes::ABORTEVENT;

  /// Run Acts layer builder
  buildActsSurfaces();

  /// Create a map of sensor TGeoNode pointers using the TrkrDefs:: hitsetkey as the key
  makeTGeoNodeMap(topNode);

  /// Export the new geometry to a root file for examination
  if(Verbosity() > 3)
    {
      PHGeomUtility::ExportGeomtry(topNode, "sPHENIXActsGeom.root"); 
      PHGeomUtility::ExportGeomtry(topNode, "sPHENIXActsGeom.gdml");
    }

  return Fun4AllReturnCodes::EVENT_OK;
}

void MakeActsGeometry::editTPCGeometry(PHCompositeNode *topNode)
{
  // Because we reset and rebuild the geomNode, we do edits of the TPC and Micromegas geometry in the same module
  
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
  
  TGeoVolume *World_vol = geoManager->GetTopVolume();

  // TPC geometry edits
  //===============

  TGeoNode *tpc_envelope_node = nullptr;
  TGeoNode *tpc_gas_north_node = nullptr;

  // find tpc north gas volume at path of World*/tpc_envelope*
  if (Verbosity())
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
  if (Verbosity())
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

  if (Verbosity())
  {
    std::cout << "EditTPCGeometry - gas volume: ";
    tpc_gas_north_vol->Print();
  }

  // adds surfaces to the underlying volume, so both north and south placements get them
  addActsTpcSurfaces(tpc_gas_north_vol, geoManager);

  // Micromegas geometry edits
  // These will only be made if the Micromegas nodes are found in the node tree
  //====================
  // The micromegas detectors have both layers in the same tile. The inner and outer sides are mirrored
  // The detector details are (printed from Init method), where the thickness corresponds to the drift volume:
  // layer 55: Phi segmented, radius: 82.2565 cm, thickness: 0.3 cm, zmin: -110cm, zmax: 110cm, pitch: 0.0976562cm
  // layer: 55 volume: MICROMEGAS_55_Gas2_inner_phys
  // layer 56: radius: z segmented, 82.6998 cm, thickness: 0.3 cm, zmin: -110cm, zmax: 110cm, pitch: 0.195312cm
  // layer: 56 volume: MICROMEGAS_55_Gas2_outer_phys

  TGeoNode *micromegas_envelope_node = nullptr;

  for (int i = 0; i < World_vol->GetNdaughters(); i++)
  {
    TString node_name = World_vol->GetNode(i)->GetName();

    if (node_name.BeginsWith("MICROMEGAS"))
    {
      if (Verbosity())
        std::cout << "EditTPCGeometry - found micromegas node " << node_name << std::endl;

      micromegas_envelope_node = World_vol->GetNode(i);
      break;
    }
  }

  if(micromegas_envelope_node)
    {
      /// If the node was found, we're building the MMs
      m_buildMMs = true;

      TGeoVolume *micromegas_envelope_vol = micromegas_envelope_node->GetVolume();
      assert(micromegas_envelope_vol);
      
      // Get inner and outer volume and edit them
      for (int i = 0; i < micromegas_envelope_vol->GetNdaughters(); i++)
	{
	  TString node_name = micromegas_envelope_vol->GetNode(i)->GetName();
	  
	  // this gets both inner and outer
	  if (node_name.BeginsWith("MICROMEGAS_55_Gas2"))
	    {
	      if (Verbosity())
		std::cout << "EditTPCGeometry - found Micromegas node " << node_name << std::endl;
	      
	      TGeoNode *micromegas_node = nullptr;
	      micromegas_node = micromegas_envelope_vol->GetNode(i);
	      
	      int mm_layer;
	      if( node_name.BeginsWith("MICROMEGAS_55_Gas2_inner") )
		mm_layer = 0;
	      else
		mm_layer = 1;
	      
	      TGeoVolume *micromegas_vol = micromegas_node->GetVolume();
	      assert(micromegas_vol);
	      
	      addActsMicromegasSurfaces(mm_layer, micromegas_vol, geoManager);
	    }
	}
    }
  
  // done
  geoManager->CloseGeometry();
  
  // save the edited geometry to DST persistent IO node for downstream DST files
  PHGeomUtility::UpdateIONode(topNode);

}

void MakeActsGeometry::addActsMicromegasSurfaces(int mm_layer, 
						 TGeoVolume *micromegas_vol, 
						 TGeoManager *geoManager)
{
  // The input micromegas_vol is either the inner (mm_layer 0) or outer (mm_layer 1) drift volume

  // The surfaces for both inner (phi segmented) and outer (z segmented) Micromegas detectors are long in z and small in phi
  // ---- We cannot approximate a cylinder with surfaces that are long in phi.

  TGeoMedium *micromegas_medium = micromegas_vol->GetMedium();
  assert(micromegas_medium);

  TGeoVolume *micromegas_measurement_vol;

  // we use the same phi steps as for the TPC, for the same reasons
  // There are 12*12 phi locations for the surfaces, just as for the TPC
  double tan_half_phi = tan(m_surfStepPhi / 2.0);
  double box_z_length = 110.0 * 2 - 0.01;


  // make a box for this layer
  char bname[500];
  sprintf(bname,"micromegas_measurement_%i",mm_layer);
  
  // Because we use a box, not a section of a cylinder, we need this to prevent overlaps
  // set the nominal r*phi dimension of the box so they just touch at the inner edge when placed 
  double box_r_phi = 2.0 * tan_half_phi * (m_mmLayerRadius[mm_layer] - m_mmLayerThickness[mm_layer] / 2.0 ) - 0.0001;
  
  
  double box_thickness =  m_mmLayerThickness[mm_layer] -  0.1;  // makes it 2 mm thick inside 3 mm thick cylinder
  
  micromegas_measurement_vol = geoManager->MakeBox(bname, micromegas_medium, 
							   box_thickness / 2.0, 
							   box_r_phi / 2.0, 
							   box_z_length / 2.0);
  
  micromegas_measurement_vol->SetLineColor(kBlack);
  micromegas_measurement_vol->SetFillColor(kYellow);
  micromegas_measurement_vol->SetVisibility(kTRUE);
  
  if(Verbosity() > 1)
    {
      std::cout << Verbosity() << " Made box for Micromegas layer " 
		<< mm_layer << " with dx " << box_thickness << " dy " 
		<< box_r_phi << " ref arc " 
		<< m_surfStepPhi * m_mmLayerRadius[mm_layer] << " dz " 
		<< box_z_length << std::endl;
      micromegas_measurement_vol->Print();
    }      
  
  
  // place the boxes inside the micromegas drift volume cylinders

  int copy = 0;	      
  for (unsigned int iz = 0; iz < m_nSurfZ; ++iz)
    {
      // The (half) micromegas volume is 2 * 110 cm long and is symmetric around (x,y,z) = (0,0,0) in its frame
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
	      
	      copy++;
		  
	      // place copies of the gas volume to fill up the layer
	      double x_center = m_mmLayerRadius[mm_layer] * cos(phi_center);
	      double y_center = m_mmLayerRadius[mm_layer] * sin(phi_center);
	      
	      char rot_name[500];
	      sprintf(rot_name,"micromegas_rotation_%i", copy);
	      TGeoCombiTrans *micromegas_measurement_location = 
		new TGeoCombiTrans(x_center, y_center, z_center,
				   new TGeoRotation(rot_name,
						    phi_center_degrees, 
						    0, 0));
	      
	      micromegas_vol->AddNode(micromegas_measurement_vol, copy, 
				      micromegas_measurement_location);
	      
	      if(Verbosity() > 10) 
		{
		  std::cout << " Made copy " << copy << mm_layer << " iphi " 
			    << iphi << std::endl;
		  std::cout << "    x_center " << x_center 
			    << " y_center " << y_center << " z_center " 
			    << z_center << " phi_center_degrees " 
			    << phi_center_degrees << std::endl;
		}
	    }
	}
    }
}


void MakeActsGeometry::addActsTpcSurfaces(TGeoVolume *tpc_gas_vol, 
					  TGeoManager *geoManager)
{
  TGeoMedium *tpc_gas_medium = tpc_gas_vol->GetMedium();
  assert(tpc_gas_medium);

  TGeoVolume *tpc_gas_measurement_vol[m_nTpcLayers];
  double tan_half_phi = tan(m_surfStepPhi / 2.0);
  
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

      if(Verbosity() > 30)
	{
	  std::cout << Verbosity() << " Made box for layer " << ilayer 
		    << " with dx " << m_layerThickness[ilayer] << " dy " 
		    << box_r_phi << " ref arc " 
		    << m_surfStepPhi * m_layerRadius[ilayer] << " dz " 
		    << m_surfStepZ << std::endl;
	  tpc_gas_measurement_vol[ilayer]->Print();
	}
    }

  int copy = 0;	      
  for (unsigned int iz = 0; iz < m_nSurfZ; ++iz)
    {
      // The (half) tpc gas volume is 105.5 cm long and is symmetric around (x,y,z) = (0,0,0) in its frame
      double z_center = -105.5/2.0 + m_surfStepZ / 2.0 + (double) iz * m_surfStepZ;
      
      for (unsigned int imod = 0; imod < m_nTpcModulesPerLayer; ++imod)
	{
	  for (unsigned int iphi = 0; iphi < m_nSurfPhi; ++iphi)
	    {

	      double min_phi = m_modulePhiStart + 
		(double) imod * m_moduleStepPhi + 
		(double) iphi * m_surfStepPhi;
	      double phi_center = min_phi + m_surfStepPhi / 2.0;
	      double phi_center_degrees = phi_center * 180.0 / M_PI;
	      
	      for (unsigned int ilayer = 0; ilayer < m_nTpcLayers; ++ilayer)
		{
		  copy++;
		  
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
		  
		  if(Verbosity() > 30 && ilayer == 30) 
		    {
		      std::cout << " Made copy " << copy << " iz " << iz 
				<< " imod " << imod << " ilayer " << ilayer
				<< " iphi " << iphi << std::endl;
		      std::cout << "    x_center " << x_center 
				<< " y_center " << y_center 
				<< " z_center " << z_center 
				<< " phi_center_degrees " << phi_center_degrees 
				<< std::endl;
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
  const int argc = 14;
  char* arg[argc];
 
  if(Verbosity() > 0)
    std::cout << PHWHERE << "Magnetic field " << m_magField 
	      << " with rescale " << m_magFieldRescale << std::endl;
  
  /// If the 2d fieldmap is provided, for now we just assume a 1.4T
  /// field (which will be properly scaled by 1.4/1.5) from magFieldRescale
  if(m_magField.find(".root") != std::string::npos)
    {
      m_magField = "1.5";
    }
  
  /// The acts constant field and geant field are backwards in convention
  m_magFieldRescale *= -1;

  if(Verbosity() > 0)
    std::cout << "Mag field now " << m_magField << " with rescale "
	      << m_magFieldRescale << std::endl;

  std::string responseFile = "tgeo-sphenix.response";
  std::string materialFile = "sphenix-material.json";
  if(m_buildMMs)
    materialFile = "sphenix-mm-material.json";

  /// Check to see if files exist locally - if not, use defaults
  std::ifstream file;

  file.open(responseFile);
  if(!file)
    {
      if(m_buildMMs)
	responseFile = std::string(getenv("OFFLINE_MAIN")) +
	  std::string("/share/tgeo-sphenix-mms.response");
      else
	responseFile = std::string(getenv("OFFLINE_MAIN")) +
	  std::string("/share/tgeo-sphenix.response");
    }
    
  file.open(materialFile);
  if(!file)
    {
      std::cout << materialFile << " not found locally, use repo version" << std::endl;
      materialFile = std::string(getenv("CALIBRATIONROOT")) +
	std::string("/ACTS/sphenix-material.json");
    }
  
  if(Verbosity() > -1)
    {
      std::cout << "using Acts material file : " << materialFile 
		<< std::endl;
      std::cout << "Using Acts TGeoResponse file : " << responseFile
		<< std::endl;
    }
  

  // Response file contains arguments necessary for geometry building
  const std::string argstr[argc]{
    "-n1", "-l0", 
      "--response-file",
      responseFile,
      "--bf-values","0","0",m_magField,
      "--bf-bscalor", std::to_string(m_magFieldRescale),
      "--mat-input-type","file",
      "--mat-input-file",
      materialFile
      };

  // Set vector of chars to arguments needed
  for (int i = 0; i < argc; ++i)
    {
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

void MakeActsGeometry::makeGeometry(int argc, char* argv[], 
				    ActsExamples::IBaseDetector &detector)
{
  
  /// setup and parse options
  auto desc = ActsExamples::Options::makeDefaultOptions();
  ActsExamples::Options::addGeometryOptions(desc);
  ActsExamples::Options::addMaterialOptions(desc);
  ActsExamples::Options::addObjWriterOptions(desc);
  ActsExamples::Options::addOutputOptions(desc);
  ActsExamples::Options::addBFieldOptions(desc);

  /// Add specific options for this geometry
  detector.addOptions(desc);
  auto vm = ActsExamples::Options::parse(desc, argc, argv);

  /// The geometry, material and decoration
  auto geometry = ActsExamples::Geometry::build(vm, detector);
  /// Geometry is a pair of (tgeoTrackingGeometry, tgeoContextDecorators)
  m_tGeometry = geometry.first;
  m_contextDecorators = geometry.second;

  m_magneticField = ActsExamples::Options::readBField(vm);

  size_t ievt = 0;
  size_t ialg = 0;

  /// Setup the event and algorithm context
  ActsExamples::WhiteBoard eventStore(
            Acts::getDefaultLogger("EventStore#" + std::to_string(ievt),
				   Acts::Logging::Level::INFO));

  /// The geometry context
  ActsExamples::AlgorithmContext context(ialg, ievt, eventStore);

  m_calibContext = context.calibContext;
  m_magFieldContext = context.magFieldContext;
  m_geoCtxt = context.geoContext;
    
  unpackVolumes();
  
  return;
}

void MakeActsGeometry::unpackVolumes()
{
  /// m_tGeometry is a TrackingGeometry pointer
  /// vol is a TrackingVolume pointer  
  auto vol = m_tGeometry->highestTrackingVolume();

  if(Verbosity() > 10 )
    std::cout << "Highest Tracking Volume is "
	      << vol->volumeName() << std::endl;

  /// volumeVector is a std::vector<TrackingVolumePtrs>
  auto volumeVector = vol->confinedVolumes()->arrayObjects();

  if(m_buildMMs)
    {
      auto mmBarrel = volumeVector.at(1);
      makeMmMapPairs(mmBarrel);
    }

  /// We have several volumes to walk through with the tpc and silicon
  auto firstVolumes = volumeVector.at(0)->confinedVolumes();
  auto topVolumesVector = firstVolumes->arrayObjects();
  
  if(Verbosity() > 10 )
    {
      for(long unsigned int i = 0; i<topVolumesVector.size(); i++)
	{
	  std::cout<< "TopVolume name: " 
		   << topVolumesVector.at(i)->volumeName() 
		   << std::endl;
	}
    }

  auto siliconVolumes = topVolumesVector.at(1)->confinedVolumes();
  auto siliconVolumesVector = siliconVolumes->arrayObjects();
  if(Verbosity() > 10 )
    {
      for(long unsigned int i =0; i<siliconVolumes->arrayObjects().size(); i++){
	std::cout << "SiliconVolumeName: " 
		  << siliconVolumes->arrayObjects().at(i)->volumeName()
		  << std::endl;
      }
    }

  /// Depending on whether or not the MMs are being built, the 
  /// Silicon and TPC volumes are packed differently
  /// This actually contains the silicon volumes
  if(!m_buildMMs)
    {
      auto mvtxVolumes = siliconVolumesVector.at(0);
      auto mvtxConfinedVolumes = mvtxVolumes->confinedVolumes();
      auto mvtxBarrel = mvtxConfinedVolumes->arrayObjects().at(1);

      makeMvtxMapPairs(mvtxBarrel);

      /// INTT only has one volume, so there is not an added volume extraction
      /// like for the MVTX
      auto inttVolume =  siliconVolumesVector.at(1);

      makeInttMapPairs(inttVolume);

      /// Same for the TPC - only one volume
      auto tpcVolume = volumeVector.at(1);
  
      makeTpcMapPairs(tpcVolume);
    }
  else
    {
      /// Additional layer unpacking if MMs were built
      auto nextSiliconVolumes = siliconVolumesVector.at(0)->confinedVolumes();
      auto siliconVolume = nextSiliconVolumes->arrayObjects().at(1)->confinedVolumes();
      
      auto mvtxVolumes = siliconVolume->arrayObjects().at(0);
      auto mvtxConfinedVolumes = mvtxVolumes->confinedVolumes();
      auto mvtxBarrel = mvtxConfinedVolumes->arrayObjects().at(1);
      makeMvtxMapPairs(mvtxBarrel);
      
      /// INTT only has one volume, so there is not an added volume extraction
      /// like for the MVTX
      auto inttVolume =  siliconVolume->arrayObjects().at(1);
      makeInttMapPairs(inttVolume);
      
      /// Same for the TPC - only one volume. Buried under silicon
      /// volume array
      auto tpcVolume = siliconVolumes->arrayObjects().at(1);

      makeTpcMapPairs(tpcVolume);
    }

  return;
}

void MakeActsGeometry::makeTpcMapPairs(TrackingVolumePtr &tpcVolume)
{
  if(Verbosity() > 1)
    std::cout << "Building TPC with " << tpcVolume->volumeName() << std::endl;
   
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

	  /// If there is already an entry for this hitsetkey, add the surface
	  /// to its corresponding vector
	  std::map<TrkrDefs::hitsetkey, std::vector<Surface>>::iterator mapIter;
	  mapIter = m_clusterSurfaceMapTpcEdit.find(hitsetkey);
	  
	  if(mapIter != m_clusterSurfaceMapTpcEdit.end())
	    {
	      mapIter->second.push_back(surf);
	    }
	  else
	    {
	      /// Otherwise make a new map entry
	      std::vector<Surface> dumvec;
	      dumvec.push_back(surf);
	      std::pair<TrkrDefs::hitsetkey, std::vector<Surface>> tmp = 
		std::make_pair(hitsetkey, dumvec);
	      m_clusterSurfaceMapTpcEdit.insert(tmp);
	    }
	  
	}
    }

}

void MakeActsGeometry::makeMmMapPairs(TrackingVolumePtr &mmVolume)
{

  if(Verbosity() > 1)
    std::cout << "Building MMs with " << mmVolume->volumeName() << std::endl;
    

  auto mmLayerArray = mmVolume->confinedLayers();
  auto mmLayerVector = mmLayerArray->arrayObjects();

  /// Need to unfold each layer that Acts builds
  for(unsigned int i = 0; i < mmLayerVector.size(); i++)
    {
      auto surfaceArray = mmLayerVector.at(i)->surfaceArray();
      if(surfaceArray == NULL){
	continue;
      }
      /// surfaceVector is a vector of surfaces corresponding to the micromegas layer
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

	  TrkrDefs::hitsetkey hitsetkey = getMmHitSetKeyFromCoords(world_center);

	  /// If there is already an entry for this hitsetkey, add the surface
	  /// to its corresponding vector
	  std::map<TrkrDefs::hitsetkey, std::vector<Surface>>::iterator mapIter;
	  mapIter = m_clusterSurfaceMapMmEdit.find(hitsetkey);
	  
	  if(mapIter != m_clusterSurfaceMapMmEdit.end())
	    {
	      mapIter->second.push_back(surf);
	    }
	  else
	    {
	      /// Otherwise make a new map entry
	      std::vector<Surface> dumvec;
	      dumvec.push_back(surf);
	      std::pair<TrkrDefs::hitsetkey, std::vector<Surface>> tmp = 
		std::make_pair(hitsetkey, dumvec);
	      m_clusterSurfaceMapMmEdit.insert(tmp);
	    }
	  
	}
    }
}

void MakeActsGeometry::makeInttMapPairs(TrackingVolumePtr &inttVolume)
{
  
  if(Verbosity() > 10)
    {
      std::cout << "intt volume name: "  << inttVolume->volumeName()
		<< std::endl;
    }

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
    {
      std::cout << "MVTX Barrel name to step surfaces through is " 
	      << mvtxVolume->volumeName() << std::endl;
    }
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
		<< phi_world*180/3.14159 << " layer " << layer 
		<< " readout_mod " << readout_mod << " side " << side 
		<< " hitsetkey " << hitset_key<< std::endl;
  
  return hitset_key;
}

TrkrDefs::hitsetkey MakeActsGeometry::getMmHitSetKeyFromCoords(std::vector<double> &world)
{
  // Look up micromegas surface index values from world position of surface center

  // layer
  unsigned int layer = 999;
  double layer_rad = sqrt(pow(world[0],2) + pow(world[1],2));
  for(unsigned int ilayer=0;ilayer<m_nMmLayers;++ilayer)
    {
      double ref_radius_low = 
	m_mmLayerRadius[ilayer] - m_mmLayerThickness[ilayer] / 2.0;
      double ref_radius_high = 
	m_mmLayerRadius[ilayer] + m_mmLayerThickness[ilayer] / 2.0;
      
      if(layer_rad >= ref_radius_low && layer_rad < ref_radius_high)
	{
	  layer =  m_mmLayerNumber[ilayer];
	  break;
	}
    }
  if(layer != m_mmLayerNumber[0] && layer != m_mmLayerNumber[1]) 
    {
      std::cout << PHWHERE 
		<< "Error: undefined layer, do nothing world =  " 
		<< world[0] << "  " << world[1] << "  " << world[2] 
		<< " layer " << layer << std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }

  // we have lumped all surfacrs into one layer, so we use a fake hitsetkey with the real layer number and segmentation type, but with a dummy value for the tile
  int tile = 0;   // assign all surfaces to tile 0
  MicromegasDefs::SegmentationType segtype;
  if(layer == m_mmLayerNumber[0])
    segtype  =  MicromegasDefs::SegmentationType::SEGMENTATION_PHI;
  else
    segtype = MicromegasDefs::SegmentationType::SEGMENTATION_Z;

  /// Get the surface key to find the surface from the map
  TrkrDefs::hitsetkey hitset_key = MicromegasDefs::genHitSetKey(layer, segtype, tile);

  if(Verbosity() > 3)
    std::cout << PHWHERE << "    micromegas layer " << layer 
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
  }

  unsigned int stave = 0;
  unsigned int chip = 0;
  layergeom->get_sensor_indices_from_world_coords(world, stave, chip);

  double check_pos[3] = {0, 0, 0};
  layergeom->find_sensor_center(stave, 0, 0, chip, check_pos);

  TrkrDefs::hitsetkey mvtx_hitsetkey = MvtxDefs::genHitSetKey(layer, stave, chip);

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

  double check_pos[3] = {0, 0, 0};
  layergeom->find_segment_center(segment_z_bin, segment_phi_bin, check_pos);

  TrkrDefs::hitsetkey intt_hitsetkey = InttDefs::genHitSetKey(layer, segment_z_bin, segment_phi_bin);

  return intt_hitsetkey;
}

void MakeActsGeometry::makeTGeoNodeMap(PHCompositeNode *topNode)
{
  // This is just a diagnostic method
  // it lets you list all of the nodes by setting print_sensors = true

  if (!m_geoManager)
  {
    std::cout << PHWHERE << " Did not find TGeoManager, quit! " << std::endl;
    return;
  }
  TGeoVolume *topVol = m_geoManager->GetTopVolume();
  TObjArray *nodeArray = topVol->GetNodes();

  TIter iObj(nodeArray);
  while (TObject *obj = iObj())
  {
    TGeoNode *node = dynamic_cast<TGeoNode *>(obj);
    std::string node_str = node->GetName();

    std::string mvtx("MVTX_Wrapper");
    std::string intt("ladder");
    std::string intt_ext("ladderext");
    std::string tpc("tpc_envelope");
    std::string micromegas("MICROMEGAS_55");

    if (node_str.compare(0, mvtx.length(), mvtx) == 0)  // is it in the MVTX?
    {
      if (Verbosity() > 2) 
	std::cout << " node " << node->GetName() << " is the MVTX wrapper" 
		  << std::endl;
  
      /// The Mvtx has an additional wrapper that needs to be unpacked
      TObjArray *mvtxArray = node->GetNodes();
      TIter mvtxObj(mvtxArray);
      while(TObject *mvtx = mvtxObj())
	{
	  TGeoNode *mvtxNode = dynamic_cast<TGeoNode *>(mvtx);
	  if(Verbosity() > 2)
	    std::cout << "mvtx node name is " << mvtxNode->GetName() 
		      << std::endl;
	  std::string mvtxav1("av_1");
	  std::string mvtxNodeName = mvtxNode->GetName();
	  
	  /// We only want the av_1 nodes
	  if(mvtxNodeName.compare(0, mvtxav1.length(), mvtxav1) == 0)
	    getMvtxKeyFromNode(mvtxNode);
	}
    }
    else if (node_str.compare(0, intt.length(), intt) == 0)  // is it in the INTT?
    {
      // We do not want the "ladderext" nodes
      if (node_str.compare(0, intt_ext.length(), intt_ext) == 0)
        continue;

      if (Verbosity() > 2) 
	std::cout << " node " << node->GetName() << " is in the INTT" 
		  << std::endl;
      getInttKeyFromNode(node);
    }
    /// Put placeholders for the TPC and MMs. Because we modify the geometry
    /// within TGeoVolume, we don't need a mapping to the TGeoNode
    else if (node_str.compare(0, tpc.length(), tpc) == 0)  // is it in the TPC?
      {
	if(Verbosity() > 2)
	  std::cout << " node " << node->GetName() 
		    << " is in the TPC " << std::endl;
      }
    else if (node_str.compare(0, micromegas.length(), micromegas) == 0)  // is it in the Micromegas?
      {
	if(Verbosity() > 2)
	  std::cout << " node " << node->GetName() 
		    << " is in the MMs " << std::endl;
      }
    else
      continue;

    bool print_sensor_paths = false;  // normally false
    if (print_sensor_paths)
    {
      // Descends the node tree to find the active silicon nodes - used for information only
      std::cout << " Top Node is " << node->GetName() << " volume name is " << node->GetVolume()->GetName() << std::endl;

      int nmax_print = 20;
      isActive(node, nmax_print);
    }
  }
}

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
  TrkrDefs::hitsetkey node_key = InttDefs::genHitSetKey(layer, ladder_z, 
							ladder_phi);

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
void MakeActsGeometry::getTpcKeyFromNode(TGeoNode *gnode)
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
      TrkrDefs::hitsetkey node_key = MvtxDefs::genHitSetKey(layer, 
							    stave, i);

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

void MakeActsGeometry::isActive(TGeoNode *gnode, int nmax_print)
{
  // Not used in analysis: diagnostic, for looking at the node tree only.
  // Recursively searches gnode for silicon sensors, prints out heirarchy

  std::string node_str = gnode->GetName();

  std::string intt_refactive("siactive");
  std::string mvtx_refactive("MVTXSensor");
  std::string tpc_refactive("tpc_gas_measurement");
  std::string micromegas_refactive("MICROMEGAS_55");

  if (node_str.compare(0, intt_refactive.length(), intt_refactive) == 0)
  {
    std::cout << "          ******* Found INTT active volume,  node is " 
	      << gnode->GetName()
	      << " volume name is " << gnode->GetVolume()->GetName() 
	      << std::endl;
        
    return;
  }
  else if (node_str.compare(0, mvtx_refactive.length(), mvtx_refactive) == 0)
  {
    std::cout << "        ******* Found MVTX active volume,  node is " 
	      << gnode->GetName()
	      << " volume name is " << gnode->GetVolume()->GetName() 
	      << std::endl;
 
     return;
  }
  else if (node_str.compare(0, tpc_refactive.length(), tpc_refactive) == 0)
  {
    if(nprint_tpc < nmax_print)
      {
	std::cout << "     ******* Found TPC  active volume,  node is " 
		  << gnode->GetName()
		  << " volume name is " << gnode->GetVolume()->GetName() 
		  << std::endl;
      }
    nprint_tpc++;

    return;
  }
  else if (node_str.compare(0, micromegas_refactive.length(), micromegas_refactive) == 0)
  {
    std::cout << "     ******* Found Micromegas  active volume,  node is "
	      << gnode->GetName()
	      << " volume name is " << gnode->GetVolume()->GetName() 
	      << std::endl;

    return;
  }
  else
    {
      if(nprint_tpc < nmax_print)
	{
	  std::cout << "          ******* Found  node " 
		    << gnode->GetName()
		    << " volume name is " << gnode->GetVolume()->GetName() 
		    << std::endl;
	}
      nprint_tpc++;
      
      return;      
    }


  int ndaught = gnode->GetNdaughters();

  for (int i = 0; i < ndaught; ++i)
  {
    std::cout << "     " << gnode->GetVolume()->GetName() 
	      << "  daughter " << i << " has name " 
	      << gnode->GetDaughter(i)->GetVolume()->GetName() << std::endl;
    
    isActive(gnode->GetDaughter(i), nmax_print);
  }
}


void MakeActsGeometry::setPlanarSurfaceDivisions()
{
  /// These are arbitrary tpc subdivisions, and may change
  /// Setup how TPC boxes will be built for Acts::Surfaces
  m_surfStepZ = (m_maxSurfZ - m_minSurfZ) / (double) m_nSurfZ;
  m_moduleStepPhi = 2.0 * M_PI / 12.0;
  m_modulePhiStart = -M_PI;
  m_surfStepPhi = 2.0 * M_PI / (double) (m_nSurfPhi * m_nTpcModulesPerLayer);
  for(unsigned int isector = 0; isector < 3; ++isector)
    {
      layer_thickness_sector[isector] = 
	(m_maxRadius[isector] - m_minRadius[isector]) / 16.0;

      for(unsigned int ilayer =0; ilayer < 16; ++ilayer)
	{
	  m_layerRadius[isector*16 + ilayer] = 
	    m_minRadius[isector] + layer_thickness_sector[isector] * 
	    (double) ilayer + layer_thickness_sector[isector] / 2.0;
	  
	  m_layerThickness[isector*16 + ilayer] = 
	    layer_thickness_sector[isector];
	}
    }
}

int MakeActsGeometry::createNodes(PHCompositeNode *topNode)
{
  PHNodeIterator iter(topNode);

  /// Get the DST Node
  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
  
  /// Check that it is there
  if (!dstNode)
    {
      std::cerr << "DST Node missing, quitting" << std::endl;
      throw std::runtime_error("failed to find DST node in PHActsSourceLinks::createNodes");
    }
  
  /// Get the tracking subnode
  PHCompositeNode *svtxNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "SVTX"));
  
  /// Check that it is there
  if (!svtxNode)
    {
      svtxNode = new PHCompositeNode("SVTX");
      dstNode->addNode(svtxNode);
    }

  m_surfMaps = findNode::getClass<ActsSurfaceMaps>(topNode,
						   "ActsSurfaceMaps");
  if(!m_surfMaps)
    {
      m_surfMaps = new ActsSurfaceMaps();
      PHDataNode<ActsSurfaceMaps> *node
	= new PHDataNode<ActsSurfaceMaps>(m_surfMaps, "ActsSurfaceMaps");
      svtxNode->addNode(node);

    }

  m_actsGeometry = findNode::getClass<ActsTrackingGeometry>(topNode,
							    "ActsTrackingGeometry");
  if(!m_actsGeometry)
    {
      m_actsGeometry = new ActsTrackingGeometry();
      PHDataNode<ActsTrackingGeometry> *tGeoNode 
	= new PHDataNode<ActsTrackingGeometry>(m_actsGeometry, "ActsTrackingGeometry");
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
      findNode::getClass<PHG4CylinderCellGeomContainer>(topNode, "CYLINDERCELLGEOM_SVTX");
  if (!m_geomContainerTpc)
  {
    std::cout << PHWHERE 
	      << "ERROR: Can't find node CYLINDERCELLGEOM_SVTX" 
	      << std::endl;
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

  return Fun4AllReturnCodes::EVENT_OK;
}
