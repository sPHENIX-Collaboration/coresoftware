#include "PHTpcResiduals.h"
#include "TpcSpaceChargeMatrixContainerv1.h"

#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>
#include <phool/phool.h>
#include <phool/PHDataNode.h>
#include <phool/PHNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>
#include <phool/PHTimer.h>

#include <trackbase/TrkrCluster.h>
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase_historic/SvtxTrack.h>
#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/SvtxVertexMap.h>
#include <micromegas/MicromegasDefs.h>

#include <Acts/Surfaces/PerigeeSurface.hpp>
#include <Acts/Geometry/GeometryIdentifier.hpp>
#include <Acts/MagneticField/ConstantBField.hpp>
#include <Acts/MagneticField/InterpolatedBFieldMap.hpp>
#include <Acts/MagneticField/SharedBField.hpp>
#include <Acts/Propagator/EigenStepper.hpp>
#include <Acts/Surfaces/Surface.hpp>

#include <ActsExamples/Plugins/BField/ScalableBField.hpp>
#include <ActsExamples/EventData/TrkrClusterSourceLink.hpp>


#include <cmath>
#include <TFile.h>
#include <TTree.h>

#include <iostream>
#include <sstream>

namespace 
{
  
  // square
  template<class T> inline constexpr T square( const T& x ) { return x*x; }

  // radius
  template<class T> T get_r( const T& x, const T& y ) { return std::sqrt( square(x) + square(y) ); }

  template<class T> T deltaPhi(const T& phi)
  {
    if (phi > M_PI) 
      return phi - 2. * M_PI;
    else if (phi <= -M_PI) 
      return phi + 2.* M_PI;
    else 
      return phi;
  }
}

PHTpcResiduals::PHTpcResiduals(const std::string &name)
  : SubsysReco(name)
  , m_matrix_container( new TpcSpaceChargeMatrixContainerv1 )
{}

int PHTpcResiduals::Init(PHCompositeNode */*topNode*/)
{
  if( m_savehistograms ) makeHistograms();
  return Fun4AllReturnCodes::EVENT_OK;
}

int PHTpcResiduals::InitRun(PHCompositeNode *topNode)
{
  if(getNodes(topNode) != Fun4AllReturnCodes::EVENT_OK)
    return Fun4AllReturnCodes::ABORTEVENT;

  if(createNodes(topNode) != Fun4AllReturnCodes::EVENT_OK)
    return Fun4AllReturnCodes::ABORTEVENT;

  return Fun4AllReturnCodes::EVENT_OK;
}

int PHTpcResiduals::process_event(PHCompositeNode *topNode)
{
  if(Verbosity() > 1)
    std::cout <<"Starting PHTpcResiduals event " 
	      << m_event << std::endl;

  if(m_event % 1000 == 0)
    std::cout << "PHTpcResiduals processed " << m_event 
	      << " events" << std::endl;

  int returnVal = processTracks(topNode);

  if(Verbosity() > 1)
    std::cout <<"Finished PHTpcResiduals event " 
	      << m_event << std::endl;
  
  m_event++;

  return returnVal;
}

int PHTpcResiduals::End(PHCompositeNode */*topNode*/)
{
  std::cout << "PHTpcResiduals::End - writing matrices to " << m_outputfile << std::endl;

  if(Verbosity() > 0)
    std::cout << "Number of bad SL propagations " 
	      << m_nBadProps << std::endl;
      
  // save matrix container in output file
  if( m_matrix_container )
  {
    std::unique_ptr<TFile> outputfile( TFile::Open( m_outputfile.c_str(), "RECREATE" ) );
    outputfile->cd();
    m_matrix_container->Write( "TpcSpaceChargeMatrixContainer" );
  }

  return Fun4AllReturnCodes::EVENT_OK;
}


int PHTpcResiduals::processTracks(PHCompositeNode */*topNode*/)
{

  std::cout << "proto track size " << m_trackMap->size()
	    <<std::endl;

  for(auto &[trackKey, track] : *m_trackMap)
    {
      if(Verbosity() > 1)
	std::cout << "Processing track key " << trackKey
		  << std::endl;
      if(checkTrack(track))
	processTrack(track);
    }
  
  return Fun4AllReturnCodes::EVENT_OK;
}

bool PHTpcResiduals::checkTrack(SvtxTrack* track)
{
 
  if(track->get_pt() < 0.5)
    return false;

  if(Verbosity() > 2)
    std::cout << "Track has pt " << track->get_pt() << std::endl;

  int nMvtxHits = 0;
  int nInttHits = 0;
  int nMMHits = 0;

  for (SvtxTrack::ConstClusterKeyIter clusIter = track->begin_cluster_keys();
       clusIter != track->end_cluster_keys();
       ++clusIter)
    {
      auto key = *clusIter;
      auto trkrId = TrkrDefs::getTrkrId(key);
      if(trkrId == TrkrDefs::TrkrId::mvtxId)
	nMvtxHits++;
      else if(trkrId == TrkrDefs::TrkrId::inttId)
	nInttHits++;
      else if(trkrId == TrkrDefs::TrkrId::micromegasId)
	nMMHits++;
    }

  if(Verbosity() > 2)
    std::cout << "Number of mvtx/intt/MM hits "
	      << nMvtxHits << "/" << nInttHits << "/" 
	      << nMMHits << std::endl;

  // Require at least 2 hits in each detector
  if(nMvtxHits < 2 or nInttHits < 2 or nMMHits < 2)
    return false;

  return true;

}

SourceLink PHTpcResiduals::makeSourceLink(TrkrCluster* cluster)
{
  auto key = cluster->getClusKey();
  auto subsurfkey = cluster->getSubSurfKey();
      
  // Make a safety check for clusters that couldn't be attached
  // to a surface
  auto surf = getSurface(key, subsurfkey);
  if(!surf)
    return SourceLink();
  
  Acts::BoundVector loc = Acts::BoundVector::Zero();
  loc[Acts::eBoundLoc0] = cluster->getLocalX() * Acts::UnitConstants::cm;
  loc[Acts::eBoundLoc1] = cluster->getLocalY() * Acts::UnitConstants::cm;
  
  Acts::BoundMatrix cov = Acts::BoundMatrix::Zero();
  cov(Acts::eBoundLoc0, Acts::eBoundLoc0) = 
    cluster->getActsLocalError(0,0) * Acts::UnitConstants::cm2;
  cov(Acts::eBoundLoc0, Acts::eBoundLoc1) =
    cluster->getActsLocalError(0,1) * Acts::UnitConstants::cm2;
  cov(Acts::eBoundLoc1, Acts::eBoundLoc0) = 
    cluster->getActsLocalError(1,0) * Acts::UnitConstants::cm2;
  cov(Acts::eBoundLoc1, Acts::eBoundLoc1) = 
    cluster->getActsLocalError(1,1) * Acts::UnitConstants::cm2;
  
  SourceLink sl(key, surf, loc, cov);
  
  return sl;
}

//___________________________________________________________________________________
Surface PHTpcResiduals::getSurface(TrkrDefs::cluskey cluskey, TrkrDefs::subsurfkey surfkey)
{
  const auto trkrid = TrkrDefs::getTrkrId(cluskey);
  const auto hitsetkey = TrkrDefs::getHitSetKeyFromClusKey(cluskey);

  switch( trkrid )
  {
    case TrkrDefs::TrkrId::micromegasId: return getMMSurface( hitsetkey );
    case TrkrDefs::TrkrId::tpcId: return getTpcSurface(hitsetkey, surfkey);
    case TrkrDefs::TrkrId::mvtxId:
    case TrkrDefs::TrkrId::inttId:
    {
      return getSiliconSurface(hitsetkey);
    }
  }
  
  // unreachable
  return nullptr;
  
}

//___________________________________________________________________________________
Surface PHTpcResiduals::getSiliconSurface(TrkrDefs::hitsetkey hitsetkey)
{
  auto surfMap = m_surfMaps->siliconSurfaceMap;
  auto iter = surfMap.find(hitsetkey);
  if(iter != surfMap.end())
    {
      return iter->second;
    }
  
  // If it can't be found, return nullptr
  return nullptr;

}

//___________________________________________________________________________________
Surface PHTpcResiduals::getTpcSurface(TrkrDefs::hitsetkey hitsetkey, TrkrDefs::subsurfkey surfkey)
{
  const auto iter = m_surfMaps->tpcSurfaceMap.find(hitsetkey);
  if(iter != m_surfMaps->tpcSurfaceMap.end())
  {
    auto surfvec = iter->second;
    return surfvec.at(surfkey);
  }
  
  // If it can't be found, return nullptr to skip this cluster
  return nullptr;
}

//___________________________________________________________________________________
Surface PHTpcResiduals::getMMSurface(TrkrDefs::hitsetkey hitsetkey)
{
  const auto iter = m_surfMaps->mmSurfaceMap.find( hitsetkey );
  return (iter == m_surfMaps->mmSurfaceMap.end()) ? nullptr:iter->second;
}

//___________________________________________________________________________________
Acts::Vector3D PHTpcResiduals::getVertex(SvtxTrack *track)
{
  auto vertexId = track->get_vertex_id();
  const SvtxVertex* svtxVertex = m_vertexMap->get(vertexId);
  if(!svtxVertex)
    {
      if(Verbosity() > 0)
	std::cout << " Warning: No SvtxVertex for track found. Using (0,0,0)" 
		  << std::endl;
      return Acts::Vector3D(0,0,0);
    }

  Acts::Vector3D vertex(svtxVertex->get_x() * Acts::UnitConstants::cm, 
			svtxVertex->get_y() * Acts::UnitConstants::cm, 
			svtxVertex->get_z() * Acts::UnitConstants::cm);
  return vertex;
}
Acts::BoundTrackParameters PHTpcResiduals::makeTrackParams(SvtxTrack* track)
{
  Acts::Vector3D momentum(track->get_px(), 
			  track->get_py(), 
			  track->get_pz());
  auto vertex = getVertex(track);
  Acts::Vector4D actsFourPos(vertex(0), vertex(1), vertex(2),
			     10 * Acts::UnitConstants::ns);
  double trackQ = track->get_charge() * Acts::UnitConstants::e;
  double p = track->get_p();
  
  Acts::BoundSymMatrix cov;
  for(int i =0; i<6; i++)
    {
      for(int j =0; j<6; j++)
	{
	  cov(i,j) = track->get_acts_covariance(i, j);
	}
    }

  auto perigee = Acts::Surface::makeShared<Acts::PerigeeSurface>(
				        Acts::Vector3D(vertex(0), 
						       vertex(1), 
						       vertex(2)));

  Acts::BoundTrackParameters seed(perigee, m_tGeometry->geoContext,
				  actsFourPos, momentum,
				  p, trackQ, cov);
 
  return seed;

}

void PHTpcResiduals::processTrack(SvtxTrack* track)
{

  if(Verbosity() > 1)
    std::cout << "Propagating silicon+MM fit params momentum: " 
	      << track->get_p() << " and position " 
	      << track->get_x() << ", " << track->get_y() 
	      << ", " << track->get_z() << " cm "
	      << std::endl;

  auto trackParams = makeTrackParams(track);

  int initNBadProps = m_nBadProps;
  for (SvtxTrack::ConstClusterKeyIter clusIter = track->begin_cluster_keys();
       clusIter != track->end_cluster_keys();
       ++clusIter)
    {
      auto cluskey = *clusIter;
      
      // only propagate to tpc surfaces
      if(TrkrDefs::getTrkrId(cluskey) != TrkrDefs::TrkrId::tpcId)
	continue;;

      auto cluster = m_clusterContainer->findCluster(cluskey);

      auto sl = makeSourceLink(cluster);
      
      auto result = propagateTrackState(trackParams, sl);

      if(result.ok())
	{	  
	  auto trackStateParams = std::move(**result);
	  if(Verbosity() > 1)
	    std::cout << "Silicon+MM momentum : " 
		      << trackParams.momentum()
		      << std::endl
		      << "Propagator momentum : " 
		      << trackStateParams.momentum()
		      << std::endl;
	  

	  calculateTpcResiduals(trackStateParams, cluster);
	}
      else
	{
	  m_nBadProps++;
	
	  continue;
	}
    } 

  if(m_nBadProps > initNBadProps && Verbosity() > 1)
    std::cout << "Starting track params position/momentum: "
	      << trackParams.position(m_tGeometry->geoContext).transpose()
	      << std::endl << trackParams.momentum().transpose() 
	      << std::endl
	      << "Track params phi/eta " 
	      << std::atan2(trackParams.momentum().y(), 
			    trackParams.momentum().x())
	      << " and " 
	      << std::atanh(trackParams.momentum().z() / 
			    trackParams.momentum().norm())
	      << std::endl;
        
}

BoundTrackParamPtrResult PHTpcResiduals::propagateTrackState(
			   const Acts::BoundTrackParameters& params,
			   const SourceLink& sl)
{
  /*
  std::cout << "Propagating to geo id " << sl.referenceSurface().geometryId() << std::endl;
  if(sl.referenceSurface().associatedDetectorElement() != nullptr)
    std::cout << " which has associated detector element " << sl.referenceSurface().associatedDetectorElement()->thickness() << std::endl;
  */
  return std::visit([params, sl, this]
		    (auto && inputField) -> BoundTrackParamPtrResult {
      using InputMagneticField = 
	typename std::decay_t<decltype(inputField)>::element_type;
      using MagneticField      = Acts::SharedBField<InputMagneticField>;
      using Stepper            = Acts::EigenStepper<MagneticField>;
      using Propagator         = Acts::Propagator<Stepper>;

      MagneticField field(inputField);
      Stepper stepper(field);
      Propagator propagator(stepper);

      Acts::Logging::Level logLevel = Acts::Logging::FATAL;
      if(Verbosity() > 10)
	logLevel = Acts::Logging::VERBOSE;

      auto logger = Acts::getDefaultLogger("PHTpcResiduals", logLevel);
      
      Acts::PropagatorOptions<> options(m_tGeometry->geoContext,
					m_tGeometry->magFieldContext,
					Acts::LoggerWrapper{*logger});
     
      auto result = propagator.propagate(params, sl.referenceSurface(), 
					 options);
   
      if(result.ok())
	return std::move((*result).endParameters);
      else
	return result.error();
   },
     std::move(m_tGeometry->magField));

}
void PHTpcResiduals::calculateTpcResiduals(
  const Acts::BoundTrackParameters &params,
  const TrkrCluster* cluster)
{
  
  cluskey = cluster->getClusKey();
  // Get all the relevant information for residual calculation
  clusR = std::sqrt(square(cluster->getX()) + square(cluster->getY()));
  clusPhi = std::atan2(cluster->getY(), cluster->getX());
  clusZ = cluster->getZ();

  clusRPhiErr = cluster->getRPhiError();
  clusZErr = cluster->getZError();
 
  if(Verbosity() > 3)
    {
      std::cout << "cluster key is " << cluskey <<std::endl;
      std::cout << "Cluster r phi and z " << clusR << "  " 
		<< clusPhi << "+/-" << clusRPhiErr
		<<" and " << clusZ << "+/-" << clusZErr << std::endl;
    }
  
  if(clusRPhiErr < 0.015)
    return;
  if(clusZErr < 0.05)
    return;

  const auto globalStatePos = params.position(m_tGeometry->geoContext);
  const auto globalStateMom = params.momentum();
  const auto globalStateCov = *params.covariance();

  stateRPhiErr = sqrt(globalStateCov(Acts::eBoundLoc0,
				     Acts::eBoundLoc0))
    / Acts::UnitConstants::cm;
  stateZErr = sqrt(globalStateCov(Acts::eBoundLoc1,
				  Acts::eBoundLoc1))
    / Acts::UnitConstants::cm;
 
  stateZ = globalStatePos.z() / Acts::UnitConstants::cm;

  const auto globStateX = globalStatePos.x() / Acts::UnitConstants::cm;
  const auto globStateY = globalStatePos.y() / Acts::UnitConstants::cm;
  const auto globStateZ = stateZ;

  stateR = std::sqrt(square(globStateX) + square(globStateY));
  
  const auto dr = clusR - stateR;
  const auto trackDrDt = (globStateX * globalStateMom(0) +
			  globStateY * globalStateMom(1)) / stateR;
  const auto trackDxDr = globalStateMom(0) / trackDrDt;
  const auto trackDyDr = globalStateMom(1) / trackDrDt;
  const auto trackDzDr = globalStateMom(2) / trackDrDt;
  
  const auto trackX = globStateX + dr * trackDxDr;
  const auto trackY = globStateY + dr * trackDyDr;
  const auto trackZ = globStateZ + dr * trackDzDr;
  
  if(Verbosity() > 2)
    std::cout << "State Calculations: " << stateR << ", " 
	      << dr << ", " << trackDrDt << ", " << trackDxDr
	      << ", " << trackDyDr << ", " << trackDzDr
	      <<" , " << trackX << ", " << trackY << ", "
	      << trackZ << std::endl;

  statePhi = std::atan2(trackY, trackX);
  stateZ = trackZ;

  if(Verbosity() > 3)
    std::cout << "State r phi and z " 
	      << stateR
	      << "   " << statePhi << "+/-" << stateRPhiErr
	      << " and " << stateZ << "+/-" << stateZErr << std::endl;

  const auto erp = square(clusRPhiErr);
  const auto ez = square(clusZErr);

  const auto dPhi = clusPhi - statePhi;

  // Calculate residuals
  drphi = clusR * deltaPhi(dPhi);
  dz  = clusZ - stateZ;

  if(Verbosity() > 3)
    std::cout << "TPC residuals " << drphi << "   " << dz << std::endl;
  
  const auto trackEta 
    = std::atanh(params.momentum().z() / params.absoluteMomentum());
  const auto clusEta = std::atanh(clusZ / std::sqrt(
    square(cluster->getX()) +
    square(cluster->getY()) +
    square(cluster->getZ())));

  const auto trackPPhi = -params.momentum()(0) * std::sin(statePhi) +
    params.momentum()(1) * std::cos(statePhi);
  const auto trackPR = params.momentum()(0) * std::cos(statePhi) +
    params.momentum()(1) * std::sin(statePhi);
  const auto trackPZ    = params.momentum()(2);
  const auto trackAlpha = -trackPPhi / trackPR;
  const auto trackBeta  = -trackPZ / trackPR;

  if(Verbosity() > 3)
    std::cout << "Track angles " << trackPPhi << ", " << trackPR
	      << ", " << trackPZ << ", " << trackAlpha << ", " << trackBeta
	      << std::endl;

  if(std::abs(trackAlpha) > m_maxTAlpha
     or std::abs(drphi) > m_maxResidualDrphi)
    return;

  if(std::abs(trackBeta) > m_maxTBeta
     or std::abs(dz) > m_maxResidualDz)
    return;
    
  tanBeta = trackBeta;
  tanAlpha = trackAlpha;
  
  Acts::Vector3D globClus(cluster->getX(), 
			  cluster->getY(), 
			  cluster->getZ());
  const auto index = getCell(globClus);
  
  if(Verbosity() > 3)
    std::cout << "Bin index found is " << index << std::endl;

  if(index < 0 ) return;

  if( m_savehistograms )
  {
    h_index->Fill(index);
    h_alpha->Fill(tanAlpha, drphi);
    h_beta->Fill(tanBeta, dz);
    h_rphiResid->Fill(clusR , drphi);
    h_zResid->Fill(stateZ , dz);
    h_etaResid->Fill(trackEta, clusEta - trackEta);
    h_zResidLayer->Fill(clusR , dz);
    h_etaResidLayer->Fill(clusR , clusEta - trackEta);
    
    residTup->Fill();
  }
  
  // Fill distortion matrices
  m_matrix_container->add_to_lhs(index, 0, 0, 1./erp );
  m_matrix_container->add_to_lhs(index, 0, 1, 0 );
  m_matrix_container->add_to_lhs(index, 0, 2, trackAlpha/erp );
  
  m_matrix_container->add_to_lhs(index, 1, 0, 0 );
  m_matrix_container->add_to_lhs(index, 1, 1, 1./ez );
  m_matrix_container->add_to_lhs(index, 1, 2, trackBeta/ez );
  
  m_matrix_container->add_to_lhs(index, 2, 0, trackAlpha/erp );
  m_matrix_container->add_to_lhs(index, 2, 1, trackBeta/ez );
  m_matrix_container->add_to_lhs(index, 2, 2, square(trackAlpha)/erp + square(trackBeta)/ez );
  
  m_matrix_container->add_to_rhs(index, 0, drphi/erp );
  m_matrix_container->add_to_rhs(index, 1, dz/ez );
  m_matrix_container->add_to_rhs(index, 2, trackAlpha*drphi/erp + trackBeta*dz/ez );
  
  // update entries in cell
  m_matrix_container->add_to_entries(index);
  
  return;
}

int PHTpcResiduals::getCell(const Acts::Vector3D& loc)
{

  // get grid dimensions from matrix container
  int phibins = 0;
  int rbins = 0;
  int zbins = 0;
  m_matrix_container->get_grid_dimensions( phibins, rbins, zbins );
  
  // phi
  float phi = std::atan2(loc(1), loc(0));
  while( phi < m_phiMin ) phi += 2.*M_PI;
  while( phi >= m_phiMax ) phi -= 2.*M_PI;
  const int iphi = phibins * (phi - m_phiMin) / (m_phiMax - m_phiMin);
  
  // r
  const auto r = get_r( loc(0), loc(1) );
  if( r < m_rMin || r >= m_rMax ) return -1;
  const int ir = rbins * (r - m_rMin) / (m_rMax - m_rMin);
  
  // z
  const auto z = loc(2);
  if( z < m_zMin || z >= m_zMax ) return -1;
  const int iz = zbins * (z - m_zMin) / (m_zMax - m_zMin);

  // get index from matrix container
  return m_matrix_container->get_cell_index( iphi, ir, iz );

}

int PHTpcResiduals::createNodes(PHCompositeNode */*topNode*/)
{

  return Fun4AllReturnCodes::EVENT_OK;
}

int PHTpcResiduals::getNodes(PHCompositeNode *topNode)
{
  m_vertexMap = findNode::getClass<SvtxVertexMap>(topNode,
						  "SvtxVertexMap");
  if(!m_vertexMap)
    {
      std::cout << PHWHERE << "No vertex map, exiting." 
		<< std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }

  m_surfMaps = findNode::getClass<ActsSurfaceMaps>(topNode,
						   "ActsSurfaceMaps");
  if(!m_surfMaps)
    {
      std::cout << PHWHERE << "No Acts surface maps, exiting"
		<< std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }
 
  m_clusterContainer = findNode::getClass<TrkrClusterContainer>(topNode,
								"TRKR_CLUSTER");
  if(!m_clusterContainer)
    {
      std::cout << PHWHERE << "No TRKR_CLUSTER node on node tree. Exiting."
		<< std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }

  m_tGeometry = findNode::getClass<ActsTrackingGeometry>(topNode, "ActsTrackingGeometry");
  if(!m_tGeometry)
    {
      std::cout << "ActsTrackingGeometry not on node tree. Exiting."
		<< std::endl;
      
      return Fun4AllReturnCodes::ABORTEVENT;
    }

  m_trackMap = findNode::getClass<SvtxTrackMap>(topNode, "SvtxSiliconMMTrackMap");
  
  if (!m_trackMap)
    {
      std::cout << PHWHERE << "SvtxSiliconMMTrackMap not on node tree. Exiting."
		<< std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }
  

  return Fun4AllReturnCodes::EVENT_OK;
}

void PHTpcResiduals::makeHistograms()
{  
  
  m_histogramfile.reset( new TFile(m_histogramfilename.c_str(), "RECREATE") );
  m_histogramfile->cd();

  const auto total_bins = m_matrix_container->get_grid_size();  
  h_beta = new TH2F("betadz",";tan#beta; #Deltaz [cm]",100,-0.5,0.5,100,-0.5,0.5);
  h_alpha = new TH2F("alphardphi",";tan#alpha; r#Delta#phi [cm]", 100,-0.5,0.5,100,-0.5,0.5);
  h_index = new TH1F("index",";index",total_bins, 0, total_bins);
  h_rphiResid = new TH2F("rphiResid", ";r [cm]; #Deltar#phi [cm]",
			 60, 20, 80, 500, -2, 2);
  h_zResid = new TH2F("zResid", ";z [cm]; #Deltaz [cm]",
		      200, -100, 100, 1000, -2, 2);
  h_etaResid = new TH2F("etaResid", ";#eta;#Delta#eta",
			20, -1, 1, 500, -0.2, 0.2);
  h_etaResidLayer = new TH2F("etaResidLayer", ";r [cm]; #Delta#eta",
			     60, 20, 80, 500, -0.2, 0.2);
  h_zResidLayer = new TH2F("zResidLayer", ";r [cm]; #Deltaz [cm]",
			   60, 20, 80, 1000, -2, 2);
  residTup = new TTree("residTree","tpc residual info");
  residTup->Branch("tanAlpha",&tanAlpha,"tanAlpha/D");
  residTup->Branch("tanBeta",&tanBeta,"tanBeta/D");
  residTup->Branch("drphi",&drphi,"drphi/D");
  residTup->Branch("dz",&dz,"dz/D");
//   residTup->Branch("cell",&cell,"cell/I");
  residTup->Branch("clusR",&clusR,"clusR/D");
  residTup->Branch("clusPhi",&clusPhi,"clusPhi/D");
  residTup->Branch("clusZ",&clusZ,"clusZ/D");
  residTup->Branch("statePhi",&statePhi,"statePhi/D");
  residTup->Branch("stateZ",&stateZ,"stateZ/D");
  residTup->Branch("stateR",&stateR,"stateR/D");
  residTup->Branch("stateRPhiErr",&stateRPhiErr,"stateRPhiErr/D");
  residTup->Branch("stateZErr",&stateZErr,"stateZErr/D");
  residTup->Branch("clusRPhiErr",&clusRPhiErr,"clusRPhiErr/D");
  residTup->Branch("clusZErr",&clusZErr,"clusZErr/D");
//   residTup->Branch("ir",&ir,"ir/I");
//   residTup->Branch("iz",&iz,"iz/I");
//   residTup->Branch("iphi",&iphi,"iphi/I");
  residTup->Branch("cluskey",&cluskey,"cluskey/I");
  residTup->Branch("event",&m_event,"event/I");

}

void PHTpcResiduals::setGridDimensions(const int phiBins, const int rBins, const int zBins)
{ m_matrix_container->set_grid_dimensions( phiBins, rBins, zBins ); }
