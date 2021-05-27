#include "PHTpcResiduals.h"

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
{}

PHTpcResiduals::~PHTpcResiduals()
{}

int PHTpcResiduals::Init(PHCompositeNode *topNode)
{
  m_rhs = std::vector<Acts::Vector3D>(m_totalBins,
  				      Acts::Vector3D::Zero());
  m_lhs = std::vector<Acts::SymMatrix3D>(m_totalBins,
					 Acts::SymMatrix3D::Zero());
  m_clusterCount = std::vector<int>(m_totalBins, 0);

  makeHistograms();

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

int PHTpcResiduals::End(PHCompositeNode *topNode)
{
  if(Verbosity() > 0)
    std::cout << "Number of bad SL propagations " 
	      << m_nBadProps << std::endl;

  calculateDistortions(topNode);

  return Fun4AllReturnCodes::EVENT_OK;
}


int PHTpcResiduals::processTracks(PHCompositeNode *topNode)
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

  /// Require at least 2 hits in each detector
  if(nMvtxHits < 2 or nInttHits < 2 or nMMHits < 2)
    return false;

  return true;

}

SourceLink PHTpcResiduals::makeSourceLink(TrkrCluster* cluster)
{
  auto key = cluster->getClusKey();
  auto subsurfkey = cluster->getSubSurfKey();
      
  /// Make a safety check for clusters that couldn't be attached
  /// to a surface
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
  
  /// If it can't be found, return nullptr
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
  
  /// If it can't be found, return nullptr to skip this cluster
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
      
      /// only propagate to tpc surfaces
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
  /// Get all the relevant information for residual calculation
  clusR = sqrt(pow(cluster->getX(), 2) +
	       pow(cluster->getY(), 2));
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

  stateR = sqrt(pow(globStateX, 2) +
		pow(globStateY, 2) );
  
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

  const auto erp = pow(clusRPhiErr, 2);
  const auto ez = pow(clusZErr, 2);

  const auto dPhi = clusPhi - statePhi;

  /// Calculate residuals
  drphi = clusR * deltaPhi(dPhi);
  dz  = clusZ - stateZ;

  if(Verbosity() > 3)
    std::cout << "TPC residuals " << drphi << "   " << dz << std::endl;
  
  const auto trackEta 
    = std::atanh(params.momentum().z() / params.absoluteMomentum());
  const auto clusEta = std::atanh(clusZ / sqrt(pow(cluster->getX(), 2) +
					       pow(cluster->getY(), 2) +
					       pow(cluster->getZ(), 2)));

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
  
  ir = m_rBins * (clusR - m_rMin) / (m_rMax - m_rMin);
  iphi = m_phiBins * (clusPhi - m_phiMin) / (m_phiMax - m_phiMin);
  iz = m_zBins * (clusZ - m_zMin) / (m_zMax - m_zMin);
  
  tanBeta = trackBeta;
  tanAlpha = trackAlpha;
  
  Acts::Vector3D globClus(cluster->getX(), 
			  cluster->getY(), 
			  cluster->getZ());
  const auto index = getCell(globClus);
  cell = index;
  
  if(Verbosity() > 3)
    std::cout << "Bin index found is " << index << std::endl;

  if(index < 0 || index > m_totalBins)
    return;

  if(index == 0)
    std::cout << "Values are : " <<iphi<<","<<ir<<","<<iz<<std::endl;

  h_index->Fill(cell);
  h_alpha->Fill(tanAlpha, drphi);
  h_beta->Fill(tanBeta, dz);
  h_rphiResid->Fill(clusR , drphi);
  h_zResid->Fill(stateZ , dz);
  h_etaResid->Fill(trackEta, clusEta - trackEta);
  h_zResidLayer->Fill(clusR , dz);
  h_etaResidLayer->Fill(clusR , clusEta - trackEta);

  residTup->Fill();

  /// Fill distortion matrices
  m_lhs[index](0,0) += 1. / erp;
  m_lhs[index](0,1) += 0;
  m_lhs[index](0,2) += trackAlpha / erp;
  
  m_lhs[index](1,0) += 0;
  m_lhs[index](1,1) += 1. / ez;
  m_lhs[index](1,2) += trackBeta / ez;
  
  m_lhs[index](2,0) += trackAlpha / erp;
  m_lhs[index](2,1) += trackBeta / ez;
  m_lhs[index](2,2) += pow(trackAlpha, 2) / erp 
                     + pow(trackBeta, 2) / ez;
  
  m_rhs[index](0,0) += drphi / erp;
  m_rhs[index](1,0) += dz / ez;
  m_rhs[index](2,0) += trackAlpha * drphi / erp + trackBeta * dz / ez;

  m_clusterCount[index]++;

  return;
}

void PHTpcResiduals::calculateDistortions(PHCompositeNode *topNode)
{
  auto hentries = new TH3F( "hentries_rec", "hentries_rec", m_phiBins, 
			    m_phiMin, m_phiMax, m_rBins, m_rMin, m_rMax, 
			    m_zBins, m_zMin, m_zMax );
  auto hphi = new TH3F( "hDistortionP_rec", "hDistortionP_rec", m_phiBins, 
			m_phiMin, m_phiMax, m_rBins, m_rMin, m_rMax, 
			m_zBins, m_zMin, m_zMax );
  auto hz = new TH3F( "hDistortionZ_rec", "hDistortionZ_rec", m_phiBins, 
		      m_phiMin, m_phiMax, m_rBins, m_rMin, m_rMax, m_zBins,
		      m_zMin, m_zMax );
  auto hr = new TH3F( "hDistortionR_rec", "hDistortionR_rec", m_phiBins, 
		      m_phiMin, m_phiMax, m_rBins, m_rMin, m_rMax, m_zBins,
		      m_zMin, m_zMax );

  for( auto h : { hentries, hphi, hz, hr } )
    {
      h->GetXaxis()->SetTitle( "#phi [rad]" );
      h->GetYaxis()->SetTitle( "r [cm]" );
      h->GetZaxis()->SetTitle( "z [cm]" );
    }
  
  for(int iphi = 0; iphi < m_phiBins; ++iphi) {
    for(int ir = 0; ir < m_rBins; ++ir) {
      for(int iz = 0; iz < m_zBins; ++iz) {

	const auto cell = getCell(iz, ir, iphi);	

	if(m_clusterCount.at(cell) < m_minClusCount) {
	  if(Verbosity() > 10)
	    std::cout << "Num clusters in bin " << cell 
		      << " is " << m_clusterCount.at(cell) 
		      << std::endl;
	  continue;
	}
	  
	const auto cov = m_lhs.at(cell).inverse();
	auto partialLu = m_lhs.at(cell).partialPivLu();
	const auto result = partialLu.solve(m_rhs.at(cell));

	// fill histograms
	hentries->SetBinContent( iphi+1, ir+1, iz+1, 
				 m_clusterCount.at(cell) );
	
	hphi->SetBinContent( iphi+1, ir+1, iz+1, result(0) );
	hphi->SetBinError( iphi+1, ir+1, iz+1, std::sqrt( cov(0,0) ) );
	
	hz->SetBinContent( iphi+1, ir+1, iz+1, result(1) );
	hz->SetBinError( iphi+1, ir+1, iz+1, std::sqrt( cov(1,1) ) );
	
	hr->SetBinContent( iphi+1, ir+1, iz+1, result(2) );
	hr->SetBinError( iphi+1, ir+1, iz+1, std::sqrt( cov(2,2) ) );

	if(Verbosity() > 10)
	  std::cout << "Bin setting for index " << cell << " with counts "
		    << m_clusterCount.at(cell) << " has settings : "
		    << " drphi:  "<<result(0) << "+/-" << std::sqrt(cov(0,0))
		    << " dz: "<<result(1) << "+/-" << std::sqrt(cov(1,1))
		    << " dr: "<<result(2) << "+/-" << std::sqrt(cov(2,2))
		    << std::endl;

      }
    }
  }
  
  /// Create output TH3s
    

  m_outputFile->cd();
  residTup->Write();

  h_rphiResid->Write();
  h_etaResid->Write();
  h_zResidLayer->Write();
  h_etaResidLayer->Write();
  h_zResid->Write();
  h_index->Write();
  h_alpha->Write();
  h_beta->Write();

  for(const auto& h : {hentries, hphi, hr, hz})
    h->Write();


  createHistogram(hentries, "hentries")->Write();
  createHistogram(hphi, "hIntDistortionP")->Write();
  createHistogram(hr, "hIntDistortionR")->Write();
  createHistogram(hz, "hIntDistortionZ")->Write();

  m_outputFile->Close();

}

int PHTpcResiduals::getCell(const Acts::Vector3D& loc)
{
  const float r = sqrt(loc(0) * loc(0) + loc(1) * loc(1));
  const auto clusPhi = deltaPhi(std::atan2(loc(1), loc(0)));
  const float z = loc(2);
  
  const int ir = m_rBins * (r - m_rMin) / (m_rMax - m_rMin);
  const int iphi = m_phiBins * (clusPhi - m_phiMin) / (m_phiMax - m_phiMin);
  const int iz = m_zBins * (z - m_zMin) / (m_zMax - m_zMin);
  
  return getCell(iz, ir, iphi);
}

int PHTpcResiduals::getCell(const int iz, const int ir, 
			    const int iphi)
{
  if( ir < 0 || ir >= m_rBins ) return -1;
  if( iphi < 0 || iphi >= m_phiBins ) return -1;
  if( iz < 0 || iz >= m_zBins ) return -1;

  return iz + m_zBins * ( ir + m_rBins * iphi );
}

int PHTpcResiduals::createNodes(PHCompositeNode *topNode)
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
 
  m_outputFile = new TFile(m_outputfile.c_str(), "RECREATE");
  m_outputFile->cd();

  h_beta = new TH2F("betadz",";tan#beta; #Deltaz [cm]",100,-0.5,0.5,100,-0.5,0.5);
  h_alpha = new TH2F("alphardphi",";tan#alpha; r#Delta#phi [cm]", 100,-0.5,0.5,100,-0.5,0.5);
  h_index = new TH1F("index",";index",m_totalBins, 0, m_totalBins);
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
  residTup->Branch("cell",&cell,"cell/I");
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
  residTup->Branch("ir",&ir,"ir/I");
  residTup->Branch("iz",&iz,"iz/I");
  residTup->Branch("iphi",&iphi,"iphi/I");
  residTup->Branch("cluskey",&cluskey,"cluskey/I");
  residTup->Branch("event",&m_event,"event/I");

}

void PHTpcResiduals::setGridDimensions(const int phiBins, const int rBins,
				       const int zBins)
{
  m_zBins = zBins;
  m_phiBins = phiBins;
  m_rBins = rBins;
  m_totalBins = m_zBins * m_phiBins * m_rBins;
}


TH3* PHTpcResiduals::createHistogram(TH3* hin, const TString& name)
{

  std::array<int, 3> bins;
  std::array<double, 3> x_min;
  std::array<double, 3> x_max;
  
  int index = 0;
  for( const auto axis:{ hin->GetXaxis(), hin->GetYaxis(), hin->GetZaxis() } )
    {
      const auto bin_width = (axis->GetXmax() - axis->GetXmin())/axis->GetNbins();
      
      // increase the number of bins by two
      bins[index] = axis->GetNbins()+2;
      
      // update axis limits accordingly
      x_min[index] = axis->GetXmin()-bin_width;
      x_max[index] = axis->GetXmax()+bin_width;
      ++index;
    }

  // create new histogram
  auto hout = new TH3F( name, name,
			bins[0], x_min[0], x_max[0],
			bins[1], x_min[1], x_max[1],
			bins[2], x_min[2], x_max[2] );
  
  // update axis legend
  hout->GetXaxis()->SetTitle( hin->GetXaxis()->GetTitle() );
  hout->GetYaxis()->SetTitle( hin->GetYaxis()->GetTitle() );
  hout->GetZaxis()->SetTitle( hin->GetZaxis()->GetTitle() );
  
  // copy content
  const auto phibins = hin->GetXaxis()->GetNbins();
  const auto rbins = hin->GetYaxis()->GetNbins();
  const auto zbins = hin->GetZaxis()->GetNbins();
  
  // fill center
  for( int iphi = 0; iphi < phibins; ++iphi )
    for( int ir = 0; ir < rbins; ++ir )
      for( int iz = 0; iz < zbins; ++iz )
	{
	  hout->SetBinContent( iphi+2, ir+2, iz+2, hin->GetBinContent( iphi+1, ir+1, iz+1 ) );
	  hout->SetBinError( iphi+2, ir+2, iz+2, hin->GetBinError( iphi+1, ir+1, iz+1 ) );
	}
  
  // fill guarding phi bins
  for( int ir = 0; ir < rbins+2; ++ir )
    for( int iz = 0; iz < zbins+2; ++iz )
      {
	hout->SetBinContent( 1, ir+1, iz+1, hout->GetBinContent( 2, ir+1, iz+1 ) );
	hout->SetBinError( 1, ir+1, iz+1, hout->GetBinError( 2, ir+1, iz+1 ) );
	
	hout->SetBinContent( phibins+2, ir+1, iz+1, hout->GetBinContent( phibins+1, ir+1, iz+1 ) );
	hout->SetBinError( phibins+2, ir+1, iz+1, hout->GetBinError( phibins+1, ir+1, iz+1 ) );
      }
  
  // fill guarding r bins
  for( int iphi = 0; iphi < phibins+2; ++iphi )
    for( int iz = 0; iz < zbins+2; ++iz )
      {
	hout->SetBinContent( iphi+1, 1, iz+1, hout->GetBinContent( iphi+1, 2, iz+1 ) );
	hout->SetBinError( iphi+1, 1, iz+1, hout->GetBinError( iphi+1, 2, iz+1 ) );
	
	hout->SetBinContent( iphi+1, rbins+2, iz+1, hout->GetBinContent( iphi+1, rbins+1, iz+1 ) );
	hout->SetBinError( iphi+1, rbins+1, iz+1, hout->GetBinError( iphi+1, rbins+1, iz+1 ) );
      }
  
  // fill guarding z bins
  for( int iphi = 0; iphi < phibins+2; ++iphi )
    for( int ir = 0; ir < rbins+2; ++ir )
      {
	hout->SetBinContent( iphi+1, ir+1, 1, hout->GetBinContent( iphi+1, ir+1, 2 ) );
	hout->SetBinError( iphi+1, ir+1, 1, hout->GetBinError( iphi+1, ir+1, 2 ) );
	
	hout->SetBinContent( iphi+1, ir+1, zbins+2, hout->GetBinContent( iphi+1, ir+1, zbins+1 ) );
	hout->SetBinError( iphi+1, ir+1, zbins+2, hout->GetBinError( iphi+1, ir+1, zbins+1 ) );
      }
  
  return hout;
  
}





