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

#include <Acts/Geometry/GeometryIdentifier.hpp>
#include <Acts/MagneticField/ConstantBField.hpp>
#include <Acts/MagneticField/InterpolatedBFieldMap.hpp>
#include <Acts/MagneticField/SharedBField.hpp>
#include <Acts/Propagator/EigenStepper.hpp>
#include <Acts/Surfaces/Surface.hpp>

#include <ActsExamples/Plugins/BField/ScalableBField.hpp>

#include <cmath>
#include <TFile.h>
#include <TH3.h>

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
  , m_actsProtoTracks(nullptr)
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

  for(auto &[trackKey, track] : *m_actsProtoTracks)
    {
      if(checkTrack(track))
	processTrack(track);
      
    }
  
  return Fun4AllReturnCodes::EVENT_OK;
}

bool PHTpcResiduals::checkTrack(ActsTrack& track)
{
 
  if(track.getTrackParams().transverseMomentum() < 0.5)
    return false;

  int nMvtxHits = 0;
  int nInttHits = 0;
  int nMMHits   = 0;
 
  for(const auto sl : track.getSourceLinks())
    {
      const auto vol = sl.referenceSurface().geometryId().volume();
      if(vol == 10) 
	nMvtxHits++;
      if(vol == 12) 
	nInttHits++;
      if(vol == 16) 
	nMMHits++;
    }

  if(Verbosity() > 2)
    std::cout << "Number of mvtx/intt/MM hits "
	      << nMvtxHits << "/" << nInttHits << "/" 
	      << nMMHits << std::endl;

  if(nMvtxHits < 2 || nInttHits < 2 || nMMHits < 2)
    return false;

  return true;

}

void PHTpcResiduals::processTrack(ActsTrack& track)
{

  const auto trackParams = track.getTrackParams();
  
  int initNBadProps = m_nBadProps;
  for(const auto sl : track.getSourceLinks())
    {
      /// Only analyze TPC 
      if(sl.referenceSurface().geometryId().volume() != 14)
	continue;
      
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
	  

	  calculateTpcResiduals(trackStateParams, sl);
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
			   const ActsExamples::TrackParameters& params,
			   const SourceLink& sl)
{
  
  if(Verbosity() > 1)
    std::cout << "Propagating silicon+MM fit params momentum: " 
	      << params.momentum() << " and position " 
	      << params.position(m_tGeometry->geoContext)
	      << std::endl;

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
      if(Verbosity() > 3)
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
			  const SourceLink& sl)
{
  Acts::Vector2D local(sl.location()(0), sl.location()(1));
  
  auto globalSL = sl.referenceSurface().localToGlobal(
				         m_tGeometry->geoContext,
					 local,
					 Acts::Vector3D(1,1,1));

  /// Get all the relevant information for residual calculation
  const auto clusR = sqrt(pow(globalSL.x(), 2) +
			  pow(globalSL.y(), 2))
    / Acts::UnitConstants::cm;
  const auto clusPhi = std::atan2(globalSL.y(), globalSL.x());
  const auto clusZ = globalSL.z() / Acts::UnitConstants::cm;
  const auto clusRPhiErr = sqrt(sl.covariance()(Acts::eBoundLoc0,
						Acts::eBoundLoc0))
    / Acts::UnitConstants::cm;
  const auto clusZErr = sqrt(sl.covariance()(Acts::eBoundLoc1,
					     Acts::eBoundLoc1))
    / Acts::UnitConstants::cm;
  
  if(Verbosity() > 3)
    std::cout <<"Cluster phi and z " << clusPhi << "+/-" << clusZErr
	      <<" and " << clusZ << "+/-" << clusZErr << std::endl;


  if(clusRPhiErr < 0.015)
    continue;
  if(clusZErr < 0.015)
    continue;

  const auto globalStatePos = params.position(m_tGeometry->geoContext);
  const auto globalStateCov = *params.covariance();
  const auto stateRPhiErr = sqrt(globalStateCov(Acts::eBoundLoc0,
						Acts::eBoundLoc0))
    / Acts::UnitConstants::cm;
  const auto stateZErr = sqrt(globalStateCov(Acts::eBoundLoc1,
					     Acts::eBoundLoc1))
    / Acts::UnitConstants::cm;
 
  /// We don't have to extrapolate the track parameters to the cluster
  /// r because the Acts::Propagator already propagated the parameters
  /// to the surface where the cluster exists (e.g. the same r)
  const auto statePhi = std::atan2(globalStatePos.y(),
				   globalStatePos.x());
  const auto stateZ = globalStatePos.z() / Acts::UnitConstants::cm;
  
  if(Verbosity() > 3)
    std::cout << "State phi and z " << statePhi << "+/-" << stateRPhiErr
	      << " and " << stateZ << "+/-" << stateZErr << std::endl;

  const auto erp = pow(stateRPhiErr, 2) + pow(clusRPhiErr, 2);
  const auto ez = pow(stateZErr, 2) + pow(clusZErr, 2);

  const auto dPhi = clusPhi - statePhi;

  /// Calculate residuals
  const auto drphi = clusR * deltaPhi(dPhi);
  const auto dz  = clusZ - stateZ;

  const auto trackEta 
    = std::atanh(params.momentum().z() / params.absoluteMomentum());
  const auto clusEta = std::atanh(clusZ / (globalSL.norm() / Acts::UnitConstants::cm));

  h_rphiResid->Fill(clusR , drphi);
  h_zResid->Fill(stateZ , dz);
  h_etaResid->Fill(trackEta, clusEta - trackEta);
  h_zResidLayer->Fill(clusR , dz);
  h_etaResidLayer->Fill(clusR , clusEta - trackEta);

  const auto trackPPhi = -params.momentum()(0) * std::sin(statePhi) +
    params.momentum()(1) * std::cos(statePhi);
  const auto trackPR = params.momentum()(0) * std::cos(statePhi) +
    params.momentum()(1) * std::sin(statePhi);
  
  const auto trackPZ    = params.momentum()(1);
  const auto trackAlpha = -trackPPhi / trackPR;
  const auto trackBeta  = -trackPZ / trackPR;

  if(std::abs(trackAlpha) > m_maxTAlpha
     || std::abs(drphi) > m_maxResidualDrphi)
    return;

  if(std::abs(trackBeta) > m_maxTBeta
     || std::abs(dz) > m_maxResidualDz)
    return;

  const auto index = getCell(sl.referenceSurface().geometryId().layer(),
			     globalSL);

  if(index < 0 || index > m_totalBins)
    return;

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
	  if(Verbosity() > -1)
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

	if(Verbosity() > -1)
	  std::cout << "Bin setting for index " << cell << " with counts "
		    << m_clusterCount.at(cell) << " has settings : "
		    <<"  "<<result(0) <<"+/-" << std::sqrt(cov(0,0))
		    <<"  "<<result(1) <<"+/-" << std::sqrt(cov(1,1))
		    <<"  "<<result(2) <<"+/-" << std::sqrt(cov(2,2))
		    <<std::endl;

      }
    }
  }
  
  /// Create output TH3s
  TFile *outputFile = new TFile(m_outputfile.c_str(), "RECREATE");
  outputFile->cd();

  h_rphiResid->Write();
  h_etaResid->Write();
  h_zResidLayer->Write();
  h_etaResidLayer->Write();
  h_zResid->Write();
    

  for(const auto& h : {hentries, hphi, hr, hz})
    h->Write();
  outputFile->Close();

}

int PHTpcResiduals::getCell(const int actsLayer, 
			    const Acts::Vector3D& loc)
{
  /// Divide by two because ACTS definition is twice ours, subtract
  /// 1 to get layer number from 0-47 instead of 1-48
  const auto layer = (actsLayer / 2.) -1;
  const int ir = m_rBins * layer / m_nLayersTpc;

  auto clusPhi = deltaPhi(std::atan2(loc(1), loc(0)));
  const int iphi = m_phiBins * (clusPhi * M_PI) / (2. * M_PI);

  const int iz = m_zBins * (loc(2) - m_zMin) / (m_zMax - m_zMin);
 
  return getCell(iz, ir, iphi);
}

int PHTpcResiduals::getCell(const int iz, const int ir, 
			    const int iphi)
{
  if( ir < 0 || ir >= m_rBins ) return -1;
  if( iphi < 0 || iphi >= m_phiBins ) return -1;
  if( iz < 0 || iz >= m_zBins ) return -1;
  return iz + m_zBins*( ir + m_rBins*iphi );
}

int PHTpcResiduals::createNodes(PHCompositeNode *topNode)
{

  return Fun4AllReturnCodes::EVENT_OK;
}

int PHTpcResiduals::getNodes(PHCompositeNode *topNode)
{
  m_tGeometry = findNode::getClass<ActsTrackingGeometry>(topNode, "ActsTrackingGeometry");
  if(!m_tGeometry)
    {
      std::cout << "ActsTrackingGeometry not on node tree. Exiting."
		<< std::endl;
      
      return Fun4AllReturnCodes::ABORTEVENT;
    }

  m_actsProtoTracks = findNode::getClass<std::map<unsigned int, ActsTrack>>(topNode, "ActsTrackMap");
  
  if (!m_actsProtoTracks)
    {
      std::cout << "Acts proto tracks not on node tree. Exiting."
		<< std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }
  

  return Fun4AllReturnCodes::EVENT_OK;
}

void PHTpcResiduals::makeHistograms()
{
 
  h_rphiResid = new TH2F("rphiResid", ";r [cm]; #Deltar#phi [cm]",
			 60, 20, 80, 50, -2, 2);
  h_zResid = new TH2F("zResid", ";z [cm]; #Deltaz [cm]",
		      200, -100, 100, 100, -2, 2);
  h_etaResid = new TH2F("etaResid", ";#eta;#Delta#eta",
			20, -1, 1, 50, -0.2, 0.2);
  h_etaResidLayer = new TH2F("etaResidLayer", ";r [cm]; #Delta#eta",
			     60, 20, 80, 50, -0.2, 0.2);
  h_zResidLayer = new TH2F("zResidLayer", ";r [cm]; #Deltaz [cm]",
			   60, 20, 80, 100, -2, 2);

}

void PHTpcResiduals::setGridDimensions(const int phiBins, const int rBins,
				       const int zBins)
{
  m_zBins = zBins;
  m_phiBins = phiBins;
  m_rBins = rBins;
  m_totalBins = m_zBins * m_phiBins * m_rBins;

}
