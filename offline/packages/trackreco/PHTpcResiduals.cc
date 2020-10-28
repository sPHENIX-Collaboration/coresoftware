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

#include <g4detectors/PHG4CylinderCellGeom.h>
#include <g4detectors/PHG4CylinderCellGeomContainer.h>

#include <Acts/Geometry/GeometryIdentifier.hpp>
#include <Acts/MagneticField/ConstantBField.hpp>
#include <Acts/MagneticField/InterpolatedBFieldMap.hpp>
#include <Acts/MagneticField/SharedBField.hpp>
#include <Acts/Propagator/EigenStepper.hpp>
#include <Acts/Surfaces/Surface.hpp>

#include <ActsExamples/Plugins/BField/ScalableBField.hpp>

#include <cmath>
#include <TGraphErrors.h>

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
{
}


PHTpcResiduals::~PHTpcResiduals()
{
}

int PHTpcResiduals::Init(PHCompositeNode *topNode)
{
  m_rhs = std::vector<Acts::Vector3D>(m_totalBins,
  				      Acts::Vector3D::Zero());
  m_lhs = std::vector<Acts::SymMatrix3D>(m_totalBins,
					 Acts::SymMatrix3D::Zero());
  m_clusterCount = std::vector<int>(m_totalBins, 0);

  if(m_outputRoot)
    outfile = new TFile(std::string(Name() + ".root").c_str(), 
			  "recreate");
  
  h_rphiResid = new TH2F("rphiResid",";r [cm]; #Deltar#phi [mm]",
			 60,20,80,50,-10,10);
  h_zResid = new TH2F("zResid",";z [cm]; #Deltaz [mm]",
		      200,-100,100,100,-10,10);
  h_etaResid = new TH2F("etaResid",";#eta;#Delta#eta",
			20,-1,1,50,-0.2,0.2);
  h_etaResidLayer = new TH2F("etaResidLayer",";r [cm]; #Delta#eta",
			     60,20,80,50,-0.2,0.2);
  h_zResidLayer = new TH2F("zResidLayer",";r [cm]; #Deltaz [mm]",
			   60,20,80,100,-10,10);
  return Fun4AllReturnCodes::EVENT_OK;

}

int PHTpcResiduals::InitRun(PHCompositeNode *topNode)
{
  if(getNodes(topNode) != Fun4AllReturnCodes::EVENT_OK)
    return Fun4AllReturnCodes::ABORTEVENT;

  return Fun4AllReturnCodes::EVENT_OK;
}

int PHTpcResiduals::ResetEvent(PHCompositeNode *topNode)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

int PHTpcResiduals::process_event(PHCompositeNode *topNode)
{
  if(Verbosity() > 0)
    std::cout <<"Starting PHTpcResiduals event " 
	      << m_event << std::endl;

  int returnVal = processTracks(topNode);

  if(Verbosity() > 0)
    std::cout <<"Finished PHTpcResiduals event " 
	      << m_event << std::endl;
  
  m_event++;

  return returnVal;
}

int PHTpcResiduals::End(PHCompositeNode *topNode)
{
  if(Verbosity() > 0)
    std::cout << "Number of bad propagations " 
	      << m_nBadProps << std::endl;

  calculateDistortions(topNode);

  if(m_outputRoot)
    {
      outfile->cd();
      h_rphiResid->Write();
      h_etaResid->Write();
      h_zResidLayer->Write();
      h_etaResidLayer->Write();
      h_zResid->Write();
      outfile->Write();
      outfile->Close();
    }
  return Fun4AllReturnCodes::EVENT_OK;
}


int PHTpcResiduals::processTracks(PHCompositeNode *topNode)
{


  std::map<unsigned int, ActsTrack>::iterator trackIter;
  for(trackIter = m_actsProtoTracks->begin();
      trackIter != m_actsProtoTracks->end();
      ++trackIter)
    {
      auto track = trackIter->second;
      processTrack(track);
      
    }
  
  return Fun4AllReturnCodes::EVENT_OK;
}

void PHTpcResiduals::processTrack(ActsTrack& track)
{
 
  auto sourceLinks = track.getSourceLinks();

  /// Check that there are silicon+MM sourcelinks in this
  /// particular track
  bool MM = false;
  bool silicon = false;
  for(auto sl : sourceLinks)
    {
      auto volume = sl.referenceSurface().geometryId().volume();
      if(volume == 10 or volume == 12)
	silicon = true;
      if(volume == 16)
	MM = true;
    }

  if(!MM or !silicon)
    return;

  const auto trackParams = track.getTrackParams();
  
  int initNBadProps = m_nBadProps;
  for(auto sl : sourceLinks)
    {
      /// Only analyze TPC 
      if(sl.referenceSurface().geometryId().volume() != 14)
	continue;
      
      //std::cout << "Propagating track state" << std::endl;
      auto result = propagateTrackState(trackParams, sl);
      //std::cout << "finished propagating"<<std::endl;
      if(result.ok())
	{	  
	  auto trackStateParams = std::move(**result);
	  if(Verbosity() > 1)
	    {
	      std::cout << "Silicon+MM momentum : " 
			<< trackParams.momentum()
			<< std::endl
			<< "Propagator momentum : " 
			<< trackStateParams.momentum()
			<< std::endl;
	    }

	  calculateTpcResiduals(trackStateParams, sl);
	}
      else
	{
	  m_nBadProps++;
	
	  continue;
	}
    } 

  if(m_nBadProps > initNBadProps && Verbosity() > 0)
    {
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
			  pow(globalSL.y(), 2));
  const auto clusPhi = std::atan2(globalSL.y(), globalSL.x());
  const auto clusZ = globalSL.z();
  const auto clusRPhiErr = sqrt(sl.covariance()(Acts::eBoundLoc0,
						Acts::eBoundLoc0));
  const auto clusZErr = sqrt(sl.covariance()(Acts::eBoundLoc1,
					     Acts::eBoundLoc1));
  

  const auto globalStatePos = params.position(m_tGeometry->geoContext);
  const auto globalStateCov = *params.covariance();
  const auto stateRPhiErr = sqrt(globalStateCov(Acts::eBoundLoc0,
						Acts::eBoundLoc0));
  const auto stateZErr = sqrt(globalStateCov(Acts::eBoundLoc1,
					     Acts::eBoundLoc1));
 
  const auto statePhi = std::atan2(globalStatePos.y(),
				   globalStatePos.x());
  const auto stateZ = globalStatePos.z();
  
  const auto erp = pow(stateRPhiErr, 2) + pow(clusRPhiErr, 2);
  const auto ez = pow(stateZErr, 2) + pow(clusZErr, 2);

  const auto dPhi = clusPhi - statePhi;

  /// Calculate residuals
  const auto drphi = clusR * deltaPhi(dPhi);
  const auto dz  = clusZ - stateZ;

  const auto trackEta 
    = std::atanh(params.momentum().z() / params.absoluteMomentum());
  const auto clusEta = std::atanh(clusZ / globalSL.norm());

  h_rphiResid->Fill(clusR / Acts::UnitConstants::cm, drphi);
  h_zResid->Fill(stateZ / Acts::UnitConstants::cm, dz);
  h_etaResid->Fill(trackEta, clusEta - trackEta);
  h_zResidLayer->Fill(clusR / Acts::UnitConstants::cm, dz);
  h_etaResidLayer->Fill(clusR / Acts::UnitConstants::cm, 
			clusEta - trackEta);

  const auto trackPPhi = -params.momentum()(0) * std::sin(statePhi) +
    params.momentum()(1) * std::cos(statePhi);
  const auto trackPR = params.momentum()(0) * std::cos(statePhi) +
    params.momentum()(1) * std::sin(statePhi);
  
  const auto trackPZ    = params.momentum()(1);
  const auto trackAlpha = -trackPPhi / trackPR;
  const auto trackBeta  = -trackPZ / trackPR;

  if(std::abs(trackAlpha) > m_maxTAlpha
     || std::abs(drphi) > m_maxResidual)
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
  auto *geomContainer = findNode::getClass<PHG4CylinderCellGeomContainer>(topNode, "CYLINDERCELLGEOM_SVTX");
  if(!geomContainer)
    {
      std::cout << PHWHERE << "No CYLINDERCELLGEOM_SVTX node, exiting."
		<< std::endl;
    }

  std::vector<Acts::Vector3D> delta(m_totalBins);
  std::vector<Acts::SymMatrix3D> cov(m_totalBins);

  for(int i = 0; i < m_totalBins; ++i)
    {
      cov[i] = m_lhs[i].inverse();
      delta[i] = m_lhs[i].partialPivLu().solve(m_rhs[i]);

      if(m_lhs[i].sum() != 0 || m_rhs[i].sum() !=0){
      std::cout << "m_lhs " <<i << std::endl;
      for(int j =0; j<m_lhs[i].rows(); j++)
	{
	  for(int k=0; k< m_lhs[i].cols(); k++)
	    {
	      std::cout << m_lhs[i](j,k) << ", ";
	    }
	  std::cout<<std::endl;
	}
      std::cout << "cov " << std::endl;
      for(int j =0; j<cov[i].rows(); j++)
	{
	  for(int k=0; k< cov[i].cols(); k++)
	    {
	      std::cout << cov[i](j,k) << ", ";
	    }
	  std::cout<<std::endl;
	}
      std::cout << "m_rhs " <<std::endl;
      for(int j =0; j<m_rhs[i].rows(); j++)
	{
	  for(int k=0; k< m_rhs[i].cols(); k++)
	    {
	      std::cout << m_rhs[i](j,k) << ", ";
	    }
	  std::cout<<std::endl;
	}
      std::cout << "delta " << std::endl;
      for(int j =0; j<delta[i].rows(); j++)
	{
	  for(int k=0; k< delta[i].cols(); k++)
	    {
	      std::cout << delta[i](j,k) << ", ";
	    }
	  std::cout<<std::endl;
	}
      }
    }
    

  /// Three dimensions for the matrices
  std::vector<std::unique_ptr<TGraphErrors>> 
    graphs(m_zBins * m_phiBins * m_nCoord);

  for(int iz = 0; iz < m_zBins; ++iz) {
    for(int iphi = 0; iphi < m_phiBins; ++iphi) {
      for(int icoord = 0; icoord < m_nCoord; ++icoord) {
	const int tgrIndex = iz + m_zBins * ( iphi + m_phiBins * icoord);
	
	graphs[tgrIndex].reset(new TGraphErrors());
	graphs[tgrIndex]->SetName(Form("tg_%i_%i_%i", iz, iphi, icoord));
	
	for(int ir = 0; ir < m_rBins; ++ir) {
	  /// Get TPC layers in sPHENIX, not Acts, coordinates, hence
	  /// add 7 to the value
	  const int innerLayer = 7 + m_nLayersTpc * ir / m_rBins;
	  const int outerLayer = 7 + m_nLayersTpc * (ir+1) / m_rBins - 1;
	  
	  const auto innerRadius = geomContainer->GetLayerCellGeom(innerLayer)->get_radius();
	  const auto outerRadius = geomContainer->GetLayerCellGeom(outerLayer)->get_radius();
	  const float r = (innerRadius + outerRadius) / 2.;
	  
	  int index = getCell(iz, ir, iphi);
	  if(!std::isnan(delta[index](icoord,0)) or
	     !std::isnan(std::sqrt(cov[index](icoord,icoord))))
	    std::cout << iz<< "  " << iphi << "  " << icoord << "  "<<ir << "  " << r << "  " << delta[index](icoord,0) << "  " << std::sqrt(cov[index](icoord, icoord)) << std::endl; 
	  graphs[tgrIndex]->SetPoint(ir, r, delta[index](icoord,0));
	  graphs[tgrIndex]->SetPointError(ir, 0, 
					 std::sqrt(cov[index](icoord,icoord)));

	}
      }
    }
  }

  TFile *outputFile = new TFile((Name() + "_distortions.root").c_str(), 
				"RECREATE");
  outputFile->cd();
  for(auto&& gr : graphs)
    gr->Write();
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

int PHTpcResiduals::getCell(const int iz, const int ir, const int iphi)
{
  if( ir < 0 || ir >= m_rBins ) return -1;
  if( iphi < 0 || iphi >= m_phiBins ) return -1;
  if( iz < 0 || iz >= m_zBins ) return -1;
  return iz + m_zBins*( ir + m_rBins*iphi );
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

