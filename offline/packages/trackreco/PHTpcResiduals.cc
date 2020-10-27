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

PHTpcResiduals::PHTpcResiduals(const std::string &name)
  : SubsysReco(name)
{
}


PHTpcResiduals::~PHTpcResiduals()
{
}

int PHTpcResiduals::Init(PHCompositeNode *topNode)
{
  outfile = new TFile(std::string(Name() + ".root").c_str(), 
		      "recreate");
  
  h_rphiResid = new TH2F("rphiResid",";r [cm]; #Deltar#phi [mm]",
			 60,20,80,50,-10,10);
  h_zResid = new TH2F("zResid",";z [cm]; #Deltaz [mm]",
		      200,-100,100,100,-100,100);
  h_etaResid = new TH2F("etaResid",";#eta;#Delta#eta",
			20,-1,1,50,-0.2,0.2);
  h_etaResidLayer = new TH2F("etaResidLayer",";r [cm]; #Delta#eta",
			     60,20,80,50,-0.2,0.2);
  h_zResidLayer = new TH2F("zResidLayer",";r [cm]; #Deltaz [mm]",
			   60,20,80,100,-100,100);
  return Fun4AllReturnCodes::EVENT_OK;

}

int PHTpcResiduals::InitRun(PHCompositeNode *topNode)
{
  if(getNodes(topNode) != Fun4AllReturnCodes::EVENT_OK)
    return Fun4AllReturnCodes::ABORTEVENT;

  return Fun4AllReturnCodes::EVENT_OK;
}

int PHTpcResiduals::process_event(PHCompositeNode *topNode)
{

  int returnVal = processTracks(topNode);

  return returnVal;
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
  const auto trackParams = track.getTrackParams();
  
  for(auto sl : sourceLinks)
    {
      /// Only analyze TPC 
      if(sl.referenceSurface().geometryId().volume() != 14)
	continue;
      
      auto result = propagateTrackState(trackParams, sl);
      if(result.ok())
	{
	  
	  auto trackStateParams = std::move(**result);
	  
	  calculateTpcResiduals(trackStateParams, sl);
	}
      else
	{
	  if(Verbosity() > 2)
	    std::cout << PHWHERE << "Propagation failed, check log"
		      << std::endl;
	  continue;
	}
    } 
  
}

BoundTrackParamPtrResult PHTpcResiduals::propagateTrackState(
			   const ActsExamples::TrackParameters& params,
			   const SourceLink& sl)
{
  
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
      
      auto logger = Acts::getDefaultLogger("PHTpcResiduals", 
					    Acts::Logging::INFO);
      
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

  

  return;
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

int PHTpcResiduals::End(PHCompositeNode *topNode)
{
  outfile->cd();
  h_rphiResid->Write();
  h_etaResid->Write();
  h_zResidLayer->Write();
  h_etaResidLayer->Write();
  h_zResid->Write();
  outfile->Write();
  outfile->Close();

  return Fun4AllReturnCodes::EVENT_OK;
}
