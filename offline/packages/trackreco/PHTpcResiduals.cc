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

  int returnVal = getTpcResiduals(topNode);

  return returnVal;
}

int PHTpcResiduals::getTpcResiduals(PHCompositeNode *topNode)
{

  std::map<unsigned int, ActsTrack>::iterator trackIter;
  for(trackIter = m_actsProtoTracks->begin();
      trackIter != m_actsProtoTracks->end();
      ++trackIter)
    {
      auto track = trackIter->second;
      auto sourceLinks = track.getSourceLinks();
      const auto trackParams = track.getTrackParams();
      const auto momentumVec = trackParams.momentum();
      
      calculateTpcResiduals(sourceLinks, momentumVec);
      

    }
  
  return Fun4AllReturnCodes::EVENT_OK;
}

void PHTpcResiduals::calculateTpcResiduals(const std::vector<SourceLink> sourceLinks,
					   const Acts::Vector3D momentum)
{

  for(auto sl : sourceLinks)
    {
      auto volume = sl.referenceSurface().geometryId().volume();
	 
      /// Only care about TPC source links
      if(volume != 14)
	continue;
      
      /// Source link returns a 3D vector, and we need a 2D vector for 
      /// the acts transformations
      auto localPos  = sl.location();
      Acts::Vector2D local(localPos.x(), localPos.y());

      
      auto globalPos = sl.referenceSurface().localToGlobal(
                       m_tGeometry->geoContext,
		       local,
		       momentum);
  
      double measPhi = atan2(globalPos.y(), globalPos.x());
      double measR = sqrt(globalPos.x() * globalPos.x() + 
			       globalPos.y() * globalPos.y());
      double measRPhi = measPhi * measR;
      double siMMFitPhi = atan2(momentum.y(), momentum.x());
      
      /// Measurement is confined to a surface in R, so evaluate
      /// the rphi of the track at this radius
      double fitRPhi = measR * siMMFitPhi;
   
      double measZ = globalPos.z();
      
      double siMMFitTheta = acos(momentum.z() / momentum.norm());
      
      /// Need to convert theta from beam axis to theta from beam
      /// perpendicular axis
      siMMFitTheta -= M_PI / 2.;

      /// Need negative sign to correct for flipped theta convention
      double fitZ = -1 * measR * tan(siMMFitTheta);
      double rphiResid = measRPhi - fitRPhi;
      double zResid = measZ - fitZ;

      double measEta = atanh(globalPos.z() / globalPos.norm());
      double fitEta = atanh(momentum.z() / momentum.norm());

      double etaResid = measEta - fitEta;

      if(Verbosity() > 0)
	{
	  std::cout << "Fit eta/phi : " 
		    << fitEta << "  " << siMMFitPhi
		    << std::endl;
	  std::cout << "Fit theta : " << siMMFitTheta << std::endl;
	  std::cout << "Fit eta from theta : " 
		    << 2. * atan(exp(-fitEta)) << std::endl;
	  std::cout << "Meas, fit Z : " << measZ << ", " << fitZ 
		    << std::endl;
	  std::cout << "rphi, z, eta residuals: " << rphiResid << ", "
		    << zResid << " mm" << "  " << etaResid << std::endl;
	}

      h_rphiResid->Fill(measR / Acts::UnitConstants::cm, rphiResid);
      h_zResid->Fill(measZ / Acts::UnitConstants::cm, zResid);
      h_etaResid->Fill(measEta, etaResid);
      h_etaResidLayer->Fill(measR / Acts::UnitConstants::cm, etaResid);
      h_zResidLayer->Fill(measR / Acts::UnitConstants::cm, zResid);
    }

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
