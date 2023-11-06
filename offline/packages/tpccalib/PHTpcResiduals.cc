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

#include <tpc/TpcDistortionCorrectionContainer.h>

#include <trackbase/TpcDefs.h>
#include <trackbase/TrkrCluster.h>
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase_historic/SvtxTrack.h>
#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/SvtxTrackState_v1.h>

#include <trackreco/ActsPropagator.h>

#include <micromegas/MicromegasDefs.h>

#include <Acts/Surfaces/PerigeeSurface.hpp>
#include <Acts/Geometry/GeometryIdentifier.hpp>
#include <Acts/MagneticField/ConstantBField.hpp>
#include <Acts/MagneticField/InterpolatedBFieldMap.hpp>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#include <Acts/Propagator/EigenStepper.hpp>
#pragma GCC diagnostic pop

#include <Acts/Propagator/Navigator.hpp>
#include <Acts/Surfaces/Surface.hpp>

#include <Acts/MagneticField/MagneticFieldProvider.hpp>


#include <cmath>
#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TH2.h>

#include <iostream>
#include <sstream>

namespace 
{
  
  // square
  template<class T> inline constexpr T square( const T& x ) { return x*x; }

  // radius
  template<class T> T get_r( const T& x, const T& y ) { return std::sqrt( square(x) + square(y) ); }

  template<class T> inline constexpr T deltaPhi(const T& phi)
  {
    if (phi > M_PI) 
      return phi - 2. * M_PI;
    else if (phi <= -M_PI) 
      return phi + 2.* M_PI;
    else 
      return phi;
  }
  
  /// get sector median angle associated to a given index
  /** this assumes that sector 0 is centered on phi=0, then numbered along increasing phi */
  inline constexpr double get_sector_phi( int isec ) 
  { return isec*M_PI/6; }
  
  // specify bins for which one will save histograms
  static const std::vector<float> phi_rec = { get_sector_phi(9) };
  static const std::vector<float> z_rec = { 5. };

  //! get cluster keys from a given track
  std::vector<TrkrDefs::cluskey> get_cluster_keys( SvtxTrack* track )
  {
    std::vector<TrkrDefs::cluskey> out;
    for( const auto& seed: { track->get_silicon_seed(), track->get_tpc_seed() } )
    {
      if( seed )
      { std::copy( seed->begin_cluster_keys(), seed->end_cluster_keys(), std::back_inserter( out ) ); }
    }
    return out;
  }

  /// return number of clusters of a given type that belong to a tracks
  template<int type>
    int count_clusters( const std::vector<TrkrDefs::cluskey>& keys )
  {
    return std::count_if( keys.begin(), keys.end(),
      []( const TrkrDefs::cluskey& key ) { return TrkrDefs::getTrkrId(key) == type; } );
  }

}

//___________________________________________________________________________________
PHTpcResiduals::PHTpcResiduals(const std::string &name)
  : SubsysReco(name)
  , m_matrix_container( new TpcSpaceChargeMatrixContainerv1 )
{}

//___________________________________________________________________________________
int PHTpcResiduals::Init(PHCompositeNode */*topNode*/)
{

  // configuration printout
  std::cout << "PHTpcResiduals::Init - m_maxTAlpha: " << m_maxTAlpha << std::endl;
  std::cout << "PHTpcResiduals::Init - m_maxTBeta: " << m_maxTBeta << std::endl;
  std::cout << "PHTpcResiduals::Init - m_maxResidualDrphi: " << m_maxResidualDrphi << " cm" << std::endl;
  std::cout << "PHTpcResiduals::Init - m_maxResidualDz: " << m_maxResidualDz << " cm" << std::endl;
  std::cout << "PHTpcResiduals::Init - m_minPt: " << m_minPt << " GeV/c" << std::endl;
  
  // reset counters
  m_total_tracks = 0;
  m_accepted_tracks = 0;

  m_total_clusters = 0;
  m_accepted_clusters = 0;

  return Fun4AllReturnCodes::EVENT_OK;
}

//___________________________________________________________________________________
int PHTpcResiduals::InitRun(PHCompositeNode *topNode)
{
  if(getNodes(topNode) != Fun4AllReturnCodes::EVENT_OK)
  { return Fun4AllReturnCodes::ABORTEVENT; }

  if(createNodes(topNode) != Fun4AllReturnCodes::EVENT_OK)
  { return Fun4AllReturnCodes::ABORTEVENT; }

  return Fun4AllReturnCodes::EVENT_OK;
}

//___________________________________________________________________________________
int PHTpcResiduals::process_event(PHCompositeNode *topNode)
{
  const auto returnVal = processTracks(topNode);  
  ++m_event;

  return returnVal;
}

//___________________________________________________________________________________
int PHTpcResiduals::End(PHCompositeNode */*topNode*/)
{
  std::cout << "PHTpcResiduals::End - writing matrices to " << m_outputfile << std::endl;
      
  // save matrix container in output file
  if( m_matrix_container )
  {
    std::unique_ptr<TFile> outputfile( TFile::Open( m_outputfile.c_str(), "RECREATE" ) );
    outputfile->cd();
    m_matrix_container->Write( "TpcSpaceChargeMatrixContainer" );
  }
  
  // print counters
  std::cout
    << "PHTpcResiduals::End -"
    << " track statistics total: " << m_total_tracks
    << " accepted: " << m_accepted_tracks
    << " fraction: " << 100.*m_accepted_tracks/m_total_tracks << "%"
    << std::endl;

  std::cout
    << "PHTpcResiduals::End -"
    << " cluster statistics total: " << m_total_clusters
    << " accepted: " << m_accepted_clusters << " fraction: "
    << 100.*m_accepted_clusters/m_total_clusters << "%"
    << std::endl;

  return Fun4AllReturnCodes::EVENT_OK;
}

//___________________________________________________________________________________
int PHTpcResiduals::processTracks(PHCompositeNode */*topNode*/)
{

  if( Verbosity() )
  { std::cout << "PHTpcResiduals::processTracks - proto track size " << m_trackMap->size() <<std::endl; }

  for(const auto &[trackKey, track] : *m_trackMap)
  {
    if(Verbosity() > 1) 
    { std::cout << "PHTpcResiduals::processTracks - Processing track key " << trackKey << std::endl; }
    
    ++m_total_tracks;
    if(checkTrack(track))
    {
      ++m_accepted_tracks;
      processTrack(track);
    }
  }
  
  return Fun4AllReturnCodes::EVENT_OK;
}

//___________________________________________________________________________________
bool PHTpcResiduals::checkTrack(SvtxTrack* track) const
{
  if(Verbosity() > 2)
  { std::cout << "PHTpcResiduals::checkTrack - pt: " << track->get_pt() << std::endl; }
 
  if(track->get_pt() < m_minPt)
  { return false; }

  // ignore tracks with too few mvtx, intt and micromegas hits
  const auto cluster_keys( get_cluster_keys( track ) );
  if( count_clusters<TrkrDefs::mvtxId>(cluster_keys) < 2 ) return false;
  if( count_clusters<TrkrDefs::inttId>(cluster_keys) < 2 ) return false;
  if( m_useMicromegas && count_clusters<TrkrDefs::micromegasId>(cluster_keys) < 2 ) return false;


  return true;

}

//___________________________________________________________________________________
Acts::BoundTrackParameters PHTpcResiduals::makeTrackParams(SvtxTrack* track) const
{
  Acts::Vector3 momentum(track->get_px(), 
			 track->get_py(), 
			 track->get_pz());
  double trackQ = track->get_charge() * Acts::UnitConstants::e;
  double p = track->get_p();

  /* get acts covariance matrix from track parameters */
  const auto cov = m_transformer.rotateSvtxTrackCovToActs( track );
  
  /* get position from  track parameters */
  const Acts::Vector3 position(track->get_x() * Acts::UnitConstants::cm,
    track->get_y() * Acts::UnitConstants::cm,
    track->get_z() * Acts::UnitConstants::cm);

  const auto perigee = Acts::Surface::makeShared<Acts::PerigeeSurface>(position);
  const auto actsFourPos = Acts::Vector4(position(0), position(1), position(2), 10 * Acts::UnitConstants::ns);

  return Acts::BoundTrackParameters::create(perigee, m_tGeometry->geometry().getGeoContext(),
					    actsFourPos, momentum,
					    trackQ/p, cov,
					    Acts::ParticleHypothesis::pion()).value();
 
}

void PHTpcResiduals::processTrack(SvtxTrack* track)
{

  if(Verbosity() > 1)
  {
    std::cout << "PHTpcResiduals::processTrack -" 
      << " track momentum: " << track->get_p() 
      << " position: (" << track->get_x() << ", " << track->get_y() << ", " << track->get_z() << ")"
      << std::endl;
  }
  ActsPropagator propagator(m_tGeometry);

  // create ACTS parameters from track parameters at origin
  auto trackParams = makeTrackParams(track);
  
  // store crossing. It is used in calculating cluster's global position
  m_crossing = track->get_crossing();
  assert( m_crossing != SHRT_MAX );

  for( const auto& cluskey:get_cluster_keys( track ) )
  {
    // increment counter
    ++m_total_clusters;
    
    // make sure cluster is from TPC
    const auto detId = TrkrDefs::getTrkrId(cluskey);
    if(detId != TrkrDefs::tpcId) continue;  
    
    const auto cluster = m_clusterContainer->findCluster(cluskey);      
    const auto surface = m_tGeometry->maps().getSurface( cluskey, cluster );   
    auto result = propagator.propagateTrack(trackParams, surface);
    
    // skip if propagation failed
    if(!result.ok())
    {
      if( Verbosity() > 1 )
      {
        std::cout << "Starting track params position/momentum: "
          << trackParams.position(m_tGeometry->geometry().geoContext).transpose()
          << std::endl
          << std::endl
          << "Track params phi/eta "
          << std::atan2(trackParams.momentum().y(),
          trackParams.momentum().x())
          << " and "
          << std::atanh(trackParams.momentum().z()/trackParams.momentum().norm())
          << std::endl;
      }
      
      continue;
    }
    
    // get extrapolated track state, convert to sPHENIX and add to track
    auto& [pathLength, trackStateParams] = result.value();    
    pathLength /= Acts::UnitConstants::cm;

    if(Verbosity() > 1)
    {
      std::cout << "PHTpcResiduals::processTrack -"
        << " path length: " << pathLength
        << " track momentum : "
        << trackParams.momentum()
        << " propagator momentum : "
        << trackStateParams.momentum()
        << std::endl;
    }


    addTrackState( track, cluskey, pathLength, trackStateParams );
    
    // calculate residuals with respect to cluster
    // Get all the relevant information for residual calculation
    const auto globClusPos = getGlobalPosition(cluskey, cluster, m_crossing);
    const double clusR = get_r(globClusPos(0),globClusPos(1));
    const double clusPhi = std::atan2(globClusPos(1), globClusPos(0));
    const double clusZ = globClusPos(2);
    
    // cluster errors
    double clusRPhiErr = 0;
    double clusZErr = 0;
 
    clusRPhiErr = cluster->getRPhiError();
    clusZErr = cluster->getZError();
    
    
    if(Verbosity() > 3)
    {
      std::cout << "PHTpcResiduals::processTrack -"
        << " cluskey: " << cluskey
        << " clusR: " << clusR 
        << " clusPhi: " << clusPhi << "+/-" << clusRPhiErr
        << " clusZ: " << clusZ << "+/-" << clusZErr
        << std::endl;
    }
    
    /* 
     * as instructed by Christof, it should not be necessary to cut on small
     * cluster errors any more with clusters of version 4 or higher 
     */ 
    
    const auto globalStatePos = trackStateParams.position(m_tGeometry->geometry().getGeoContext());
    const auto globalStateMom = trackStateParams.momentum();
    const auto globalStateCov = *trackStateParams.covariance();
    
    const double trackRPhiErr = std::sqrt(globalStateCov(Acts::eBoundLoc0, Acts::eBoundLoc0))/Acts::UnitConstants::cm;
    const double trackZErr = sqrt(globalStateCov(Acts::eBoundLoc1, Acts::eBoundLoc1))/Acts::UnitConstants::cm;    
    
    const double globStateX = globalStatePos.x()/Acts::UnitConstants::cm;
    const double globStateY = globalStatePos.y()/Acts::UnitConstants::cm;
    const double globStateZ = globalStatePos.z()/Acts::UnitConstants::cm;
    
    const double trackR = std::sqrt(square(globStateX) + square(globStateY));
    
    const double dr = clusR - trackR;
    const double trackDrDt = (globStateX * globalStateMom(0) + globStateY * globalStateMom(1)) / trackR;
    const double trackDxDr = globalStateMom(0) / trackDrDt;
    const double trackDyDr = globalStateMom(1) / trackDrDt;
    const double trackDzDr = globalStateMom(2) / trackDrDt;
    
    const double trackX = globStateX + dr * trackDxDr;
    const double trackY = globStateY + dr * trackDyDr;
    const double trackZ = globStateZ + dr * trackDzDr;
    const double trackPhi = std::atan2(trackY, trackX);
    
    if(Verbosity() > 2)
    {
      std::cout << "PHTpcResiduals::processTrack -"
        << " trackR: " << trackR 
        << " dr: " << dr 
        << " trackDrDt: " << trackDrDt 
        << " trackDxDr: " << trackDxDr
        << " trackDyDr: " << trackDyDr 
        << " trackDzDr: " << trackDzDr
        << " trackPhi: " << trackPhi << "+/-" << trackRPhiErr
        << " track position: (" << trackX << ", " << trackY << ", " << trackZ  << ")"
        << std::endl;
    }
        
    const double erp = square(clusRPhiErr) + square(trackRPhiErr);
    const double ez = square(clusZErr) + square(trackZErr);
    
    // Calculate residuals
    const double drphi = clusR * deltaPhi(clusPhi - trackPhi);
    const double dz  = clusZ - trackZ;
    
    if(Verbosity() > 3)
    {
      std::cout << "PHTpcResiduals::processTrack -"
        << " drphi: " << drphi 
        << " dz: " << dz 
        << std::endl;
    }
        
    const double trackPPhi = -trackStateParams.momentum()(0) * std::sin(trackPhi) + trackStateParams.momentum()(1) * std::cos(trackPhi);
    const double trackPR = trackStateParams.momentum()(0) * std::cos(trackPhi) + trackStateParams.momentum()(1) * std::sin(trackPhi);
    const double trackPZ = trackStateParams.momentum()(2);
    
    const double trackAlpha = -trackPPhi / trackPR;
    const double trackBeta = -trackPZ / trackPR;
    
    if(Verbosity() > 3)
    {
      std::cout 
        << "PHTpcResiduals::processTrack -"
        << " trackPPhi: " << trackPPhi 
        << " trackPR: " << trackPR
        << " trackPZ: " << trackPZ 
        << " trackAlpha: " << trackAlpha 
        << " trackBeta: " << trackBeta
        << std::endl;
    }
        
    // get cell index
    const auto index = getCell(globClusPos);
    if(Verbosity() > 3)
    { std::cout << "Bin index found is " << index << std::endl; }
    
    if(index < 0 ) continue;
    
    if(Verbosity() > 3)
    {
      std::cout << "PHTpcResiduals::processTrack - layer: " << (int) TrkrDefs::getLayer(cluskey) << std::endl;
      std::cout << "PHTpcResiduals::processTrack -"
        << " cluster: (" << clusR << ", " << clusR*clusPhi << ", " << clusZ << ")"
        << " (" << clusRPhiErr << ", " << clusZErr << ")"
        << std::endl;
      
      std::cout << "PHTpcResiduals::processTrack -"
        << " track: (" << trackR << ", " << clusR*trackPhi << ", " << trackZ << ")"
        << " (" << trackAlpha << ", " << trackBeta << ")"
        << " (" << trackRPhiErr << ", " << trackZErr << ")"
        << std::endl;
      std::cout << std::endl;
    }
    
    // check track angles and residuals agains cuts
    if(std::abs(trackAlpha) > m_maxTAlpha || std::abs(drphi) > m_maxResidualDrphi)
    { continue; }
    
    if(std::abs(trackBeta) > m_maxTBeta || std::abs(dz) > m_maxResidualDz)
    { continue; }
    
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
    
    // increment number of accepted clusters
    ++m_accepted_clusters;
      
  }
        
}

//_______________________________________________________________________________________________________
void PHTpcResiduals::addTrackState( SvtxTrack* track, TrkrDefs::cluskey key, float pathlength, const Acts::BoundTrackParameters& params )
{

  /* this is essentially a copy of the code from trackbase_historic/ActsTransformations::fillSvtxTrackStates */
  
  // create track state
  SvtxTrackState_v1 state( pathlength );

  // save global position
  const auto global = params.position(m_tGeometry->geometry().getGeoContext());
  state.set_x(global.x() / Acts::UnitConstants::cm);
  state.set_y(global.y() / Acts::UnitConstants::cm);
  state.set_z(global.z() / Acts::UnitConstants::cm);

  // save momentum
  const auto momentum = params.momentum();
  state.set_px(momentum.x());
  state.set_py(momentum.y());
  state.set_pz(momentum.z());

  // covariance
  const auto globalCov = m_transformer.rotateActsCovToSvtxTrack(params);
  for (int i = 0; i < 6; ++i)
    for (int j = 0; j < 6; ++j)
  { state.set_error(i, j, globalCov(i,j)); }
 
  state.set_name(std::to_string((TrkrDefs::cluskey) key));

  track->insert_state(&state);
}

//_______________________________________________________________________________
int PHTpcResiduals::getCell(const Acts::Vector3& loc)
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

//_______________________________________________________________________________
int PHTpcResiduals::createNodes(PHCompositeNode */*topNode*/)
{ return Fun4AllReturnCodes::EVENT_OK; }

//_______________________________________________________________________________
int PHTpcResiduals::getNodes(PHCompositeNode *topNode)
{
  m_clusterContainer = findNode::getClass<TrkrClusterContainer>(topNode,"TRKR_CLUSTER");
  if(!m_clusterContainer)
  {
    std::cout << PHWHERE << "No TRKR_CLUSTER node on node tree. Exiting." << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  m_tGeometry = findNode::getClass<ActsGeometry>(topNode, "ActsGeometry");
  if(!m_tGeometry)
  {
    std::cout << "ActsTrackingGeometry not on node tree. Exiting." << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  m_trackMap = findNode::getClass<SvtxTrackMap>(topNode, "SvtxSiliconMMTrackMap");  
  if (!m_trackMap)
  {
    std::cout << PHWHERE << "SvtxSiliconMMTrackMap not on node tree. Exiting." << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  // tpc distortion corrections
  m_dcc_static = findNode::getClass<TpcDistortionCorrectionContainer>(topNode,"TpcDistortionCorrectionContainerStatic");
  m_dcc_average = findNode::getClass<TpcDistortionCorrectionContainer>(topNode,"TpcDistortionCorrectionContainerAverage");
  m_dcc_fluctuation = findNode::getClass<TpcDistortionCorrectionContainer>(topNode,"TpcDistortionCorrectionContainerFluctuation");

  return Fun4AllReturnCodes::EVENT_OK;
}

//_________________________________________________________________________________
Acts::Vector3 PHTpcResiduals::getGlobalPosition( TrkrDefs::cluskey key, TrkrCluster* cluster, short int crossing ) const
{
  // get global position from Acts transform
  auto globalPosition = m_tGeometry->getGlobalPosition(key, cluster);
  
  // for the TPC calculate the proper z based on crossing and side
  const auto trkrid = TrkrDefs::getTrkrId(key);
  if(trkrid ==  TrkrDefs::tpcId)
  {	 
    const auto side = TpcDefs::getSide(key);
    globalPosition.z() = m_clusterCrossingCorrection.correctZ(globalPosition.z(), side, crossing);    
    
    // apply distortion corrections
    if(m_dcc_static) 
    {
      globalPosition = m_distortionCorrection.get_corrected_position( globalPosition, m_dcc_static ); 
    }
    
    if(m_dcc_average) 
    { 
      globalPosition = m_distortionCorrection.get_corrected_position( globalPosition, m_dcc_average ); 
    }
    
    if(m_dcc_fluctuation) 
    { 
      globalPosition = m_distortionCorrection.get_corrected_position( globalPosition, m_dcc_fluctuation ); 
    }
  }
    
  return globalPosition;
}

void PHTpcResiduals::setGridDimensions(const int phiBins, const int rBins, const int zBins)
{ m_matrix_container->set_grid_dimensions( phiBins, rBins, zBins ); }

