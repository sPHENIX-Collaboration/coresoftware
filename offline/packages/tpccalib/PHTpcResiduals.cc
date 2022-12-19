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
  if( m_savehistograms ) makeHistograms();
  
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

  // save histograms
  if( m_savehistograms && m_histogramfile )
  {
    m_histogramfile->cd();
    for( const auto o:std::initializer_list<TObject*>({ h_rphiResid, h_zResid, h_etaResidLayer, h_zResidLayer, h_etaResid, h_index, h_alpha, h_beta, h_deltarphi_layer, h_deltaz_layer, residTup }) )
    { if( o ) o->Write(); }
    
    // also write histograms from vectors
    for( const auto& [cell,h]:h_drphi ) { if(h) h->Write(); }
    for( const auto& [cell,h]:h_dz ) { if(h) h->Write(); }
    for( const auto& [cell,h]:h_drphi_alpha ) { if(h) h->Write(); }
    for( const auto& [cell,h]:h_dz_beta ) { if(h) h->Write(); }
    
    m_histogramfile->Close();
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
 
  if(track->get_pt() < 0.5)
    return false;

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
					    trackQ/p, cov).value();
 
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
  
  // create ACTS parameters from track parameters at origin
  auto trackParams = makeTrackParams(track);

  // store crossing. It is used in calculating cluster's global position
  m_crossing = track->get_crossing();
  assert( m_crossing != SHRT_MAX );
  
  for( const auto& key:get_cluster_keys( track ) )
  {
    
    // assign to static cluskey
    cluskey = key;

    // increment counter
    ++m_total_clusters;
    
    // make sure cluster is from TPC
    const auto detId = TrkrDefs::getTrkrId(cluskey);
    if(detId != TrkrDefs::tpcId) continue;  
    
    const auto cluster = m_clusterContainer->findCluster(cluskey);      
    const auto surf = m_tGeometry->maps().getSurface( cluskey, cluster );
    auto result = propagateTrackState(trackParams, surf);
    
    // skip if propagation failed
    if(std::isnan(result.first))
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
    auto pathLength = result.first;
    auto trackStateParams = result.second;
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
    
    addTrackState( track, pathLength, trackStateParams );
    
    // calculate residuals with respect to cluster
    // Get all the relevant information for residual calculation
    const auto globClusPos = getGlobalPosition(cluskey, cluster, m_crossing);
    clusR = get_r(globClusPos(0),globClusPos(1));
    clusPhi = std::atan2(globClusPos(1), globClusPos(0));
    clusZ = globClusPos(2);
    
    // cluster errors
    if( m_cluster_version >= 4 )
    {
      const auto errors_square = m_cluster_error_parametrization.get_cluster_error( track->get_tpc_seed(), cluster, clusR, cluskey ); 
      clusRPhiErr = std::sqrt( errors_square.first );
      clusZErr = std::sqrt( errors_square.second );
    } else {
      clusRPhiErr = cluster->getRPhiError();
      clusZErr = cluster->getZError();
    }
    
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
    if( m_cluster_version < 4 )
    {
      /*
       * remove clusters with too small errors since they are likely pathological
       * and have a large contribution to the chisquare
       * TODO: make these cuts configurable
       */
      if(clusRPhiErr < 0.015) continue;
      if(clusZErr < 0.05) continue;
    }
    
    const auto globalStatePos = trackStateParams.position(m_tGeometry->geometry().getGeoContext());
    const auto globalStateMom = trackStateParams.momentum();
    const auto globalStateCov = *trackStateParams.covariance();
    
    stateRPhiErr = std::sqrt(globalStateCov(Acts::eBoundLoc0, Acts::eBoundLoc0))/Acts::UnitConstants::cm;
    stateZErr = sqrt(globalStateCov(Acts::eBoundLoc1, Acts::eBoundLoc1))/Acts::UnitConstants::cm;
    
    stateZ = globalStatePos.z()/Acts::UnitConstants::cm;
    
    const auto globStateX = globalStatePos.x()/Acts::UnitConstants::cm;
    const auto globStateY = globalStatePos.y()/Acts::UnitConstants::cm;
    const auto globStateZ = stateZ;
    
    stateR = std::sqrt(square(globStateX) + square(globStateY));
    
    const auto dr = clusR - stateR;
    const auto trackDrDt = (globStateX * globalStateMom(0) + globStateY * globalStateMom(1)) / stateR;
    const auto trackDxDr = globalStateMom(0) / trackDrDt;
    const auto trackDyDr = globalStateMom(1) / trackDrDt;
    const auto trackDzDr = globalStateMom(2) / trackDrDt;
    
    const auto trackX = globStateX + dr * trackDxDr;
    const auto trackY = globStateY + dr * trackDyDr;
    const auto trackZ = globStateZ + dr * trackDzDr;
    
    if(Verbosity() > 2)
    {
      std::cout << "PHTpcResiduals::processTrack -"
        << " stateR: " << stateR 
        << " dr: " << dr 
        << " trackDrDt: " << trackDrDt 
        << " trackDxDr: " << trackDxDr
        << " trackDyDr: " << trackDyDr 
        << " trackDzDr: " << trackDzDr
        << " track position: (" << trackX << ", " << trackY << ", " << trackZ  << ")"
        << std::endl;
    }
    
    statePhi = std::atan2(trackY, trackX);
    stateZ = trackZ;
    
    if(Verbosity() > 3)
    {
      std::cout << "PHTpcResiduals::processTrack -" 
        << " stateR: " << stateR
        << " statePhi: " << statePhi << "+/-" << stateRPhiErr
        << " stateZ: " << stateZ << "+/-" << stateZErr 
        << std::endl;
    }
    
    const auto erp = square(clusRPhiErr) + square(stateRPhiErr);
    const auto ez = square(clusZErr) + square(stateZErr);
    
    // Calculate residuals
    drphi = clusR * deltaPhi(clusPhi - statePhi);
    dz  = clusZ - stateZ;
    
    if(Verbosity() > 3)
    {
      std::cout << "PHTpcResiduals::processTrack -"
        << " drphi: " << drphi 
        << " dz: " << dz 
        << std::endl;
    }
    
    const auto trackEta = std::atanh(trackStateParams.momentum().z() / trackStateParams.absoluteMomentum());
    const auto clusEta = std::atanh(clusZ / std::sqrt(
      square(globClusPos(0)) +
      square(globClusPos(1)) +
      square(globClusPos(2))));
    
    const auto trackPPhi = -trackStateParams.momentum()(0) * std::sin(statePhi) + trackStateParams.momentum()(1) * std::cos(statePhi);
    const auto trackPR = trackStateParams.momentum()(0) * std::cos(statePhi) + trackStateParams.momentum()(1) * std::sin(statePhi);
    const auto trackPZ = trackStateParams.momentum()(2);
    
    const auto trackAlpha = -trackPPhi / trackPR;
    const auto trackBeta = -trackPZ / trackPR;
    
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
    
    tanBeta = trackBeta;
    tanAlpha = trackAlpha;
    
    // get cell index
    const auto index = getCell(globClusPos);
    if(Verbosity() > 3)
    { std::cout << "Bin index found is " << index << std::endl; }
    
    if(index < 0 ) continue;
    
    if(Verbosity() > 3)
    {
      std::cout << "PHTpcResiduals::processTrack - layer: " << (int) TrkrDefs::getLayer(key) << std::endl;
      std::cout << "PHTpcResiduals::processTrack -"
        << " cluster: (" << clusR << ", " << clusR*clusPhi << ", " << clusZ << ")"
        << " (" << clusRPhiErr << ", " << clusZErr << ")"
        << std::endl;
      
      std::cout << "PHTpcResiduals::processTrack -"
        << " track: (" << stateR << ", " << clusR*statePhi << ", " << stateZ << ")"
        << " (" << tanAlpha << ", " << tanBeta << ")"
        << " (" << stateRPhiErr << ", " << stateZErr << ")"
        << std::endl;
      std::cout << std::endl;
    }
    
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
      
      const auto layer =  TrkrDefs::getLayer(key);
      h_deltarphi_layer->Fill( layer, drphi );
      h_deltaz_layer->Fill( layer, dz );
      
      { const auto iter = h_drphi.find( index ); if( iter != h_drphi.end() ) iter->second->Fill( drphi ); }
      { const auto iter = h_drphi_alpha.find( index ); if( iter != h_drphi_alpha.end() ) iter->second->Fill( tanAlpha, drphi ); }
      { const auto iter = h_dz.find( index ); if( iter != h_dz.end() ) iter->second->Fill( dz ); }
      { const auto iter = h_dz_beta.find( index ); if( iter != h_dz_beta.end() ) iter->second->Fill( tanBeta, dz ); }
      
      residTup->Fill();
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

//__________________________________________________________________________________________________________________________________________
PHTpcResiduals::BoundTrackParamPair  PHTpcResiduals::propagateTrackState( const Acts::BoundTrackParameters& params, const Surface& surf) const
{
 

  using Stepper = Acts::EigenStepper<>;
  using Propagator = Acts::Propagator<Stepper, Acts::Navigator>;

  Stepper stepper(m_tGeometry->geometry().magField);
  Acts::Navigator::Config cfg{m_tGeometry->geometry().tGeometry};
  Acts::Navigator navigator(cfg);
  Propagator propagator(stepper, navigator);

  Acts::Logging::Level logLevel = Acts::Logging::FATAL;
  if(Verbosity() > 10)
    logLevel = Acts::Logging::VERBOSE;

  auto logger = Acts::getDefaultLogger("PHTpcResiduals", logLevel);
      
  Acts::PropagatorOptions<> options(m_tGeometry->geometry().getGeoContext(),
				    m_tGeometry->geometry().magFieldContext,
				    Acts::LoggerWrapper{*logger});
     
  auto result = propagator.propagate(params, *surf, options);
   
  if(result.ok())
  {
    const Acts::BoundTrackParameters params = *result.value().endParameters;
    double pathlength = result.value().pathLength / Acts::UnitConstants::cm;
    // return both path length and extrapolated parameters
    return std::make_pair( pathlength, params );
  } else {
    return std::make_pair(NAN, Acts::BoundTrackParameters(nullptr, Acts::BoundVector::Zero(),-1));

  }
}

//_______________________________________________________________________________________________________
void PHTpcResiduals::addTrackState( SvtxTrack* track, float pathlength, const Acts::BoundTrackParameters& params )
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

//_______________________________________________________________________________
void PHTpcResiduals::makeHistograms()
{  
  
  m_histogramfile.reset( new TFile(m_histogramfilename.c_str(), "RECREATE") );
  m_histogramfile->cd();

  const auto total_bins = m_matrix_container->get_grid_size();  
  h_beta = new TH2F("betadz",";tan#beta; #Deltaz [cm]",100,-0.5,0.5,100,-0.5,0.5);
  h_alpha = new TH2F("alphardphi",";tan#alpha; r#Delta#phi [cm]", 100,-0.5,0.5,100,-0.5,0.5);
  h_index = new TH1F("index",";index",total_bins, 0, total_bins);
  h_rphiResid = new TH2F("rphiResid", ";r [cm]; #Deltar#phi [cm]", 60, 20, 80, 500, -2, 2);
  h_zResid = new TH2F("zResid", ";z [cm]; #Deltaz [cm]", 200, -100, 100, 1000, -2, 2);
  h_etaResid = new TH2F("etaResid", ";#eta;#Delta#eta",	20, -1, 1, 500, -0.2, 0.2);
  h_etaResidLayer = new TH2F("etaResidLayer", ";r [cm]; #Delta#eta", 60, 20, 80, 500, -0.2, 0.2);
  h_zResidLayer = new TH2F("zResidLayer", ";r [cm]; #Deltaz [cm]", 60, 20, 80, 1000, -2, 2);
  h_deltarphi_layer = new TH2F( "deltarphi_layer", ";layer; r.#Delta#phi_{track-cluster} (cm)", 57, 0, 57, 500, -2, 2 );
  h_deltaz_layer = new TH2F( "deltaz_layer", ";layer; #Deltaz_{track-cluster} (cm)", 57, 0, 57, 100, -2, 2 );

  {

    // get grid dimensions from matrix container
    int phibins = 0;
    int rbins = 0;
    int zbins = 0;
    m_matrix_container->get_grid_dimensions( phibins, rbins, zbins );
    
    // get bins corresponding to selected angles
    std::set<int> phibin_rec;
    std::transform( phi_rec.begin(), phi_rec.end(), std::inserter( phibin_rec, phibin_rec.end() ), [&]( const float& phi ) { return phibins*(phi-m_phiMin)/(m_phiMax-m_phiMin); } );
    
    std::set<int> zbin_rec;
    std::transform( z_rec.begin(), z_rec.end(), std::inserter( zbin_rec, zbin_rec.end() ), [&]( const float& z ) { return zbins*(z-m_zMin)/(m_zMax-m_zMin); } );

    // keep track of all cell ids that match selected histograms
    for( int iphi = 0; iphi < phibins; ++iphi )
      for( int ir = 0; ir < rbins; ++ir )
      for( int iz = 0; iz < zbins; ++iz )
    {
    
      if( phibin_rec.find( iphi ) == phibin_rec.end() || zbin_rec.find( iz ) == zbin_rec.end() ) continue;
      const auto icell = m_matrix_container->get_cell_index( iphi, ir, iz );
      
      {
        // rphi residuals
        const auto hname = Form( "residual_drphi_p%i_r%i_z%i", iphi, ir, iz );
        auto h = new TH1F( hname, hname, 100, -m_maxResidualDrphi, +m_maxResidualDrphi );
        h->GetXaxis()->SetTitle( "r.#Delta#phi_{cluster-track} (cm)" );
        h_drphi.insert( std::make_pair( icell, h ) );
      }

      {
        // 2D histograms
        const auto hname = Form( "residual_2d_drphi_p%i_r%i_z%i", iphi, ir, iz );
        auto h = new TH2F( hname, hname, 100, -m_maxTAlpha, m_maxTAlpha, 100, -m_maxResidualDrphi, +m_maxResidualDrphi );
        h->GetXaxis()->SetTitle( "tan#alpha" );
        h->GetYaxis()->SetTitle( "r.#Delta#phi_{cluster-track} (cm)" );
        h_drphi_alpha.insert( std::make_pair( icell, h ) );
      }

      {
        // z residuals
        const auto hname = Form( "residual_dz_p%i_r%i_z%i", iphi, ir, iz );
        auto h = new TH1F( hname, hname, 100, -m_maxResidualDz, +m_maxResidualDz );
        h->GetXaxis()->SetTitle( "#Deltaz_{cluster-track} (cm)" );
        h_dz.insert( std::make_pair( icell, h ) );
      }

      {
        // 2D histograms
        static constexpr double maxTBeta = 0.5;
        const auto hname = Form( "residual_2d_dz_p%i_r%i_z%i", iphi, ir, iz );
        auto h = new TH2F( hname, hname, 100, -maxTBeta, maxTBeta, 100, -m_maxResidualDz, +m_maxResidualDz );
        h->GetXaxis()->SetTitle( "tan#beta" );
        h->GetYaxis()->SetTitle( "#Deltaz_{cluster-track} (cm)" );
        h_dz_beta.insert( std::make_pair( icell, h ) );
      }     
    }
    
  }
  
  residTup = new TTree("residTree","tpc residual info");
  residTup->Branch("tanAlpha",&tanAlpha,"tanAlpha/D");
  residTup->Branch("tanBeta",&tanBeta,"tanBeta/D");
  residTup->Branch("drphi",&drphi,"drphi/D");
  residTup->Branch("dz",&dz,"dz/D");
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
  residTup->Branch("cluskey",&cluskey,"cluskey/l");
  residTup->Branch("event",&m_event,"event/I");

}

void PHTpcResiduals::setGridDimensions(const int phiBins, const int rBins, const int zBins)
{ m_matrix_container->set_grid_dimensions( phiBins, rBins, zBins ); }

