#include "KshortReconstruction.h"

#include <ffaobjects/EventHeader.h>

#include <trackbase/TrkrDefs.h>
#include <trackbase/ActsGeometry.h>

#include <trackbase_historic/ActsTransformations.h>
#include <trackbase_historic/SvtxTrackMap_v2.h>

#include <trackreco/ActsPropagator.h>

#include <globalvertex/SvtxVertexMap.h>

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>

#include <Acts/Surfaces/CylinderSurface.hpp>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wshadow"
#include <ActsExamples/EventData/Trajectories.hpp>
#pragma GCC diagnostic pop

#include <TFile.h>
#include <TH1.h>
#include <TLorentzVector.h>
#include <TNtuple.h>
#include <TSystem.h>

#include <cmath>
#include <utility>

using BoundTrackParam = const Acts::BoundTrackParameters;
using BoundTrackParamResult = Acts::Result<BoundTrackParam>;
using SurfacePtr = std::shared_ptr<const Acts::Surface>;
using Trajectory = ActsExamples::Trajectories;


int KshortReconstruction::process_event(PHCompositeNode* topNode)
{

  EventHeader* evtHeader = findNode::getClass<EventHeader>(topNode, "EventHeader");

  int m_runNumber;
  int m_evtNumber;
  if (evtHeader)
  {
    m_runNumber = evtHeader->get_RunNumber();
    m_evtNumber = evtHeader->get_EvtSequence();
  }
  else
  {
    m_runNumber = m_evtNumber = -1;
  }

  // Loop over tracks and check for close DCA match with all other tracks
  for (auto tr1_it = m_svtxTrackMap->begin(); tr1_it != m_svtxTrackMap->end(); ++tr1_it)
  {
    auto id1 = tr1_it->first;
    auto *tr1 = tr1_it->second;
    if (tr1->get_quality() > _qual_cut)
    {
      continue;
    }
    if (tr1->get_pt() < track_pt_cut)
    {
      continue;
    }

    short int crossing1 = tr1->get_crossing();

    // calculate number silicon tracks
    double this_dca_cut = track_dca_cut;
    TrackSeed* siliconseed = tr1->get_silicon_seed();

    if (!siliconseed)
    {
      this_dca_cut *= 5;
      if (Verbosity() > 2)
      {
        std::cout << "silicon seed not found" << std::endl;
      }
      if (_require_mvtx)
      {
        continue;
      }
    }

    std::vector<unsigned int> nstates1 = getTrackStates(tr1);
    unsigned int track1_mvtx_state_size = nstates1[0];
    unsigned int track1_intt_state_size = nstates1[1];
    // unsigned int track1_tpc_state_size = nstates1[2];
    // unsigned int track1_mms_state_size = nstates1[3];
    
    unsigned int track1_silicon_cluster_size = std::numeric_limits<unsigned int>::quiet_NaN();
    if (siliconseed)
    {
      track1_silicon_cluster_size = siliconseed->size_cluster_keys();
    }

    std::vector<TrkrDefs::cluskey> ckeys1;
    if (siliconseed)
    {
      ckeys1.insert(ckeys1.end(), siliconseed->begin_cluster_keys(), siliconseed->end_cluster_keys());
    }

    unsigned int track1_mvtx_cluster_size = 0;
    unsigned int track1_intt_cluster_size = 0;
    for (const auto& ckey : ckeys1)
    {
      auto detid = TrkrDefs::getTrkrId(ckey);
      if (detid == TrkrDefs::TrkrId::mvtxId)
      {
        track1_mvtx_cluster_size++;
      }
      else if (detid == TrkrDefs::TrkrId::inttId)
      {
        track1_intt_cluster_size++;
      }
    }

    Acts::Vector3 pos1(tr1->get_x(), tr1->get_y(), tr1->get_z());
    Acts::Vector3 mom1(tr1->get_px(), tr1->get_py(), tr1->get_pz());
    Acts::Vector3 dcaVals1 = calculateDca(tr1, mom1, pos1);
    // first dca cuts
    if (fabs(dcaVals1(0)) < this_dca_cut or fabs(dcaVals1(1)) < this_dca_cut)
    {
      continue;
    }

    // look for close DCA matches with all other such tracks
    for (auto tr2_it = std::next(tr1_it); tr2_it != m_svtxTrackMap->end(); ++tr2_it)
    {
      auto id2 = tr2_it->first;
      auto *tr2 = tr2_it->second;
      if (tr2->get_quality() > _qual_cut)
      {
        continue;
      }
      if (tr2->get_pt() < track_pt_cut)
      {
        continue;
      }

      short int crossing2 = tr2->get_crossing();

      // calculate number silicon tracks
      TrackSeed* siliconseed2 = tr2->get_silicon_seed();

      double this_dca_cut2 = track_dca_cut;

      if (!siliconseed2)
      {
        this_dca_cut2 *= 5;
        if (Verbosity() > 2)
        {
          std::cout << "silicon seed not found" << std::endl;
        }
        if (_require_mvtx)
        {
          continue;
        }
      }

      std::vector<unsigned int> nstates2 = getTrackStates(tr2);
      unsigned int track2_mvtx_state_size = nstates2[0];
      unsigned int track2_intt_state_size = nstates2[1];
      // unsigned int track2_tpc_state_size = nstates2[2];
      // unsigned int track2_mms_state_size = nstates2[3];
      
      unsigned int track2_silicon_cluster_size = std::numeric_limits<unsigned int>::quiet_NaN();
      if (siliconseed2)
      {
        track2_silicon_cluster_size = siliconseed2->size_cluster_keys();
      }

      std::vector<TrkrDefs::cluskey> ckeys2;
      if (siliconseed2)
      {
        ckeys2.insert(ckeys2.end(), siliconseed2->begin_cluster_keys(), siliconseed2->end_cluster_keys());
      }

      unsigned int track2_mvtx_cluster_size = 0;
      unsigned int track2_intt_cluster_size = 0;
      for (const auto& ckey : ckeys2)
      {
         auto detid = TrkrDefs::getTrkrId(ckey);
        if (detid == TrkrDefs::TrkrId::mvtxId)
        {
          track2_mvtx_cluster_size++;
        }
        else if (detid == TrkrDefs::TrkrId::inttId)
        {
          track2_intt_cluster_size++;
        }
      }

      // dca xy and dca z cut here compare to track dca cut
      Acts::Vector3 pos2(tr2->get_x(), tr2->get_y(), tr2->get_z());
      Acts::Vector3 mom2(tr2->get_px(), tr2->get_py(), tr2->get_pz());
      Acts::Vector3 dcaVals2 = calculateDca(tr2, mom2, pos2);

      if (fabs(dcaVals2(0)) < this_dca_cut2 or fabs(dcaVals2(1)) < this_dca_cut2)
      {
        continue;
      }

      // find DCA of these two tracks
      if (Verbosity() > 3)
      {
        std::cout << "Check DCA for tracks " << id1 << " and  " << id2 << std::endl;
      }

      if (tr1->get_charge() == tr2->get_charge())
      {
        // continue;
      }

      // declare these variables to pass into findPCAtwoTracks and fillHistogram by reference
      double pair_dca;
      Acts::Vector3 pca_rel1;
      Acts::Vector3 pca_rel2;
      double invariantMass;
      double invariantPt;
      float invariantPhi;
      float rapidity;
      float pseudorapidity;

      // Initial calculation of point of closest approach between the two tracks
      // This presently assumes straight line tracks to get a rough answer
      // Should update to use circles instead?
      findPcaTwoTracks(pos1, pos2, mom1, mom2, pca_rel1, pca_rel2, pair_dca);

      // tracks with small relative pca are k short candidates
      if (abs(pair_dca) < pair_dca_cut)
      {
        // Pair pca and dca were calculated with nominal track parameters and are approximate
        // Project tracks to this rough pca
        Eigen::Vector3d projected_pos1;
        Eigen::Vector3d projected_mom1;
        Eigen::Vector3d projected_pos2;
        Eigen::Vector3d projected_mom2;

        bool ret1 = projectTrackToPoint(tr1, pca_rel1, projected_pos1, projected_mom1);
        bool ret2 = projectTrackToPoint(tr2, pca_rel2, projected_pos2, projected_mom2);

        if (!ret1 or !ret2)
        {
          continue;
        }

        // recalculate pca starting with projected position and momentum
        double pair_dca_proj;
        Acts::Vector3 pca_rel1_proj;
        Acts::Vector3 pca_rel2_proj;
        findPcaTwoTracks(projected_pos1, projected_pos2, projected_mom1, projected_mom2, pca_rel1_proj, pca_rel2_proj, pair_dca_proj);

        // if(pair_dca_proj > pair_dca_cut) continue;

        // invariant mass is calculated in this method
        fillHistogram(projected_mom1, projected_mom2, recomass, invariantMass, invariantPt, invariantPhi, rapidity, pseudorapidity);
        fillNtp(tr1, tr2, dcaVals1, dcaVals2, pca_rel1, pca_rel2, pair_dca, invariantMass, invariantPt, invariantPhi, rapidity, pseudorapidity, projected_pos1, projected_pos2, projected_mom1, projected_mom2, pca_rel1_proj, pca_rel2_proj, pair_dca_proj, track1_silicon_cluster_size, track2_silicon_cluster_size, track1_mvtx_cluster_size, track1_mvtx_state_size, track1_intt_cluster_size, track1_intt_state_size, track2_mvtx_cluster_size, track2_mvtx_state_size, track2_intt_cluster_size, track2_intt_state_size, m_runNumber, m_evtNumber);

        if (Verbosity() > 1)
        {
          std::cout << " Accepted Track Pair" << std::endl;
	  std::cout << " id1 " << id1 << " id2 " << id2 << std::endl;
	  std::cout << " crossing1 " << crossing1 << " crossing2 " << crossing2 << std::endl;
          std::cout << " invariant mass: " << invariantMass << std::endl;
          std::cout << " track1 dca_cut: " << this_dca_cut << " track2 dca_cut: " << this_dca_cut2 << std::endl;
          std::cout << " dca3dxy1,dca3dz1,phi1: " << dcaVals1 << std::endl;
          std::cout << " dca3dxy2,dca3dz2,phi2: " << dcaVals2 << std::endl;
          std::cout << "Initial:  pca_rel1: " << pca_rel1 << " pca_rel2: " << pca_rel2 << std::endl;
          std::cout << " Initial: mom1: " << mom1 << " mom2: " << mom2 << std::endl;
          std::cout << "Proj_pca_rel:  proj_pos1: " << projected_pos1 << " proj_pos2: " << projected_pos2 << " proj_mom1: " << projected_mom1 << " proj_mom2: " << projected_mom2 << std::endl;
          std::cout << " Relative PCA = " << abs(pair_dca) << " pca_cut = " << pair_dca_cut << std::endl;
          std::cout << " charge 1: " << tr1->get_charge() << " charge2: " << tr2->get_charge() << std::endl;
          std::cout << "found viable projection" << std::endl;
          std::cout << "Final: pca_rel1_proj: " << pca_rel1_proj << " pca_rel2_proj: " << pca_rel2_proj << " mom1: " << projected_mom1 << " mom2: " << projected_mom2 << std::endl
                    << std::endl;
        }

        if (m_save_tracks)
        {
          m_output_trackMap = findNode::getClass<SvtxTrackMap>(topNode, m_output_trackMap_node_name);
          m_output_trackMap->insertWithKey(tr1, tr1->get_id());
          m_output_trackMap->insertWithKey(tr2, tr2->get_id());
        }

      }
    }
  }
  return 0;
}

std::vector<unsigned int> KshortReconstruction::getTrackStates(SvtxTrack *track)
{
  std::vector<unsigned int> nstates;
  unsigned int nmapsstate = 0;
  unsigned int ninttstate = 0;
  unsigned int ntpcstate = 0;
  unsigned int nmmsstate = 0;
  
  // the track states from the Acts fit are fitted to fully corrected clusters, and are on the surface
  for (auto state_iter = track->begin_states();
       state_iter != track->end_states();
       ++state_iter)
    {
      SvtxTrackState* tstate = state_iter->second;
      auto stateckey = tstate->get_cluskey();

      switch (TrkrDefs::getTrkrId(stateckey))
	{
	case TrkrDefs::mvtxId:
	  nmapsstate++;
	  break;
	case TrkrDefs::inttId:
	  ninttstate++;
	  break;
	case TrkrDefs::tpcId:
	  ntpcstate++;
	  break;
	case TrkrDefs::micromegasId:
	  nmmsstate++;
	  break;
	default:
	  std::cout << PHWHERE << " unknown key " << stateckey << std::endl;
	  gSystem->Exit(1);
	  exit(1);
	}
    }
  nstates.push_back(nmapsstate);
  nstates.push_back(ninttstate);
  nstates.push_back(ntpcstate);
  nstates.push_back(nmmsstate);

  return nstates;
}

void KshortReconstruction::fillNtp(SvtxTrack* track1, SvtxTrack* track2, Acts::Vector3 dcavals1, Acts::Vector3 dcavals2, Acts::Vector3 pca_rel1, Acts::Vector3 pca_rel2, double pair_dca, double invariantMass, double invariantPt, float invariantPhi, float rapidity, float pseudorapidity, Eigen::Vector3d projected_pos1, Eigen::Vector3d projected_pos2, Eigen::Vector3d projected_mom1, Eigen::Vector3d projected_mom2, Acts::Vector3 pca_rel1_proj, Acts::Vector3 pca_rel2_proj, double pair_dca_proj, unsigned int track1_silicon_cluster_size, unsigned int track2_silicon_cluster_size, unsigned int track1_mvtx_cluster_size,  unsigned int track1_mvtx_state_size, unsigned int track1_intt_cluster_size,  unsigned int track1_intt_state_size, unsigned int track2_mvtx_cluster_size,  unsigned int track2_mvtx_state_size, unsigned int track2_intt_cluster_size,  unsigned int track2_intt_state_size, int runNumber, int eventNumber)
{
  double px1 = track1->get_px();
  double py1 = track1->get_py();
  double pz1 = track1->get_pz();
  auto *tpcSeed1 = track1->get_tpc_seed();
  size_t tpcClusters1 = tpcSeed1->size_cluster_keys();
  double eta1 = asinh(pz1 / sqrt(pow(px1, 2) + pow(py1, 2)));

  double px2 = track2->get_px();
  double py2 = track2->get_py();
  double pz2 = track2->get_pz();
  auto *tpcSeed2 = track2->get_tpc_seed();
  size_t tpcClusters2 = tpcSeed2->size_cluster_keys();
  double eta2 = asinh(pz2 / sqrt(pow(px2, 2) + pow(py2, 2)));

  auto vtxid = track1->get_vertex_id();

  Acts::Vector3 vertex(0, 0, track1->get_z());  // fake primary vertex
  auto *svtxVertex = m_vertexMap->get(vtxid);
  if (svtxVertex)
  {
    vertex(0) = svtxVertex->get_x();
    vertex(1) =  svtxVertex->get_y();
    vertex(2) = svtxVertex->get_z(); 
  }

  Acts::Vector3 pathLength = (pca_rel1 + pca_rel2) * 0.5 - vertex;
  Acts::Vector3 pathLength_proj = (pca_rel1_proj + pca_rel2_proj) * 0.5 - vertex;

  float mag_pathLength = sqrt(pow(pathLength(0), 2) + pow(pathLength(1), 2) + pow(pathLength(2), 2));
  float mag_pathLength_proj = sqrt(pow(pathLength_proj(0), 2) + pow(pathLength_proj(1), 2) + pow(pathLength_proj(2), 2));

  Acts::Vector3 projected_momentum = projected_mom1 + projected_mom2;
  float cos_theta_reco = pathLength_proj.dot(projected_momentum) / (projected_momentum.norm() * pathLength_proj.norm());


  float reco_info[] = {(float) track1->get_id(), (float) track1->get_crossing(), track1->get_x(), track1->get_y(), track1->get_z(), track1->get_px(), track1->get_py(), track1->get_pz(), (float) dcavals1(0), (float) dcavals1(1), (float) dcavals1(2), (float) pca_rel1(0), (float) pca_rel1(1), (float) pca_rel1(2), (float) eta1, (float) track1->get_charge(), (float) tpcClusters1, (float) track2->get_id(), (float) track2->get_crossing(), track2->get_x(), track2->get_y(), track2->get_z(), track2->get_px(), track2->get_py(), track2->get_pz(), (float) dcavals2(0), (float) dcavals2(1), (float) dcavals2(2), (float) pca_rel2(0), (float) pca_rel2(1), (float) pca_rel2(2), (float) eta2, (float) track2->get_charge(), (float) tpcClusters2, (float) vertex(0), (float) vertex(1), (float) vertex(2), (float) pair_dca, (float) invariantMass, (float) invariantPt, invariantPhi, (float) pathLength(0), (float) pathLength(1), (float) pathLength(2), mag_pathLength, rapidity, pseudorapidity, (float) projected_pos1(0), (float) projected_pos1(1), (float) projected_pos1(2), (float) projected_pos2(0), (float) projected_pos2(1), (float) projected_pos2(2), (float) projected_mom1(0), (float) projected_mom1(1), (float) projected_mom1(2), (float) projected_mom2(0), (float) projected_mom2(1), (float) projected_mom2(2), (float) pca_rel1_proj(0), (float) pca_rel1_proj(1), (float) pca_rel1_proj(2), (float) pca_rel2_proj(0), (float) pca_rel2_proj(1), (float) pca_rel2_proj(2), (float) pair_dca_proj, (float) pathLength_proj(0), (float) pathLength_proj(1), (float) pathLength_proj(2), mag_pathLength_proj, track1->get_quality(), track2->get_quality(), cos_theta_reco, (float) track1_silicon_cluster_size, (float) track2_silicon_cluster_size, (float) track1_mvtx_cluster_size, (float) track1_mvtx_state_size, (float) track1_intt_cluster_size,  (float) track1_intt_state_size, (float) track2_mvtx_cluster_size, (float) track2_mvtx_state_size,  (float) track2_intt_cluster_size, (float) track2_intt_state_size, (float) runNumber, (float) eventNumber};

  ntp_reco_info->Fill(reco_info);
}

void KshortReconstruction::fillHistogram(Eigen::Vector3d mom1, Eigen::Vector3d mom2, TH1* massreco, double& invariantMass, double& invariantPt, float& invariantPhi, float& rapidity, float& pseudorapidity)
{
  double E1 = sqrt(pow(mom1(0), 2) + pow(mom1(1), 2) + pow(mom1(2), 2) + pow(decaymass, 2));
  double E2 = sqrt(pow(mom2(0), 2) + pow(mom2(1), 2) + pow(mom2(2), 2) + pow(decaymass, 2));

  TLorentzVector v1(mom1(0), mom1(1), mom1(2), E1);
  TLorentzVector v2(mom2(0), mom2(1), mom2(2), E2);

  TLorentzVector tsum;
  tsum = v1 + v2;

  rapidity = tsum.Rapidity();
  pseudorapidity = tsum.Eta();
  invariantMass = tsum.M();
  invariantPt = tsum.Pt();
  invariantPhi = tsum.Phi();

  if (Verbosity() > 2)
  {
    std::cout << "px1: " << mom1(0) << " py1: " << mom1(1) << " pz1: " << mom1(2) << " E1: " << E1 << std::endl;
    std::cout << "px2: " << mom2(0) << " py2: " << mom2(1) << " pz2: " << mom2(2) << " E2: " << E2 << std::endl;
    std::cout << "tsum: " << tsum(0) << " " << tsum(1) << " " << tsum(2) << " " << tsum(3) << std::endl;
    std::cout << "invariant mass: " << invariantMass << " invariant Pt: " << invariantPt << " invariantPhi: " << invariantPhi << std::endl;
  }

  if (invariantPt > invariant_pt_cut)
  {
    massreco->Fill(invariantMass);
  }
}

bool KshortReconstruction::projectTrackToPoint(SvtxTrack* track, Eigen::Vector3d PCA, Eigen::Vector3d& pos, Eigen::Vector3d& mom)
{
  bool ret = true;

  /// create perigee surface
  ActsPropagator actsPropagator(_tGeometry);
  auto perigee = actsPropagator.makeVertexSurface(PCA);  // PCA is in cm here
  auto params = actsPropagator.makeTrackParams(track, m_vertexMap);
  if (!params.ok())
  {
    return false;
  }
  auto result = actsPropagator.propagateTrack(params.value(), perigee);

  if (result.ok())
  {
    auto projectionPos = result.value().second.position(_tGeometry->geometry().getGeoContext());
    const auto momentum = result.value().second.momentum();
    pos(0) = projectionPos.x() / Acts::UnitConstants::cm;
    pos(1) = projectionPos.y() / Acts::UnitConstants::cm;
    pos(2) = projectionPos.z() / Acts::UnitConstants::cm;

    if (Verbosity() > 2)
    {
      std::cout << "                 Input PCA " << PCA << "  projection out " << pos << std::endl;
    }

    mom(0) = momentum.x();
    mom(1) = momentum.y();
    mom(2) = momentum.z();
  }
  else
  {
    pos(0) = track->get_x();
    pos(1) = track->get_y();
    pos(2) = track->get_z();

    mom(0) = track->get_px();
    mom(1) = track->get_py();
    mom(2) = track->get_pz();

    if(Verbosity() > 0)
      {
	std::cout << result.error() << std::endl;
	std::cout << result.error().message() << std::endl;
	std::cout << " Failed projection of track with: " << std::endl;
	std::cout << " x,y,z = " << track->get_x() << "  " << track->get_y() << "  " << track->get_z() << std::endl;
	std::cout << " px,py,pz = " << track->get_px() << "  " << track->get_py() << "  " << track->get_pz() << std::endl;
	std::cout << " to point (x,y,z) = " << PCA(0) / Acts::UnitConstants::cm << "  " << PCA(1) / Acts::UnitConstants::cm << "  " << PCA(2) / Acts::UnitConstants::cm << std::endl;
      }

    //    ret = false;
  }

  return ret;
}

bool KshortReconstruction::projectTrackToCylinder(SvtxTrack* track, double Radius, Eigen::Vector3d& pos, Eigen::Vector3d& mom)
{
  // Make a cylinder surface at the radius and project the track to that
  bool ret = true;
  const double eta = 2.0;
  const double theta = 2. * atan(exp(-eta));
  const double halfZ = Radius / tan(theta) * Acts::UnitConstants::cm;
  Radius *= Acts::UnitConstants::cm;

  /// Make a cylindrical surface at (0,0,0) aligned along the z axis
  auto transform = Acts::Transform3::Identity();

  std::shared_ptr<Acts::CylinderSurface> cylSurf =
      Acts::Surface::makeShared<Acts::CylinderSurface>(transform,
                                                       Radius,
                                                       halfZ);
  ActsPropagator actsPropagator(_tGeometry);
  auto params = actsPropagator.makeTrackParams(track, m_vertexMap);
  if (!params.ok())
  {
    return false;
  }

  auto result = actsPropagator.propagateTrack(params.value(), cylSurf);
  if (result.ok())
  {
    auto projectionPos = result.value().second.position(_tGeometry->geometry().getGeoContext());
    const auto momentum = result.value().second.momentum();
    pos(0) = projectionPos.x() / Acts::UnitConstants::cm;
    pos(1) = projectionPos.y() / Acts::UnitConstants::cm;
    pos(2) = projectionPos.z() / Acts::UnitConstants::cm;

    mom(0) = momentum.x();
    mom(1) = momentum.y();
    mom(2) = momentum.z();
  }
  else
  {
    ret = false;
  }

  return ret;
}

Acts::Vector3 KshortReconstruction::getVertex(SvtxTrack* track)
{
  auto vertexId = track->get_vertex_id();
  const SvtxVertex* svtxVertex = m_vertexMap->get(vertexId);
  Acts::Vector3 vertex = Acts::Vector3::Zero();
  if (svtxVertex)
  {
    vertex(0) = svtxVertex->get_x() * Acts::UnitConstants::cm;
    vertex(1) = svtxVertex->get_y() * Acts::UnitConstants::cm;
    vertex(2) = svtxVertex->get_z() * Acts::UnitConstants::cm;
  }

  return vertex;
}

void KshortReconstruction::findPcaTwoTracks(const Acts::Vector3& pos1, const Acts::Vector3& pos2, Acts::Vector3 mom1, Acts::Vector3 mom2, Acts::Vector3& pca1, Acts::Vector3& pca2, double& dca) const
{
  TLorentzVector v1;
  TLorentzVector v2;

  double px1 = mom1(0);
  double py1 = mom1(1);
  double pz1 = mom1(2);
  double px2 = mom2(0);
  double py2 = mom2(1);
  double pz2 = mom2(2);

  Float_t E1 = sqrt(pow(px1, 2) + pow(py1, 2) + pow(pz1, 2) + pow(decaymass, 2));
  Float_t E2 = sqrt(pow(px2, 2) + pow(py2, 2) + pow(pz2, 2) + pow(decaymass, 2));

  v1.SetPxPyPzE(px1, py1, pz1, E1);
  v2.SetPxPyPzE(px2, py2, pz2, E2);

  // calculate lorentz vector
  const Eigen::Vector3d& a1 = pos1;
  const Eigen::Vector3d& a2 = pos2;

  Eigen::Vector3d b1(v1.Px(), v1.Py(), v1.Pz());
  Eigen::Vector3d b2(v2.Px(), v2.Py(), v2.Pz());

  // The shortest distance between two skew lines described by
  //  a1 + c * b1
  //  a2 + d * b2
  // where a1, a2, are vectors representing points on the lines, b1, b2 are direction vectors, and c and d are scalars
  // dca = (b1 x b2) .(a2-a1) / |b1 x b2|

  // bcrossb/mag_bcrossb is a unit vector perpendicular to both direction vectors b1 and b2
  auto bcrossb = b1.cross(b2);
  auto mag_bcrossb = bcrossb.norm();
  // a2-a1 is the vector joining any arbitrary points on the two lines
  auto aminusa = a2 - a1;

  // The DCA of these two lines is the projection of a2-a1 along the direction of the perpendicular to both
  // remember that a2-a1 is longer than (or equal to) the dca by definition
  dca = 999;
  if (mag_bcrossb != 0)
  {
    dca = bcrossb.dot(aminusa) / mag_bcrossb;
  }
  else
  {
    return;  // same track, skip combination
  }

  // get the points at which the normal to the lines intersect the lines, where the lines are perpendicular
  double X = b1.dot(b2) - (b1.dot(b1) * b2.dot(b2) / b2.dot(b1));
  double Y = (a2.dot(b2) - a1.dot(b2)) - ((a2.dot(b1) - a1.dot(b1)) * b2.dot(b2) / b2.dot(b1));
  double c = Y / X;

  double F = b1.dot(b1) / b2.dot(b1);
  double G = -(a2.dot(b1) - a1.dot(b1)) / b2.dot(b1);
  double d = (c * F) + G;

  // then the points of closest approach are:
  pca1 = a1 + c * b1;
  pca2 = a2 + d * b2;

  return;
}

KshortReconstruction::KshortReconstruction(const std::string& name)
  : SubsysReco(name)
{
}

Acts::Vector3 KshortReconstruction::calculateDca(SvtxTrack* track, const Acts::Vector3& momentum, Acts::Vector3 position)
{
  // For the purposes of this module, we set default values to prevent this track from being rejected if the dca calc fails
  Acts::Vector3 r = momentum.cross(Acts::Vector3(0., 0., 1.));
  float phi = atan2(r(1), r(0));
  Acts::Vector3 outVals(track_dca_cut*1.1, track_dca_cut*1.1, phi);
  auto vtxid = track->get_vertex_id();
  if (!m_vertexMap)
  {
    //std::cout << "Could not find m_vertexmap " << std::endl;
    return outVals;
  }
  auto *svtxVertex = m_vertexMap->get(vtxid);
  if (!svtxVertex)
  {
    //std::cout << "Could not find vtxid in m_vertexMap " << vtxid << std::endl;
    return outVals;
  }
  Acts::Vector3 vertex(svtxVertex->get_x(), svtxVertex->get_y(), svtxVertex->get_z());
  position -= vertex;

  Acts::RotationMatrix3 rot;
  rot(0, 0) = std::cos(phi);
  rot(0, 1) = -std::sin(phi);
  rot(0, 2) = 0;
  rot(1, 0) = std::sin(phi);
  rot(1, 1) = std::cos(phi);
  rot(1, 2) = 0;
  rot(2, 0) = 0;
  rot(2, 1) = 0;
  rot(2, 2) = 1;

  Acts::Vector3 pos_R = rot * position;
  double dca3dxy = pos_R(0);
  double dca3dz = pos_R(2);

  outVals(0) = abs(dca3dxy);
  outVals(1) = abs(dca3dz);
  outVals(2) = phi;

  if (Verbosity() > 4)
  {
    std::cout << " pre-position: " << position << std::endl;
    std::cout << " vertex: " << vertex << std::endl;
    std::cout << " vertex subtracted-position: " << position << std::endl;
  }

  return outVals;
}

int KshortReconstruction::InitRun(PHCompositeNode* topNode)
{
  const char* cfilepath = filepath.c_str();
  fout = new TFile(cfilepath, "recreate");
  ntp_reco_info = new TNtuple("ntp_reco_info", "decay_pairs", "id1:crossing1:x1:y1:z1:px1:py1:pz1:dca3dxy1:dca3dz1:phi1:pca_rel1_x:pca_rel1_y:pca_rel1_z:eta1:charge1:tpcClusters_1:id2:crossing2:x2:y2:z2:px2:py2:pz2:dca3dxy2:dca3dz2:phi2:pca_rel2_x:pca_rel2_y:pca_rel2_z:eta2:charge2:tpcClusters_2:vertex_x:vertex_y:vertex_z:pair_dca:invariant_mass:invariant_pt:invariantPhi:pathlength_x:pathlength_y:pathlength_z:pathlength:rapidity:pseudorapidity:projected_pos1_x:projected_pos1_y:projected_pos1_z:projected_pos2_x:projected_pos2_y:projected_pos2_z:projected_mom1_x:projected_mom1_y:projected_mom1_z:projected_mom2_x:projected_mom2_y:projected_mom2_z:projected_pca_rel1_x:projected_pca_rel1_y:projected_pca_rel1_z:projected_pca_rel2_x:projected_pca_rel2_y:projected_pca_rel2_z:projected_pair_dca:projected_pathlength_x:projected_pathlength_y:projected_pathlength_z:projected_pathlength:quality1:quality2:cosThetaReco:track1_silicon_clusters:track2_silicon_clusters:track1_mvtx_clusters:track1_mvtx_states:track1_intt_clusters:track1_intt_states:track2_mvtx_clusters:track2_mvtx_states:track2_intt_clusters:track2_intt_states:runNumber:eventNumber");

  getNodes(topNode);

  recomass = new TH1D("recomass", "recomass", 1000, 0.0, 1);  // root histogram arguments: name,title,bins,minvalx,maxvalx

  //Add new track map to save selected tracks
  if (m_save_tracks)
  {
    PHNodeIterator nodeIter(topNode);

    PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode *>(nodeIter.findFirst("PHCompositeNode", "DST"));
    if (!dstNode)
    {
      dstNode = new PHCompositeNode("DST");
      topNode->addNode(dstNode);
      std::cout << "DST node added" << std::endl;
    }

    m_output_trackMap = new SvtxTrackMap_v2();
    PHIODataNode<PHObject> *outputTrackNode = new PHIODataNode<PHObject>(m_output_trackMap, m_output_trackMap_node_name, "PHObject");
    dstNode->addNode(outputTrackNode);
    if (Verbosity() > 1) { std::cout << m_output_trackMap_node_name << " node added" << std::endl; }
  }

  return 0;
}

int KshortReconstruction::End(PHCompositeNode* /**topNode*/)
{
  fout->cd();
  ntp_reco_info->Write();
  recomass->Write();
  fout->Close();

  return 0;
}

int KshortReconstruction::getNodes(PHCompositeNode* topNode)
{
  m_svtxTrackMap = findNode::getClass<SvtxTrackMap>(topNode, "SvtxTrackMap");
  if (!m_svtxTrackMap)
  {
    std::cout << PHWHERE << "No SvtxTrackMap on node tree, exiting." << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  m_vertexMap = findNode::getClass<SvtxVertexMap>(topNode, "SvtxVertexMap");
  if (!m_vertexMap)
  {
    std::cout << PHWHERE << "No vertexMap on node tree, exiting." << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  _tGeometry = findNode::getClass<ActsGeometry>(topNode, "ActsGeometry");
  if (!_tGeometry)
  {
    std::cout << PHWHERE << "Error, can't find acts tracking geometry" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}
