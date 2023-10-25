#include "SecondaryVertexFinder.h"
#include "ActsPropagator.h"

/// Tracking includes

#include <trackbase/TrkrCluster.h>            // for TrkrCluster
#include <trackbase/TrkrDefs.h>               // for cluskey, getLayer, TrkrId
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrackFitUtils.h>
#include <trackbase/ActsGeometry.h>
#include <trackbase/TpcDefs.h>

#include <trackbase_historic/SvtxTrack.h>     // for SvtxTrack, SvtxTrack::C...
#include <trackbase_historic/SvtxTrackMap_v2.h>
#include <trackbase_historic/SvtxVertexMap.h>
#include <trackbase_historic/SvtxVertex_v1.h>
#include <trackbase_historic/SvtxVertexMap_v1.h>
#include <trackbase_historic/TrackSeed_v1.h>
#include <trackbase_historic/ActsTransformations.h>

#include <Acts/Geometry/GeometryIdentifier.hpp>
#include <Acts/Surfaces/PerigeeSurface.hpp>
#include <Acts/MagneticField/MagneticFieldProvider.hpp>

#include <phool/PHCompositeNode.h>
#include <phool/PHNode.h>
#include <phool/PHNodeIterator.h>
#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/getClass.h>
#include <phool/phool.h>

#include <cmath>                              // for sqrt, fabs, atan2, cos
#include <iostream>                           // for operator<<, basic_ostream
#include <map>                                // for map
#include <set>                                // for _Rb_tree_const_iterator
#include <utility>                            // for pair, make_pair
#include <iomanip>

#include <vector>
#include <cassert>
#include <numeric>
#include <iostream>
#include <algorithm>
#include <functional>

#include <Eigen/Dense>

#include <TLorentzVector.h>
#include <TNtuple.h>
#include <TH2D.h>
#include <TH1D.h>
#include <TFile.h>

//____________________________________________________________________________..
SecondaryVertexFinder::SecondaryVertexFinder(const std::string &name)
  : SubsysReco(name)
{

}

//____________________________________________________________________________..
SecondaryVertexFinder::~SecondaryVertexFinder()
{

}

//____________________________________________________________________________..
int SecondaryVertexFinder::InitRun(PHCompositeNode *topNode)
{
  if(_write_electrons_node)
    {
      CreateOutputNode(topNode);
    }

  int ret = GetNodes(topNode);
  if (ret != Fun4AllReturnCodes::EVENT_OK) return ret;

  if(_write_ntuple)
    {
      recomass = new TH2D("recomass", "invariant mass vs pT", 1000, 0, 5, 4000,0,2);
      hdecaypos = new TH2D("hdecaypos","decay radius",200, -40, 40, 200,-40,40);
      hdecaypos->GetXaxis()->SetTitle("decay X (cm)");
      hdecaypos->GetYaxis()->SetTitle("decay Y (cm)");
      hdecay_radius = new TH1D("hdecay_radius", "Decay Radius", 200, 0, 40.0);
      
      ntp = new TNtuple("ntp","decay_pairs","x1:y1:z1:px1:py1:pz1:dca3dxy1:dca3dz1:vposx1:vposy1:vposz1:vmomx1:vmomy1:vmomz1:pca_relx_1:pca_rely_1:pca_relz_1:eta1:charge1:tpcClusters_1:quality1:eta1:x2:y2:z2:px2:py2:pz2:dca3dxy2:dca3dz2:vposx2:vposy2:vposz2:vmomx2:vmomy2:vmomz2:pca_relx_2:pca_rely_2:pca_relz_2:eta2:charge2:tpcClusters_2:quality2:eta2:vertex_x:vertex_y:vertex_z:pair_dca:invariant_mass:invariant_pt:path:has_silicon1:has_silicon2");
    }
  
  return ret;
}

void SecondaryVertexFinder::fillNtp(SvtxTrack *track1, SvtxTrack *track2, double dca3dxy1, double dca3dz1, double dca3dxy2, double dca3dz2,  Eigen::Vector3d vpos1,  Eigen::Vector3d vmom1, Eigen::Vector3d vpos2, Eigen::Vector3d vmom2, Acts::Vector3 pca_rel1, Acts::Vector3 pca_rel2, double pair_dca, double invariantMass, double invariantPt, double path, int has_silicon_1, int has_silicon_2)
{
  double px1          = track1->get_px();
  double py1          = track1->get_py();
  double pz1          = track1->get_pz();
  auto tpcSeed1       = track1->get_tpc_seed();
  size_t tpcClusters1 = tpcSeed1->size_cluster_keys();
  double eta1         = asinh(pz1 / sqrt(pow(px1, 2) + pow(py1, 2)));

  double px2          = track2->get_px();
  double py2          = track2->get_py();
  double pz2          = track2->get_pz();
  auto tpcSeed2       = track2->get_tpc_seed();
  size_t tpcClusters2 = tpcSeed2->size_cluster_keys();
  double eta2         = asinh(pz2 / sqrt(pow(px2, 2) + pow(py2, 2)));

  auto vtxid      = track1->get_vertex_id();
  auto svtxVertex = _svtx_vertex_map->get(vtxid);
  if(!svtxVertex){ return; }
  
  float reco_info[] = {
    track1->get_x(), track1->get_y(), track1->get_z(), 
    track1->get_px(), track1->get_py(), track1->get_pz(), 
    (float) dca3dxy1, (float) dca3dz1, (float) vpos1(0), (float) vpos1(1), (float) vpos1(2),
    (float) vmom1(0), (float) vmom1(1), (float) vmom1(2),
    (float) pca_rel1(0), (float) pca_rel1(1), (float) pca_rel1(2), 
    (float) eta1,  (float) track1->get_charge(), (float) tpcClusters1, 
    (float) track1->get_quality(), (float) eta1,
    track2->get_x(), track2->get_y(), track2->get_z(),  
    track2->get_px(), track2->get_py(), track2->get_pz(), 
    (float) dca3dxy2, (float) dca3dz2, (float) vpos2(0), (float) vpos2(1), (float) vpos2(2),
    (float) vmom2(0), (float) vmom2(1), (float) vmom2(2),
    (float) pca_rel2(0), (float) pca_rel2(1), (float) pca_rel2(2), 
    (float) eta2, (float) track2->get_charge(), (float) tpcClusters2, 
    (float) track2->get_quality(), (float) eta2,
    svtxVertex->get_x(), svtxVertex->get_y(), svtxVertex->get_z(), 
    (float) pair_dca,(float) invariantMass, (float) invariantPt, (float) path,
    (float) has_silicon_1, (float)  has_silicon_2};


  ntp->Fill(reco_info);
}

//____________________________________________________________________________..
int SecondaryVertexFinder::process_event(PHCompositeNode */*topNode*/)
{
  if(Verbosity() > 0)
    std::cout << PHWHERE << " track map size " << _track_map->size()  << std::endl;

  // Loop over tracks and check for close DCA match with all other tracks
  for(auto tr1_it = _track_map->begin(); tr1_it != _track_map->end(); ++tr1_it)
    {
      auto id1 = tr1_it->first;
      auto tr1 = tr1_it->second;

      auto vertexId = tr1->get_vertex_id();
      const SvtxVertex* svtxVertex = _svtx_vertex_map->get(vertexId);
      if(!svtxVertex)
	{
	  if(Verbosity() > 1)
	    { std::cout << PHWHERE << " Failed to find vertex id " << vertexId << " skip this track " << std::endl; }
	  continue;
	}
      if(svtxVertex->size_tracks() == 0) continue;  // event did not have a reconstructed vertex, vertex is bogus 

      // Reverse or remove this to consider TPC-only tracks too
      if(_require_mvtx && !hasSiliconSeed(tr1)) continue;

      int has_silicon_1 = 0;
      if(hasSiliconSeed(tr1)) has_silicon_1 = 1;


      if(Verbosity() > 3)
	{
	  std::cout << "Track1 " << id1 << " details: " << std::endl;
	  outputTrackDetails(tr1);
	}
      
      if(tr1->get_quality() > _qual_cut) continue;

      auto tpc_seed1 = tr1->get_tpc_seed();
      int ntpc1 = tpc_seed1->size_cluster_keys();
      if(ntpc1 < 20) continue;

      float dca3dxy1, dca3dz1,dca3dxysigma1, dca3dzsigma1;
      get_dca(tr1, dca3dxy1, dca3dz1, dca3dxysigma1, dca3dzsigma1);
      //std::cout << "   tr1 " << id1 << " dca3dxy1 = " << dca3dxy1 << " dca3dz1 = " << dca3dz1 << std::endl;

      if(!dca3dxy1)
	{
	  std::cout << " get_dca returned NAN " << std::endl;
	  continue;
	}
      if(fabs(dca3dxy1) < _track_dcaxy_cut) continue;
      if(fabs(dca3dz1) < _track_dcaz_cut) continue;
      
      // look for close DCA matches with all other such tracks
      for(auto tr2_it = std::next(tr1_it); tr2_it != _track_map->end(); ++tr2_it)
	{
	  auto id2 = tr2_it->first;
	  auto tr2 = tr2_it->second;

	  // Reverse or remove this to consider TPC tracks too
	  if(_require_mvtx && !hasSiliconSeed(tr2)) continue;

	  int has_silicon_2 = 0;
	  if(hasSiliconSeed(tr2)) has_silicon_2 = 1;

	  if(Verbosity() > 3)
	    {
	      std::cout << "Track2 " << id2 << " details: " << std::endl;
	      outputTrackDetails(tr2);
	    }

	  if(tr2->get_charge() == tr1->get_charge()) continue;

	  if(tr2->get_quality() > _qual_cut) continue;

	  auto tpc2_seed = tr2->get_tpc_seed();
	  int ntpc2 = tpc2_seed->size_cluster_keys();
	  if(ntpc2 < 20) continue;

	  float dca3dxy2, dca3dz2,dca3dxysigma2, dca3dzsigma2;
	  get_dca(tr2, dca3dxy2, dca3dz2, dca3dxysigma2, dca3dzsigma2);
	  //std::cout << " tr2 " << id2 << " dca3dxy2 = " << dca3dxy2 << " dca3dz2 = " << dca3dz2 << std::endl;
	  if(!dca3dxy2)
	    {
	      std::cout << " get_dca returned NAN " << std::endl;
	      continue;
	    }
	  if(fabs(dca3dxy2) < _track_dcaxy_cut) continue;
	  if(fabs(dca3dz2) < _track_dcaz_cut) continue;

	  // find DCA and PCA of these two tracks
	  if(Verbosity() > 3) 
	  { std::cout << "Check pair DCA for tracks " << id1 << " and  " << id2 << std::endl;}

	  double R1;
	  Eigen::Vector2d center1;
	  getCircleXYTrack(tr1, R1, center1);
	  if(Verbosity() > 2)
	    {
	      std::cout << std::endl << "Track 1 circle: " << std::endl;
	      std::cout <<  "  tr1.x " << tr1->get_x() << " tr1.y " << tr1->get_y() << " tr1.z " << tr1->get_z() << " charge " << tr1->get_charge() << std::endl;
	      std::cout << "   tr1.px " << tr1->get_px() << " tr1.py " << tr1->get_py() << " tr1.pz " << tr1->get_pz() << std::endl;	  
	      std::cout << "         track1 " << tr1->get_id() << " circle R " << R1 << " (x, y)  " << center1(0) << "  " << center1(1) << std::endl;
	    }
	  double R2;
	  Eigen::Vector2d center2;
	  getCircleXYTrack(tr2, R2, center2);	  
	  if(Verbosity() > 2)
	    {
	      std::cout << "Track 2 circle: " << std::endl;
	      std::cout << "   tr2.x " << tr2->get_x() << " tr2.y " << tr2->get_y() << " tr2.z " << tr2->get_z() << " charge " << tr2->get_charge() << std::endl;
	      std::cout << "   tr2.px " << tr2->get_px() << " tr2.py " << tr2->get_py() << " tr2.pz " << tr2->get_pz() << std::endl;
	      std::cout << "         track2 " << tr2->get_id() << " circle R " << R2 << " (x, y)  " << center2(0) << "  " << center2(1) << std::endl;
	    }

	  std::vector<double> intersections;
	  if( !circle_circle_intersection(R1, center1(0), center1(1), R2, center2(0), center2(1),  intersections) ) 
	    {
	      if(Verbosity() > 2) std::cout << "    - no intersections, skip this pair" << std::endl;
	      continue;
	    }

	  Eigen::Vector2d intersection[2];
	  if(intersections.size() == 2)
	    {
	      intersection[0](0) = intersections[0];
	      intersection[0](1) = intersections[1];
	    }	      
	  if(intersections.size() == 4)
	    {
	      intersection[1](0) = intersections[2];
	      intersection[1](1) = intersections[3];
	    }
	  
	  // process both intersections
	  for(int i=0; i<2; ++i)
	    {
	      if(intersection[i].norm() == 0) continue;

	      double vradius = sqrt(intersection[i](0)*intersection[i](0) + intersection[i](1) * intersection[i](1));
	      if(vradius > _max_intersection_radius) continue;

	      double z1 = getZFromIntersectionXY(tr1, R1, center1, intersection[i]);
	      double z2 = getZFromIntersectionXY(tr2, R2, center2, intersection[i]);
	      if(Verbosity() > 2)
		std::cout << " track intersection " << i << " at (x,y) " << intersection[i](0) << "  " << intersection[i](1) 
			  << " radius " << vradius << " est. z1 " << z1 << " est. z2 " << z2 << std::endl;

	      if( fabs(z1-z2) > 2.0)   // skip this intersection, it is not good
		{
		  if(Verbosity() > 2) std::cout << "    z-mismatch - wrong intersection, skip it " << std::endl;
		  continue;  
		}

	      Eigen::Vector3d vpos1(0,0,0), vmom1(0,0,0);
	      Eigen::Vector3d vpos2(0,0,0), vmom2(0,0,0);

	      // Project the tracks to the intersection
	      Eigen::Vector3d intersect1(intersection[i](0), intersection[i](1), z1);
	      if(!projectTrackToPoint(tr1, intersect1, vpos1, vmom1)) continue;
	      if(Verbosity() > 2) 
		{
		  std::cout << "  Projected track 1 to point " << intersect1(0) << "  " << intersect1(1) << "  " << intersect1(2) << std::endl; 
		  std::cout <<  "                                 has vpos " << vpos1(0) << "  " << vpos1(1) << "  " << vpos1(2) << std::endl;	      
		}
	      Eigen::Vector3d intersect2(intersection[i](0), intersection[i](1), z2);
	      if(!projectTrackToPoint(tr2, intersect2, vpos2, vmom2)) continue;
	      if(Verbosity() > 2) 
		{
		 std::cout << "  Projected track 2 to point " << intersect2(0) << "  " << intersect2(1) << "  " << intersect2(2) << std::endl; 
		 std::cout << "                                 has vpos " << vpos2(0) << "  " << vpos2(1) << "  " << vpos2(2) << std::endl;	      
		 }
	      
	      // check that the z positions are close
	      if(fabs(vpos1(2) - vpos2(2)) > _projected_track_z_cut) 
		{
		  if(Verbosity() > 0) 
		    { std::cout << "    Warning: projected z positions are screwed up, should not be" << std::endl; }
		  continue;
		}

	      if(Verbosity() > 2)
		{
		  std::cout << "Summary for projected pair:" << std::endl;
		  std::cout << " Fitted tracks: " << std::endl;
		  std::cout << "   tr1.x " << tr1->get_x() << " tr1.y " << tr1->get_y() << " tr1.z " << tr1->get_z() << std::endl;
		  std::cout << "   tr2.x " << tr2->get_x() << " tr2.y " << tr2->get_y() << " tr2.z " << tr2->get_z() << std::endl;
		  std::cout << "   tr1.px " << tr1->get_px() << " tr1.py " << tr1->get_py() << " tr1.pz " << tr1->get_pz() << std::endl;
		  std::cout << "   tr2.px " << tr2->get_px() << " tr2.py " << tr2->get_py() << " tr2.pz " << tr2->get_pz() << std::endl;
		  std::cout << " Projected tracks: " << std::endl;
		  std::cout << "   pos1.x " << vpos1(0) << " pos1.y " << vpos1(1) << " pos1.z " << vpos1(2) << std::endl; 
		  std::cout << "   pos2.x " << vpos2(0) << " pos2.y " << vpos2(1) << " pos2.z " << vpos2(2) << std::endl; 
		  std::cout << "   mom1.x " << vmom1(0) << " mom1.y " << vmom1(1) << " mom1.z " << vmom1(2) << std::endl; 
		  std::cout << "   mom2.x " << vmom2(0) << " mom2.y " << vmom2(1) << " mom2.z " << vmom2(2) << std::endl; 
		}

	      // Improve the pair dca using a local straight line approximation	      	      
	      double pair_dca;
	      Eigen::Vector3d PCA1(0,0,0), PCA2(0,0,0);
	      findPcaTwoLines(vpos1, vmom1, vpos2, vmom2, pair_dca, PCA1, PCA2);  
	      if(Verbosity() > 2 ) 
		{
		  std::cout << "       pair_dca " << pair_dca << " two_track_dcacut " << _two_track_dcacut << std::endl;
		  std::cout << "       PCA1 " << PCA1(0) << "  " << PCA1(1) << "  " << PCA1(2) << std::endl;
		  std::cout << "       PCA2 " << PCA2(0) << "  " << PCA2(1) << "  " << PCA2(2) << std::endl;
		}

	      if(fabs(pair_dca) > _two_track_dcacut) { continue; }
	      
	      // calculate the invariant mass using the track states at the decay vertex

	      TLorentzVector t1;
	      Float_t E1 = sqrt(pow(vmom1(0),2) + pow(vmom1(1),2) + pow(vmom1(2),2) 
				+ pow(_decaymass,2));
	      t1.SetPxPyPzE(vmom1(0),vmom1(1),vmom1(2),E1);	
	      
	      TLorentzVector t2;
	      Float_t E2 = sqrt(pow(vmom2(0),2) + pow(vmom2(1),2) + pow(vmom2(2),2) 
				+ pow(_decaymass,2));
	      t2.SetPxPyPzE(vmom2(0),vmom2(1),vmom2(2),E2);	
	      
	      TLorentzVector tsum = t1+t2;
	      
	      // calculate the decay length
	      Eigen::Vector3d PCA = (vpos1+vpos2) / 2.0;  // average the PCA of the track pair
	      auto vtxid = tr1->get_vertex_id();
	      auto vertex1 = _svtx_vertex_map->get(vtxid);
	      Eigen::Vector3d VTX(vertex1->get_x(), vertex1->get_y(), vertex1->get_z());
	      Eigen::Vector3d path = PCA - VTX;
	      double decay_radius = sqrt( pow(PCA(0),2) + pow(PCA(1),2) );
	      
	      if(path.norm() > _min_path_cut)
		{
		  if(Verbosity() > 0)
		    {
		      std::cout << "    Pair mass " << tsum.M() << " pair pT " << tsum.Pt() 
				<< " decay length " << path.norm() << std::endl;
		    }

		  if(_write_ntuple)
		    {
		      fillNtp(tr1, tr2,  dca3dxy1, dca3dz1, dca3dxy2, dca3dz2, vpos1, vmom1, vpos2, vmom2,
			      PCA1, PCA2, pair_dca, tsum.M(), tsum.Pt(), path.norm(), has_silicon_1, has_silicon_2);
		    }

		  if(_write_electrons_node)
		    {
		      if(passConversionElectronCuts(tsum, tr1, tr2, pair_dca, PCA, VTX))
			{
			  if(Verbosity() > 0)
			    std::cout << "     **** inserting tracks " << tr1->get_id() << "  and " << tr2->get_id() << std::endl; 

			  // Add decay particles to output node
			  _track_map_electrons->insertWithKey(tr1, tr1->get_id());
			  _track_map_electrons->insertWithKey(tr2, tr2->get_id());

			  if(_write_ntuple)
			    { 
			      // these are just to check on the effect of the cuts
			      recomass->Fill(tsum.Pt(), tsum.M());  
			      hdecay_radius->Fill(decay_radius);
			      hdecaypos->Fill( PCA(0), PCA(1));
			    }
			}
		    }
		}
	    }
	}

    }

  if(Verbosity() > 0)
    {
      std::cout << PHWHERE << " electron track map size " << _track_map_electrons->size()  << std::endl;  
      for(auto tr_it = _track_map_electrons->begin(); tr_it != _track_map_electrons->end(); ++tr_it)
	{
	  auto id = tr_it->first;
	  auto tr = tr_it->second;

	  std::cout << " Electron track " << id 
		    << " x " << tr->get_x() << " y " << tr->get_y() << " z " << tr->get_z() 
		    << " px " << tr->get_px() << " py " << tr->get_py() << " pz " << tr->get_pz() 
		    << std::endl;	  
	}
    }

  return Fun4AllReturnCodes::EVENT_OK;
}

bool SecondaryVertexFinder::passConversionElectronCuts(TLorentzVector tsum, SvtxTrack* tr1, SvtxTrack* tr2, float pair_dca, Eigen::Vector3d PCA, Eigen::Vector3d VTX)
{
  bool pass = false;

  // additional single track cuts
  auto tpc_seed1 = tr1->get_tpc_seed();
  unsigned int ntpc1 = tpc_seed1->size_cluster_keys();
  auto tpc_seed2 = tr2->get_tpc_seed();
  unsigned int ntpc2 = tpc_seed2->size_cluster_keys();
  if(ntpc1 < _min_tpc_clusters || ntpc2 < _min_tpc_clusters) {return pass;}

  // pair cuts
  if(fabs(pair_dca) > _conversion_pair_dcacut) {return pass;};

  if(tsum.M() > _max_mass_cut) {return pass;}

  if(tsum.Pt() < _invariant_pt_cut) {return pass;}

  Eigen::Vector3d path = PCA - VTX;
  if(path.norm() < 0.2) {return pass;}

  // Angle between path vector and reco pair momentum vector
  Eigen::Vector3d invariant_mom(tsum.Px(), tsum.Py(), tsum.Pz());
  double costheta = path.dot(invariant_mom) / (path.norm() * invariant_mom.norm());
  //      std::cout << " costheta " << costheta << std::endl;
  if(costheta < _costheta_cut) {return pass;};
  
  // selects "decays" with zero mass
  if(fabs(tr1->get_eta() - tr2->get_eta()) > _deta_cut) {return pass;}; 


  pass = true;
  return pass;
}

double SecondaryVertexFinder::getZFromIntersectionXY(SvtxTrack *track, double& R, Eigen::Vector2d& center, Eigen::Vector2d intersection)
{
  // Starting at the track base position, the path in the XY plane and the path in the Z direction are proportional
  // They are related by:   dpath/dz = pT/pz

  // The xy path to the intersection is an arc due to a rotation of pT:
  double phi0 = atan2(track->get_y() - center(1), track->get_x() - center(0));
  double phi1 = atan2(intersection(1) - center(1), intersection(0) - center(0));
  double xypath = R * fabs(phi1-phi0);

  double zpath = xypath * track->get_pz() / track->get_pt();
  double z = track->get_z() + zpath;
  if(Verbosity() > 3)
    {
      std::cout << " Radius " << R << " center xy " << center(0) << "  " << center(1) << "  pT " << track->get_pt() << " pz " << track->get_pz() << std::endl;
      std::cout << " phi0 " << phi0 << " y0 " << track->get_y() << " x0 " << track->get_x() << std::endl; 
      std::cout << " phi1 " << phi1 << " intersection y " << intersection(1) << " intersection x " << intersection(0) << std::endl;
      std::cout << " xypath " << xypath << " zpath " << zpath << " z0 " << track->get_z() << " z " << z << std::endl; 
    }

  return z;
}

void SecondaryVertexFinder::getCircleXYTrack(SvtxTrack *track, double& R, Eigen::Vector2d& center)
{
  // Get the circle equation for the track
  double Bz = 1.4; // T
  double x = track->get_x();
  double y = track->get_y();
  double px = track->get_px();
  double py = track->get_py();
  double pt = track->get_pt();
  int charge = track->get_charge();

  // radius of curvature is from pT and B field
  R = 3.3 * pt / (Bz * (double) charge) * 100;  // convert from m to cm, sign changes for negative charge

  // the normal unit vector in (x,y) space at (x,y) is:
  Eigen::Vector2d normal(py, -px);
  normal /= normal.norm();

  // The circle center is at
  Eigen::Vector2d base(x,y);
  center = base + R*normal;

  // needed the sign to get the direction of the center, now we drop it
  R = fabs(R);

  return;

}

bool SecondaryVertexFinder::projectTrackToPoint(SvtxTrack* track, Eigen::Vector3d& PCA, Eigen::Vector3d& pos, Eigen::Vector3d& mom)
{

  bool ret = true;
  PCA *= Acts::UnitConstants::cm;  
 
  /// create perigee surface
  auto perigee = Acts::Surface::makeShared<Acts::PerigeeSurface>(PCA);

  ActsPropagator propagator(_tGeometry);
  auto params = propagator.makeTrackParams(track, _svtx_vertex_map);
  if(!params.ok())
    {
      /// Couldn't construct strict PCA track parameter requirement
      return false;
    }

  propagator.verbosity(Verbosity());
  auto result = propagator.propagateTrack(params.value(), perigee);	\
  if(result.ok())
    {
      auto resultparams = result.value().second;
      auto projectionPos = resultparams.position(_tGeometry->geometry().getGeoContext());
      const auto momentum = resultparams.momentum();
      pos(0) = projectionPos.x() / Acts::UnitConstants::cm;
      pos(1) = projectionPos.y() / Acts::UnitConstants::cm;
      pos(2) = projectionPos.z() / Acts::UnitConstants::cm;
      
      mom(0) = momentum.x();
      mom(1) = momentum.y();
      mom(2) = momentum.z();	      
    }
  else
    {
      if(Verbosity() > 0)
	{
	  std::cout << "   Failed projection of track with: " << std::endl;
	  std::cout << " x,y,z = " << track->get_x() << "  " << track->get_y() << "  " << track->get_z() << std::endl;
	  std::cout << " px,py,pz = " << track->get_px() << "  " << track->get_py() << "  " << track->get_pz() << std::endl;
	  std::cout << " to point (x,y,z) = " << PCA(0) / Acts::UnitConstants::cm << "  " 
		    << PCA(1) / Acts::UnitConstants::cm << "  " << PCA(2) / Acts::UnitConstants::cm << std::endl; 
	}

      ret = false;
    }

  return ret;

}

void SecondaryVertexFinder::outputTrackDetails(SvtxTrack *tr)
{
  auto tpc_seed = tr->get_tpc_seed();
  int ntpc = tpc_seed->size_cluster_keys();

  auto silicon_seed = tr->get_silicon_seed();

  int nsilicon = 0;
  if(silicon_seed) 
    { nsilicon = silicon_seed->size_cluster_keys(); }

  auto pt = tr->get_pt();
  auto eta = tr->get_eta();
  auto x = tr->get_x();
  auto y = tr->get_y();
  auto z = tr->get_z();
  auto qual = tr->get_quality();

  std::cout << "   ntpc " << ntpc << " nsilicon " << nsilicon << " quality " << qual 
	    << " eta " << eta << std::endl; 
  std::cout << "   pt " << pt << " x " << x << " y " << y << " z " << z << std::endl;

  auto vtxid = tr->get_vertex_id();
  auto vertex = _svtx_vertex_map->get(vtxid);
  std::cout << "    vtxid " << vtxid 
	    << " vertex x " << vertex->get_x()
	    << " vertex y " << vertex->get_y()
	    << " vertex z " << vertex->get_z()
	    << std::endl;

}

bool  SecondaryVertexFinder::hasSiliconSeed(SvtxTrack* tr) 
{
  bool ret = false;
  auto silicon_seed = tr->get_silicon_seed();
  if(silicon_seed) ret = true;

  return ret;
  }

void SecondaryVertexFinder::get_dca(SvtxTrack* track, 
			    float& dca3dxy, float& dca3dz,
			    float& dca3dxysigma, float& dca3dzsigma)
{
  dca3dxy = NAN;
  Acts::Vector3 pos(track->get_x(),
		    track->get_y(),
		    track->get_z());
  Acts::Vector3 mom(track->get_px(),
		    track->get_py(),
		    track->get_pz());

  auto vtxid = track->get_vertex_id();
  auto svtxVertex = _svtx_vertex_map->get(vtxid);
  if(!svtxVertex)
    { 
      std::cout << "   Failed to find vertex for track " << std::endl;
      return; 
    }

  Acts::Vector3 vertex(svtxVertex->get_x(),
		       svtxVertex->get_y(),
		       svtxVertex->get_z());

  if(Verbosity() > 3)
    {  
      std::cout << "   track " << track->get_id() << " vertex id is " << vtxid 
		<< " vertex is " << svtxVertex->get_x() << "  " 
		<< svtxVertex->get_y() << "  "
		<< svtxVertex->get_z() << std::endl;
    }

  pos -= vertex;

  Acts::ActsSquareMatrix<3> posCov;
  for(int i = 0; i < 3; ++i)
    {
      for(int j = 0; j < 3; ++j)
	{
	  posCov(i, j) = track->get_error(i, j);
	} 
    }
  
  Acts::Vector3 r = mom.cross(Acts::Vector3(0.,0.,1.));
  float phi = atan2(r(1), r(0));
  phi *= -1.0;  // rotates vector clockwise to x axis
  
  Acts::RotationMatrix3 rot;
  Acts::RotationMatrix3 rot_T;
  rot(0,0) = cos(phi);
  rot(0,1) = -sin(phi);
  rot(0,2) = 0;
  rot(1,0) = sin(phi);
  rot(1,1) = cos(phi);
  rot(1,2) = 0;
  rot(2,0) = 0;
  rot(2,1) = 0;
  rot(2,2) = 1;
  
  rot_T = rot.transpose();

  Acts::Vector3 pos_R = rot * pos;
  Acts::ActsSquareMatrix<3> rotCov = rot * posCov * rot_T;

  dca3dxy = pos_R(0);
  dca3dz = pos_R(2);
  dca3dxysigma = sqrt(rotCov(0,0));
  dca3dzsigma = sqrt(rotCov(2,2));
}

void SecondaryVertexFinder::findPcaTwoLines(Eigen::Vector3d pos1, Eigen::Vector3d mom1, Eigen::Vector3d pos2, Eigen::Vector3d mom2,
			double &dca, Eigen::Vector3d &PCA1, Eigen::Vector3d &PCA2)
{
  Eigen::Vector3d a1(pos1(0), pos1(1), pos1(2));
  Eigen::Vector3d b1(mom1(0) / mom1.norm(), mom1(1) / mom1.norm(), mom1(2) / mom1.norm());
  Eigen::Vector3d a2(pos2(0), pos2(1), pos2(2));
  Eigen::Vector3d b2(mom2(0) / mom2.norm(), mom2(1) / mom2.norm(), mom2(2) / mom2.norm());

  // The shortest distance between two skew lines described by
  //  a1 + c * b1,   a2 + d * b2
  // where a1, a2, are vectors representing points on the lines, b1, b2 are direction vectors, 
  //  and c and d are scalars, is:
  // dca = (b1 x b2) .(a2-a1) / |b1 x b2|

  // bcrossb/mag_bcrossb is a unit vector perpendicular to both direction vectors b1 and b2
  auto bcrossb = b1.cross(b2);
  auto mag_bcrossb = bcrossb.norm();
  // a2-a1 is the vector joining any arbitrary points on the two lines
  auto aminusa = a2-a1;

  // The DCA of these two lines is the projection of a2-a1 along the direction of the perpendicular to both 
  // remember that a2-a1 is longer than (or equal to) the dca by definition
  dca = 999;
  if( mag_bcrossb != 0)
    dca = bcrossb.dot(aminusa) / mag_bcrossb;
  else
    return;   // same track, skip combination
  
  // get the points at which the normal to the lines intersect the lines, where the lines are perpendicular

  double X =  b1.dot(b2) - b1.dot(b1) * b2.dot(b2) / b2.dot(b1);
  double Y =  (a2.dot(b2) - a1.dot(b2)) - (a2.dot(b1) - a1.dot(b1)) * b2.dot(b2) / b2.dot(b1) ;
  double c = Y/X;

  double F = b1.dot(b1) / b2.dot(b1); 
  double G = - (a2.dot(b1) - a1.dot(b1)) / b2.dot(b1);
  double d = c * F + G;

  // then the points of closest approach are:
  PCA1 = a1+c*b1;
  PCA2 = a2+d*b2;

  return;
}



//===========================
// replace with a call to TrackFitUtils
//===========================
bool SecondaryVertexFinder::circle_circle_intersection(double r0, double x0, double y0, double r1, double x1, double y1, std::vector<double>& intersectionXY)
{
  bool ret = true;
 
  Eigen::Vector2d p0(x0,y0);
  Eigen::Vector2d p1(x1,y1);

  double d = (p0 - p1).norm();

  //  std::cout << "   d " << d << " r1-r0 " << fabs(r1-r0) << std::endl;

  // no intersection possible
  if(d < fabs(r1-r0)) return false;  // one circle inside the other
  if(d > r0+r1)
    { 
      // careful: conversion electrons will look like zero mass decays
      // fluctuations may cause the circles to (just) not touch - what to do about that?
      //   if d - (r0+r1) < dr,   then there is only one PCA, and it is on the line between the two centers
      double dr = 0.2;  // 2 mm
      if( fabs(d - (r0+r1)) < dr)
	{
	  // find the closest point on circle 0 to the center of circle 1
	  Eigen::Vector2d u0 = (p1 - p0);
	  u0 /= u0.norm();
	  Eigen::Vector2d PCA0 = p0 + u0 * r0;

 	  Eigen::Vector2d u1 = (p0 - p1);
	  u1 /= u1.norm();
	  Eigen::Vector2d PCA1 = p1 + u1 * r1;

	  auto PCA = (PCA0+PCA1) / 2.0;
	  intersectionXY.push_back(PCA(0));
	  intersectionXY.push_back(PCA(1));

	  if(Verbosity() > 2) std::cout << "      *** Special case: Barely touching circles: " << " PCA.x, PCA.y " << PCA(0) << "   " << PCA(1) << std::endl; 
	  return ret;
	}
      else
	return false;  
    }

  double a=(r0*r0-r1*r1+d*d)/(2*d);
  double h = sqrt(r0*r0 - a*a);

  double x2= x0 + a*(x1-x0)/d;   
  double y2=y0+a*(y1-y0)/d;

  double x3a=x2+h*(y1-y0)/d;       // also x3=x2-h*(y1-y0)/d
  double y3a=y2-h*(x1-x0)/d;       // also y3=y2+h*(x1-x0)/d

  double x3b=x2-h*(y1-y0)/d;       // also x3=x2-h*(y1-y0)/d
  double y3b=y2+h*(x1-x0)/d;       // also y3=y2+h*(x1-x0)/d

    intersectionXY.push_back(x3a);
    intersectionXY.push_back(y3a);
    intersectionXY.push_back(x3b);
    intersectionXY.push_back(y3b);

    return ret;

}

int SecondaryVertexFinder::End(PHCompositeNode */*topNode*/)
{
  if(_write_ntuple)
    {
      TFile *fout = new TFile(outfile.c_str(),"recreate");
      recomass->Write();
      hdecaypos->Write();
      hdecay_radius->Write();
      ntp->Write();
      fout->Close();
    }

  return Fun4AllReturnCodes::EVENT_OK;
}

int SecondaryVertexFinder::CreateOutputNode(PHCompositeNode* topNode)
{

  PHNodeIterator iter(topNode);
  
  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST"));

  if (!dstNode)
  {
    std::cerr << "DST node is missing, quitting" << std::endl;
    throw std::runtime_error("Failed to find DST node in SecondaryVertexFinder");
  }
  
  PHNodeIterator dstIter(topNode);
  PHCompositeNode *svtxNode = dynamic_cast<PHCompositeNode *>(dstIter.findFirst("PHCompositeNode", "SVTX"));

  if (!svtxNode)
  {
    svtxNode = new PHCompositeNode("SVTX");
    dstNode->addNode(svtxNode);
  }

  if(_write_electrons_node)
  {
    _track_map_electrons = new SvtxTrackMap_v2;
    auto tracks_node = new PHIODataNode<PHObject>( _track_map_electrons, "SvtxTrackMapElectrons", "PHObject");
    svtxNode->addNode(tracks_node);
    if (Verbosity())
    { std::cout << "Svtx/SvtxTrackMapElectrons node added" << std::endl; }
  }

  return true;
}

int SecondaryVertexFinder::GetNodes(PHCompositeNode* topNode)
{

  _track_map = findNode::getClass<SvtxTrackMap>(topNode, "SvtxTrackMap");
  if (!_track_map)
  {
    std::cout << PHWHERE << " ERROR: Can't find SvtxTrackMap: " << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  /*
  _cluster_map = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");
  if (!_cluster_map)
    {
      std::cout << PHWHERE << " ERROR: Can't find node TRKR_CLUSTER" << std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }
  */

  _svtx_vertex_map = findNode::getClass<SvtxVertexMap>(topNode,"SvtxVertexMap");  
  if(!_svtx_vertex_map)
    {
      std::cout << PHWHERE << "No vertex node, quit!" << std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    } 

  _tGeometry = findNode::getClass<ActsGeometry>(topNode,"ActsGeometry");
  if(!_tGeometry)
    {
      std::cout << PHWHERE << "Error, can't find acts tracking geometry" << std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }

  return Fun4AllReturnCodes::EVENT_OK;
}


