#include <trackbase/TrkrDefs.h>
#include "TruthRecoTrackMatching.h"

#include <trackbase/TrkrTruthTrackContainer.h>
#include <trackbase/TrkrTruthTrack.h>

#include <trackbase_historic/SvtxTrack.h>
#include <trackbase_historic/SvtxTrackMap.h>

#include <g4detectors/PHG4TpcCylinderGeom.h>
#include <g4detectors/PHG4TpcCylinderGeomContainer.h>

#include <trackbase/TrkrCluster.h>
#include <trackbase/TrkrClusterContainer.h>

#include <trackbase/EmbRecoMatch.h>
#include <trackbase/EmbRecoMatchv1.h>

#include <trackbase/EmbRecoMatchContainer.h>
#include <trackbase/EmbRecoMatchContainerv1.h>

// not actually sure if I need all of these
#include <phool/PHCompositeNode.h>
#include <phool/PHDataNode.h>  // for PHDataNode
#include <phool/PHIODataNode.h>
#include <phool/PHNode.h>  // for PHNode
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>  // for PHObject
#include <phool/PHRandomSeed.h>
#include <phool/getClass.h>
#include <phool/phool.h>  // for PHWHERE

#include <g4main/PHG4TruthInfoContainer.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <TSystem.h>
#include <algorithm>

typedef TruthRecoTrackMatching::RECO_pair_iter RECO_pair_iter;
typedef TruthRecoTrackMatching::PossibleMatch  PossibleMatch;
typedef TruthRecoTrackMatching::CompCluster    CompCluster;

/* using vector<pair<short,short>> = VecShortPairs; */
using std::make_tuple;
using std::make_pair;
using std::vector;
using std::get;
using std::array;
using std::endl;
using std::cout;

TruthRecoTrackMatching::TruthRecoTrackMatching(
      const unsigned short _nmincluster_match,
      const float _nmincluster_ratio,
      const double _cutoff_dphi,
      const double _same_dphi,
      const double _cutoff_deta,
      const double _same_deta,
      const double _cluster_nzwidths,
      const double _cluster_nphiwidths,
      const std::string _filename )
        : m_nmincluster_match    { _nmincluster_match  } // minimum number of clusters to match, default=4
        , m_nmincluster_ratio    { _nmincluster_ratio  } // minimum ratio to match a track, default=0.
                                                         // -- Track Kinematic Cuts to match --
        , m_cutoff_dphi          { _cutoff_dphi        } // maximum |dphi|=|phi_reco-phi_truth| to search for match
        , m_same_dphi            { _same_dphi          } // all tracks in this |dphi| must be tested for matches
        , m_cutoff_deta          { _cutoff_deta        } // same as m_cutoff_dphi for deta
        , m_same_deta            { _same_deta          } // same as m_same_dphi for deta
        , m_cluster_nzwidths_2   { _cluster_nzwidths * _cluster_nzwidths } 
        // match TrkrClusters_reco within this ratio of TrkrCluster_true z widths
        , m_cluster_nphiwidths_2 { _cluster_nphiwidths * _cluster_nphiwidths } 
        // same as m_cluster_nzwidths for phi
{ 
  m_maketree = (_filename != "");
  if (m_maketree) {
    m_fileout = new TFile(_filename.c_str(),"recreate");
    m_fileout->cd();
    m_treeout = new TTree("TrackClusters", "track_clusters");

    //track kinematics
    m_treeout->Branch("reco_eta", &b_reco_eta);
    m_treeout->Branch("reco_phi", &b_reco_phi);
    m_treeout->Branch("reco_pt",  &b_reco_pt);

    m_treeout->Branch("truth_eta", &b_truth_eta);
    m_treeout->Branch("truth_phi", &b_truth_phi);
    m_treeout->Branch("truth_pt",  &b_truth_pt);

    //truth kinematics
    m_treeout->Branch("reco_cl_layer",    &b_reco_cl_layer);
    m_treeout->Branch("reco_cl_phi",      &b_reco_cl_phi);
    m_treeout->Branch("reco_cl_phiwidth", &b_reco_cl_phiwidth);
    m_treeout->Branch("reco_cl_z",        &b_reco_cl_z);
    m_treeout->Branch("reco_cl_zwidth",   &b_reco_cl_zwidth);

    m_treeout->Branch("truth_cl_layer",    &b_truth_cl_layer);
    m_treeout->Branch("truth_cl_phi",      &b_truth_cl_phi);
    m_treeout->Branch("truth_cl_phiwidth", &b_truth_cl_phiwidth);
    m_treeout->Branch("truth_cl_z",        &b_truth_cl_z);
    m_treeout->Branch("truth_cl_zwidth",   &b_truth_cl_zwidth);
  }
}

int TruthRecoTrackMatching::InitRun(PHCompositeNode *topNode) //`
{ 
  if (createNodes(topNode) != Fun4AllReturnCodes::EVENT_OK)
  {
    return Fun4AllReturnCodes::ABORTRUN;
  }
  return Fun4AllReturnCodes::EVENT_OK;
  topNode->print();
}


int TruthRecoTrackMatching::process_event(PHCompositeNode* topnode)  //`
{
  if (topnode == nullptr) {
    std::cout << " topnode is null " << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  /* topnode->print(); */

  // -------------------------------------------------------------------------------
  // Build recoData
  // ------------------------------------------------------------------------------
  // Generate a vector use for matching criteria
  if (false) std::cout << "reco tracks size: " << m_SvtxTrackMap->size() << std::endl; // FIXME

  recoData.clear();
  for (auto reco = m_SvtxTrackMap->begin(); reco != m_SvtxTrackMap->end(); ++reco) {
    auto index_reco = reco->first;
    auto track = reco->second;
    /* if (true) std::cout << Form(" Reconstructed Track in SvtxTrackMap: id(%2i-%2i)   phi(%5.2f) eta(%5.2f) pt(%5.2f)", // FIXME */
        /* reco->first, track->get_id(), track->get_phi(), track->get_eta(), track->get_pt()) << std::endl; */
    recoData.push_back( {track->get_phi(), track->get_eta(), track->get_pt(), index_reco } );
  }
  std::sort(recoData.begin(), recoData.end(), CompRECOtoPhi());

  // wrap low-to-high phi entries
  RECOvec wraps{};
  auto wrap_from_start = std::upper_bound(recoData.begin(), recoData.end(), (-M_PI+m_cutoff_dphi),  CompRECOtoPhi());
  for (auto iter = recoData.begin(); iter != wrap_from_start; ++iter) {
    auto entry = *iter; // make a new copy to wrap to the other side of recoData
    get<RECOphi>(entry) = get<RECOphi>(entry) + 2*M_PI; // put the new copy on the other end
    wraps.push_back(entry);
  }
  auto wrap_to_end = std::lower_bound(recoData.begin(), recoData.end(),  (M_PI-m_cutoff_dphi), CompRECOtoPhi());
  for (auto iter = wrap_to_end; iter != recoData.end(); ++iter) {
    auto entry = *iter;
    get<RECOphi>(entry) = get<RECOphi>(entry) - 2*M_PI;
    wraps.push_back(entry);
  }
  for (auto E : wraps) recoData.push_back(E);

  // sort again to put the sandwiching phi-wrapped entries at the beginning and end in Phi
  std::sort(recoData.begin(), recoData.end(), CompRECOtoPhi());

  if (false) {//FIXME
    cout << " ************************************* " << endl;
    cout << " This is the map of tracks to match to " << endl;
    cout << " ************************************* " << endl;
    for (auto& E : recoData) { cout << Form(" id:%2i  (phi:eta:pt)  (%5.2f:%5.2f:%5.2f)", get<RECOid>(E), get<RECOphi>(E), get<RECOeta>(E), get<RECOpt>(E)) << endl; }
    cout << " ****end*listing*********************** " << endl;
  }

  // -------------------------------------------------------------------------------
  // Loop through the truth tracks
  //   - find which ones may have matches (truthid_to_match)
  //   - get PossibleMatch values for all innterbox matches
  //   - collect outerbox matches
  // ------------------------------------------------------------------------------

  std::vector<std::pair<unsigned short, unsigned short>> outerbox_indices {};
  std::vector<std::pair<unsigned short, unsigned short>> innerbox_indices {};

  if (false) std::cout << "Number of truth tracks: " << m_TrkrTruthTrackContainer->getTruthTracks().size() << std::endl;
  unsigned short index_true {0};
  std::map<unsigned int,bool> map_matched;
  for (auto track : m_TrkrTruthTrackContainer->getTruthTracks()) {
    auto match_indices = find_box_matches(track->getPhi(), track->getPseudoRapidity(), track->getPt()); 
    if (false) { //FIXME
      std::cout << Form(" Embedded Track:  (%2i)   phi(%5.2f) eta(%5.2f) pt(%5.2f)", track->getTrackid(), track->getPhi(),
          track->getPseudoRapidity(), track->getPt()) << std::endl;
      cout << " matched: (inner box) ";
      for (auto match : match_indices.first) cout << " IM(" << match <<")";
      cout << " (outer box) ";
      for (auto match : match_indices.second) cout << " OM(" << match<<")";
      cout << endl;
    }
    map_matched[track->getTrackid()]=false;
    for (auto& id_reco : match_indices.first)  innerbox_indices.push_back( {index_true, id_reco} );
    for (auto& id_reco : match_indices.second) outerbox_indices.push_back( {index_true, id_reco} );
    ++index_true;
  }

  std::set<int> idreco_matched {};
  std::set<int> idtrue_matched {};
  
  if (innerbox_indices.size() > 0)
  match_tracks_in_box (innerbox_indices, idtrue_matched, idreco_matched, true);
  if (outerbox_indices.size() > 0)
  match_tracks_in_box (outerbox_indices, idtrue_matched, idreco_matched, false);

  m_EmbRecoMatchContainer->sort();

  // make a list of embedded tracks which are unmatched
  auto matched   = m_EmbRecoMatchContainer->ids_TruthMatched();
  for (auto &id : matched) {
    map_matched[id] = true;
  }
  auto unmatched = m_EmbRecoMatchContainer->ids_TruthUnmatched();
  for (auto entry : map_matched) {
    if (!entry.second) unmatched.push_back(entry.first);
  }

  if (true) { // FIXME
  // re-print all the tracks with the matches with the fit values
    for (auto& track : m_TrkrTruthTrackContainer->getTruthTracks()) {
      auto id = track->getTrackid();
      if (m_EmbRecoMatchContainer->hasTruthMatch(id)) {
        cout << " Matched " << id;
      }
    }
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

int TruthRecoTrackMatching::End(PHCompositeNode * /*topNode*/)
{
  if (m_maketree) {
    std::cout << "FIXME location 0 " << std::endl;
    m_fileout->cd();
    m_treeout->Write();
    m_fileout->Save();
    m_fileout->Close();
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

  //--------------------------------------------------
  // Internal functions
  //--------------------------------------------------

int TruthRecoTrackMatching::createNodes(PHCompositeNode* topNode)
{
  // Initiailize the modules
  m_PHG4TruthInfoContainer = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");
  if (!m_PHG4TruthInfoContainer)
  {
    std::cout << "Could not locate G4TruthInfo node when running \"TruthRecoTrackMatching\" module." << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  m_SvtxTrackMap = findNode::getClass<SvtxTrackMap>(topNode, "SvtxTrackMap");
  if (!m_SvtxTrackMap)
  {
    std::cout << "Could not locate SvtxTrackMap node when running \"TruthRecoTrackMatching\" module." << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }


  m_TruthClusterContainer = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_TRUTHCLUSTERCONTAINER");
  if (!m_TruthClusterContainer)
  {
    std::cout << "Could not locate TRKR_TRUTHCLUSTERCONTAINER node when running \"TruthRecoTrackMatching\" module." << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  m_RecoClusterContainer = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");
  if (!m_RecoClusterContainer)
  {
    std::cout << "Could not locate TRKR_CLUSTER node when running \"TruthRecoTrackMatching\" module." << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  m_TrkrTruthTrackContainer = findNode::getClass<TrkrTruthTrackContainer>(topNode, "TRKR_TRUTHTRACKCONTAINER");
  if (!m_TrkrTruthTrackContainer)
  {
    std::cout << "Could not locate TRKR_TRUTHTRACKCONTAINER node when running \"TruthRecoTrackMatching\" module." << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  m_PHG4TpcCylinderGeomContainer = 
      findNode::getClass<PHG4TpcCylinderGeomContainer>(topNode, "CYLINDERCELLGEOM_SVTX");
  if (!m_PHG4TpcCylinderGeomContainer)
  {
    std::cout << "Could not locate CYLINDERCELLGEOM_SVTX node when running \"TruthRecoTrackMatching\" module." << std::endl;
    /* std::cout << PHWHERE << "ERROR: Can't find node CYLINDERCELLGEOM_SVTX" << std::endl; */
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  // note that layers 0-6, and > 55, don't work
  for (int layer=7; layer<55; ++layer) {
    PHG4TpcCylinderGeom *layergeom = m_PHG4TpcCylinderGeomContainer->GetLayerCellGeom(layer);
    if (layer==7) m_zstep_2 = pow(layergeom->get_zstep(),2.);
    m_phistep_2[layer] = pow(layergeom->get_phistep(),2.);
  }

  PHNodeIterator iter(topNode);
  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    std::cout << PHWHERE << "DST Node missing, doing nothing." << std::endl;
    exit(1);
  }
  m_EmbRecoMatchContainer = findNode::getClass<EmbRecoMatchContainer>(topNode,"TRKR_EMBRECOMATCHCONTAINER");
  if (!m_EmbRecoMatchContainer) {
    PHNodeIterator dstiter(dstNode);

    auto DetNode = dynamic_cast<PHCompositeNode *>(dstiter.findFirst("PHCompositeNode", "TRKR"));
    if (!DetNode)
    {
      DetNode = new PHCompositeNode("TRKR");
      dstNode->addNode(DetNode);
    }

    m_EmbRecoMatchContainer = new EmbRecoMatchContainerv1;
    auto newNode = new PHIODataNode<PHObject>(m_EmbRecoMatchContainer, "TRKR_EMBRECOMATCHCONTAINER", "PHObject");
    DetNode->addNode(newNode);
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

std::pair<std::vector<unsigned short>, std::vector<unsigned short>> 
TruthRecoTrackMatching::find_box_matches(float truth_phi, float truth_eta, float truth_pt) {
  // sort through the recoData to find:
  //     inner_box : possible matches withing m_same_d{phi,eta,pt}
  //     outer_box : (assumed to be strictly larger than, and contain, inner_box, in m_cutoff_d{phi, eta,pt}
  // recoData is already sorted by eta, but as the boxes "zoom in" on acceptable matches, it gets
  // sorted by eta and (if applicable) pT as well. `mix_pair` keeps track of the range of recoData that has
  // been sorted so that it can be resorted to phi order.
  RECO_pair_iter mix_pair { recoData.begin(), recoData.end() };

  mix_pair.first = std::lower_bound(mix_pair.first, mix_pair.second, 
      truth_phi-m_cutoff_dphi, CompRECOtoPhi());

  mix_pair.second = std::upper_bound(mix_pair.first, mix_pair.second,
      truth_phi+m_cutoff_dphi, CompRECOtoPhi());

  if (mix_pair.first == mix_pair.second) {
    // there are no possible phi_matches; return nothing
    return { {}, {} };
  }
  
  // There are tracks within the outer_box phi range;
  // now re-shuffle the phi range for eta and see if there
  // are tracks in outer_box eta range
  std::sort(mix_pair.first, mix_pair.second, CompRECOtoEta());

  RECO_pair_iter outer_box = mix_pair;
  outer_box.first = std::lower_bound(mix_pair.first, mix_pair.second,
      truth_eta-m_cutoff_deta, CompRECOtoEta());

  outer_box.second = std::lower_bound(outer_box.first, outer_box.second,
      truth_eta+m_cutoff_deta, CompRECOtoEta());

  if (outer_box.first == outer_box.second) {
    // there are no eta_matches in outerbox; resort mix_pair and return nothing
    std::sort(mix_pair.first, mix_pair.second, CompRECOtoPhi());
    return { {}, {} };
  }

  // if pt is specified, restrict pt_box to the given range
  RECO_pair_iter inner_box = outer_box;
  const float _delta_outer_pt = delta_outer_pt(truth_pt);
  const float _delta_inner_pt = delta_inner_pt(truth_pt);
  if (_delta_outer_pt > 0) {
    std::sort(outer_box.first, outer_box.second, CompRECOtoPt());
    outer_box.first  = std::lower_bound(outer_box.first, outer_box.second, truth_pt-_delta_outer_pt, CompRECOtoPt());
    outer_box.second = std::upper_bound(outer_box.first, outer_box.second, truth_pt+_delta_outer_pt, CompRECOtoPt());

    // if not outer_box, resort and return nothing
    if (outer_box.first == outer_box.second) {
      // there are no eta_matches in outerbox; resort mix_pair and return nothing
      std::sort(mix_pair.first, mix_pair.second, CompRECOtoPhi());
      return { {}, {} };
    }

    inner_box = outer_box;
    if (_delta_inner_pt > 0) {
      // start the inner pair -- the outerbox is already sorted by pT
      inner_box.first  = std::lower_bound(inner_box.first, inner_box.second, truth_pt-_delta_inner_pt, CompRECOtoPt());
      inner_box.second = std::upper_bound(inner_box.first, inner_box.second, truth_pt+_delta_inner_pt, CompRECOtoPt());
    }

    // go back to sorted by eta
    std::sort(inner_box.first, inner_box.second, CompRECOtoEta());
  }

  // At this point we know that we have outer_box matches
  // Finish checking if there are any possible inner_box matches

  // check for inner_box eta matches
  if (inner_box.first != inner_box.second) {
    inner_box.first = std::lower_bound(inner_box.first, inner_box.second,
        truth_eta-m_same_deta, CompRECOtoEta());

    inner_box.second = std::lower_bound(inner_box.first, inner_box.second,
        truth_eta+m_same_deta, CompRECOtoEta());
  }

  // check for inner_box phi matches
  if (inner_box.first != inner_box.second) {
    std::sort(inner_box.first, inner_box.second, CompRECOtoPhi());

    inner_box.first = std::lower_bound(inner_box.first, inner_box.second,
        truth_phi-m_same_dphi, CompRECOtoPhi());

    inner_box.second = std::lower_bound(inner_box.first, inner_box.second,
        truth_phi+m_same_dphi, CompRECOtoPhi());
  }

  // At this point there are definitely outer_box matches, and
  // possible inner_box matches. 
  // Return all these possible matches, evaluating all inner_box_matches
  std::vector<unsigned short> inner_box_matches, outer_box_matches;

  // fill inner_box_matches
  for (auto iter = inner_box.first; iter != inner_box.second; ++iter) {
    inner_box_matches.push_back( get<RECOid>(*iter) );
  }
  // fill inner_box_matches
  for (auto iter = outer_box.first; iter != inner_box.first; ++iter) {
    outer_box_matches.push_back(get<RECOid>(*iter));
  }
  for (auto iter = inner_box.second; iter != outer_box.second; ++iter) {
    outer_box_matches.push_back(get<RECOid>(*iter));
  }
  
  // resort recoData back to phi order
  std::sort(mix_pair.first, mix_pair.second, CompRECOtoPhi());

  // return the box matches
  return std::make_pair(inner_box_matches, outer_box_matches);
}

unsigned int TruthRecoTrackMatching::match_tracks_in_box(
    std::vector<std::pair<unsigned short,unsigned short>>& indices,  // pairs of {id_true, id_reco}
    std::set<int>& true_matched,
    std::set<int>& reco_matched,
    bool update_matched_sets)
{

  if (false) { 
    cout << "---Verbosity() " << Verbosity() << "---in ::match_tracks_in_box" << endl
         << " True : ";
    for (auto id : true_matched) cout << " " << id;
    cout << endl;
    cout << " Reco : ";
    for (auto id : reco_matched) cout << " " << id;
    cout << endl << endl;
    cout << " Pairs: ";
    for (auto dat : indices) cout << " ("<<dat.first<<","<<dat.second<<")";
    cout << endl;
  }

  unsigned int cnt_matched_tracks {0};
  // sort by truth_id, then by reco_id
  std::sort(indices.begin(), indices.end());
  auto index0 = indices.begin();

  vector<PossibleMatch> poss_matches;
  while (index0 != indices.end()) {
    auto id_true = index0->first;

    // skip all entries for which we already have a truth track match in the track
    if (true_matched.count(id_true) != 0) {
      // skip all entries for this track
      while (index0 != indices.end()) {
        ++index0;
        if (index0->first != id_true) break;
      }
      continue;
    }

    // now make all possible reco matches matched for this track
    // fill in the keys_truth to match to all keys_reco (found in sub-function)
    array<TrkrDefs::cluskey,55> keys_truth; // only generate these once for comparisons
    bool need_keys_truth = true;

    while (index0 != indices.end()) {
      if (index0->first != id_true) break;
      if (reco_matched.count(index0->second) == 0) {
        /* if (true) cout << " aa2(id_true) : " << index0->first << " (id_reco) " << index0->second << endl; */
        if (need_keys_truth) {
          /* cout << " getting keys" << endl; // FIXME */
          need_keys_truth = false;
          keys_truth = truekey_arr55(id_true);
        }
         if (false){
           cout << " ---------------------------------------" << endl
                  << " ------- start cluster matching ------- " << endl;
         }
        
        PossibleMatch match = make_PossibleMatch(index0->first, index0->second, keys_truth);
        if (false) { // FIXME
          cout << " Got possible match (nMatch nTrue nReco idTrue idReco) : "; // FIXME
          for (auto m : match) {cout << " " << m;};
          auto truth = m_TrkrTruthTrackContainer->getTruthTracks()[match[PM_idtrue]];
          auto reco = m_SvtxTrackMap->get(index0->second);
          cout << endl << " Deltas: phi(" << truth->getPhi() <<":"<<reco->get_phi() 
            << ") eta("<< truth->getPseudoRapidity()<<":"<<reco->get_eta()
            << ") pT(" << truth->getPt()<<":"<<reco->get_pt()<<")"<<endl;
        }
        if (match_pass_cuts(match)) poss_matches.push_back(match);
      }
      ++index0;
    }
    continue;
  }

  // now all the possible matches are present, already excluded any that include pairs using
  // a truth or reco track that are already matched

  // Sort them from large to small in PM_nmatch, then small to large in PM_ntrue, and PM_nreco
  if (false) {
    cout << "- - - - - - - " << endl << " FIXME : A001, pre-sort size(" << poss_matches.size() << ")" << endl;
    int i{0};
    for (auto match : poss_matches) cout << " match(" << (i++) <<") " 
      <<get<0>(match)<<"-"
      <<get<1>(match)<<"-"
      <<get<2>(match)<<"-"
      <<get<3>(match)<<"-"
      <<get<4>(match)<<endl;
    std::sort(poss_matches.begin(), poss_matches.end(), SortPossibleMatch());
    for (auto match : poss_matches) cout << " SORTED match(" << (i++) <<") " 
      <<get<0>(match)<<"-"
      <<get<1>(match)<<"-"
      <<get<2>(match)<<"-"
      <<get<3>(match)<<"-"
      <<get<4>(match)<<endl;
  }
  // ok add all possible matches started for the largest PM_nmatch (the top)
  //  for groups of ties of PM_nmatch, sort by PM_ntrue (from smallest)
  //    for groups of ties (JNmath, PM_ntrue), go by smallest PM_nreco
  //       for groups of ties (PM_nmatch, PM_ntrue, PM_nreco), do a detailed sum diff in the deltas
  std::set<int> loc_idreco, loc_idtrue;
  auto iter = poss_matches.begin();
  while (iter != poss_matches.end()) {
    if (skip_match(*iter, loc_idtrue, loc_idreco)) { 
      ++iter;
      continue;
    }
    vector<std::pair<float,int>> sigma_metric = {{0.,0}};
    int n_sigma = 0;
    auto iter_same = iter+1;
    while (
          iter_same != poss_matches.end()
       && (*iter_same) [PM_nmatch] == (*iter)[PM_nmatch]
       && (*iter_same) [PM_ntrue]  == (*iter)[PM_ntrue]
       && (*iter_same) [PM_nreco]  == (*iter)[PM_nreco]
       && !skip_match (*iter_same, loc_idtrue, loc_idreco)
    ) {
      ++n_sigma;
      if (n_sigma==1) { sigma_metric[0].first = sigma_CompCluster(*iter); }
      sigma_metric.push_back({sigma_CompCluster(*iter_same), n_sigma});
      ++iter_same;
    }
    std::sort(sigma_metric.begin(), sigma_metric.end());

    bool first = true;
    for (auto& sigma : sigma_metric) {
      if (first) first = false;
      else if (skip_match(*(iter+sigma.second), loc_idtrue, loc_idreco)) continue;
      auto match = *(iter+sigma.second);
      auto idtrue = match[PM_idtrue];
      auto idreco = match[PM_idreco];
      auto truth_track = m_TrkrTruthTrackContainer->getTruthTracks()[idtrue];
      auto reco_track  = m_SvtxTrackMap->get(idreco);
      auto save_match = new EmbRecoMatchv1( truth_track->getTrackid(), reco_track->get_id(),
          match[PM_nmatch], match[PM_ntrue], match[PM_nreco]); // FIXME -- add the id_trackseed and id_svtxtrackseed info
      m_EmbRecoMatchContainer->addMatch(save_match);
      loc_idtrue.insert(idtrue);
      loc_idreco.insert(idreco);
      ++cnt_matched_tracks;
    }
    iter += sigma_metric.size();
  }
  if (update_matched_sets) {
    for (auto& id : loc_idtrue) { true_matched.insert(id); }
    for (auto& id : loc_idreco) { true_matched.insert(id); }
  }
  return cnt_matched_tracks;
}


// ----------------------------------------
// functions for permissible matching pt
// ----------------------------------------
float TruthRecoTrackMatching::delta_outer_pt(float pt) const {
  return 10. + 0.01*pt;
}
float TruthRecoTrackMatching::delta_inner_pt(float pt) const {
  return 10. + 0.01*pt;
}

// ----------------------------------------
// convenience function
// ----------------------------------------
inline float TruthRecoTrackMatching::abs_dphi (float phi0, float phi1) {
  float dphi = fabs(phi0-phi1);
  while (dphi > M_PI) dphi -= 2*M_PI;
  return fabs(dphi);
}

// -------------------------------------------------------------------
// functions for gettings sets of TrkrClusters associated with tracks
// -------------------------------------------------------------------
array<TrkrDefs::cluskey,55> TruthRecoTrackMatching::truekey_arr55(int id_true) {
    array<TrkrDefs::cluskey,55> keys_truth; 
    for (auto& key : keys_truth) key = 0;

  auto truth_track = m_TrkrTruthTrackContainer->getTruthTracks()[id_true];
  bool pt_gt5 = (truth_track->getPt()<2 && truth_track->getPt()>0.5);
  std::ofstream outfile;
  if (pt_gt5) { //FIXME
    outfile.open("truth_phi.txt",std::ios_base::app);
    outfile << " Truth Track phi:eta:pT : " << truth_track->getPhi()
                                           <<":"<< truth_track->getPseudoRapidity() 
                                           <<":"<< truth_track->getPt() << endl;
  }
    
    for (auto& key : truth_track->getClusters()) {
      auto layerid = static_cast<int>(TrkrDefs::getLayer(key));
      /* cout << " FIXME LAYER " << static_cast<int>(layerid) << endl; */
      keys_truth[layerid] = key;

      if (pt_gt5) { //FIXME
        auto cluster = m_TruthClusterContainer->findCluster(key); 
        outfile  << " layer[" << layerid <<"] phi: " << cluster->getPosition(0) << endl;
      }

    }
    if (pt_gt5) outfile.close();
    return keys_truth;
}

std::array<std::vector<TrkrDefs::cluskey>,55> TruthRecoTrackMatching::recokey_arr55vec (int index_reco) {
  SvtxTrack* reco_track = m_SvtxTrackMap->get(index_reco);
  std::array<std::vector<TrkrDefs::cluskey>,55> keys_reco  {};

  bool pt_gt5 = (reco_track->get_pt()<2 && reco_track->get_pt()>0.5); // FIXME
  std::ofstream outfile;
  auto tpcseed = reco_track->get_tpc_seed();
  if (tpcseed) {
    if (pt_gt5) { //FIXME
      outfile.open("reco_phi.txt",std::ios_base::app);
      outfile << " Reco Track phi:eta:pT : " << reco_track->get_phi()
        <<":"<< reco_track->get_eta() 
        <<":"<< reco_track->get_pt() << " seed phi:eta:pt " 
        << " --" 
        <<":"<< tpcseed->get_eta() 
        << ":"<< tpcseed->get_pt() << " theta: " << tpcseed->get_theta() << endl;
    }

    // alternative way to get the keys
    for (auto key = tpcseed->begin_cluster_keys(); key!=tpcseed->end_cluster_keys(); ++key) {
      auto layerid = TrkrDefs::getLayer(*key);
      if (layerid > 54) continue;
      keys_reco[layerid].push_back(*key);

      if (pt_gt5) { //FIXME
        auto cluster = m_RecoClusterContainer->findCluster(*key); 
        outfile  << " layer[" << layerid <<"] phi: " << cluster->getPosition(0) << " " << cluster->getPosition(1) << endl;
      }
    }
  }
  if (pt_gt5) outfile.close();

  // FIXME
  /* for (int i=0;i<55;++i) { */
  /*   cout << " n entries in layer (" << i << ") = " << keys_reco[i].size() << endl; */
  /* } */

  /* for (auto key = reco_track->begin_cluster_keys(); key!=reco_track->end_cluster_keys(); ++key) { */
  /*   auto layerid = TrkrDefs::getLayer(*key); */
  /*   keys_reco[layerid].push_back(*key); */
  /* } */
  return keys_reco;
}

PossibleMatch TruthRecoTrackMatching::make_PossibleMatch(
    unsigned short index_truth, 
    unsigned short index_reco, 
    array<TrkrDefs::cluskey,55>& keys_truth) {
  // make the PossibleMatch tuple using index_reco and index_truth
  // loop through all clusters and see which are matches and which are not
  std::array<std::vector<TrkrDefs::cluskey>,55> keys_reco = recokey_arr55vec(index_reco);
  
  //cluster counters
  unsigned short truth_cnt {0};
  unsigned short reco_cnt  {0};
  unsigned short match_cnt {0};

  for (int i{0};i<55;++i) {
    if (keys_reco[i].size() != 0) ++reco_cnt;
    if (keys_truth[i]    != 0) ++truth_cnt;

    if ( keys_truth[i] != 0 ) {
      for (auto& key_reco : keys_reco[i]) {
        if (compare_clusters(keys_truth[i], key_reco, true).first) {
          ++match_cnt;
          break;
        }
      }
    } 
  }

  if (m_maketree) {
    // make output diagnostic values
    b_reco_cl_phi     .clear();
    b_reco_cl_phiwidth.clear();
    b_reco_cl_z       .clear();
    b_reco_cl_zwidth  .clear();

    b_truth_cl_phi     .clear();
    b_truth_cl_phiwidth.clear();
    b_truth_cl_z       .clear();
    b_truth_cl_zwidth  .clear();

    auto& truth = m_TrkrTruthTrackContainer->getTruthTracks()[index_truth];
    auto  reco  = m_SvtxTrackMap->get(index_reco);

    b_truth_phi = truth->getPhi();
    b_truth_eta = truth->getPseudoRapidity();
    b_truth_pt  = truth->getPt();

    b_reco_phi = reco->get_phi();
    b_reco_eta = reco->get_eta();
    b_reco_pt  = reco->get_pt();

    for (int i{0};i<55;++i) {
      if ( keys_truth[i] != 0 ) {
        auto cl_true  = m_TruthClusterContainer->findCluster(keys_truth[i]);
        b_truth_cl_phi      .push_back( cl_true->getPosition(0));
        b_truth_cl_z        .push_back( cl_true->getPosition(1));
        b_truth_cl_phiwidth .push_back( cl_true->getPhiSize()*pow(m_phistep_2[i],0.5));
        b_truth_cl_zwidth   .push_back( cl_true->getZSize() * pow(m_zstep_2, 0.5));
      }

      for (auto& cl_key : keys_reco[i]) {
        auto cl_reco = m_RecoClusterContainer->findCluster(cl_key);
        b_reco_cl_phi      .push_back( cl_reco->getPosition(0));
        b_reco_cl_z        .push_back( cl_reco->getPosition(1));
        b_reco_cl_phiwidth .push_back( cl_reco->getPhiSize()*pow(m_phistep_2[i],0.5));
        b_reco_cl_zwidth   .push_back( cl_reco->getZSize() * pow(m_zstep_2, 0.5));
      }
    }
    m_fileout->cd();
    m_treeout->Fill();
  }
  /* return std::make_tuple(match_cnt, truth_cnt, reco_cnt, index_truth, index_reco); */
  return { match_cnt, truth_cnt, reco_cnt, index_truth, index_reco };
}

std::pair<bool, double> TruthRecoTrackMatching::compare_clusters(
    TrkrDefs::cluskey truthkey, TrkrDefs::cluskey recokey, bool print) {
  // return phi_reco-phi_true, phi_size_true/step, phi_size_reco/step, 
  //        ditto for Z
  if (TrkrDefs::getHitSetKeyFromClusKey(truthkey) != TrkrDefs::getHitSetKeyFromClusKey(recokey)) return {false, 0.};

  /* cout << " <<<<< Same hitset key? layer(" << ((int)TrkrDefs::getLayer(truthkey))<<":"<<((int)TrkrDefs::getLayer(recokey)) */ 
    /* << ( (TrkrDefs::getHitSetKeyFromClusKey(truthkey) == TrkrDefs::getHitSetKeyFromClusKey(recokey)) ? " yes " : " NO ") << endl; */

  auto cl_true  = m_TruthClusterContainer->findCluster(truthkey);
  auto cl_reco  = m_RecoClusterContainer-> findCluster(recokey);

  //NOTE: criteria is somewhat arbitrary
  //      Starting pick is to require that the location reconstructed be within the width of the 
  //      truth level
  auto layerid = TrkrDefs::getLayer(truthkey);

  const float phi_true      = cl_true->getPosition(0);
  const float phi_reco      = cl_reco->getPosition(0);
  const float dphi_2        = pow( abs_dphi(phi_reco, phi_true) , 2. );

  const float phi_size_true = cl_true->getPhiSize();
  const float phi_size_reco = cl_reco->getPhiSize();
  const float phi_spread_2  = (pow(phi_size_true,2.) + pow(phi_size_reco,2.)) * m_phistep_2[layerid];

  const float z_true        = cl_true ->getPosition(1);
  const float z_reco        = cl_reco ->getPosition(1);
  const float dz_2          = pow( (z_reco-z_true), 2.);

  const float z_size_true   = cl_true ->getZSize();
  const float z_size_reco   = cl_reco->getZSize();
  const float z_spread_2    = (pow(z_size_true,2.) + pow(z_size_reco,2.)) * m_zstep_2;

  if (print && false) {
    std::cout << "phi_true("<<phi_true<<") phi_reco("<<phi_reco<<")   -> testphi(" << 
      dphi_2 << " < " << (phi_spread_2*m_cluster_nphiwidths_2) << ") ? phi:" <<
    ( dphi_2      <      (phi_spread_2*m_cluster_nphiwidths_2) ) << std::endl;
    std::cout << "phi more : phi_true("<<phi_true<<") phi_reco("<<phi_reco <<") " 
      << " phi_delta_t("<< phi_size_true <<")  phi_size_reco(" << phi_size_reco<<") " 
      << " --step(" << pow(m_phistep_2[layerid],0.5) <<") " << std::endl;
    std::cout << "z_true("<<z_true<<") z_reco("<<z_reco<<")   -> testZ(" << 
      dz_2 << " < " << (z_spread_2*m_cluster_nzwidths_2) << ") ? z:" <<
    ( dz_2      <      (z_spread_2*m_cluster_nzwidths_2) ) << std::endl;
    std::cout << "z more : z_true("<<z_true<<") z_reco("<<z_reco <<") " 
      << " z_delta_t("<< z_size_true <<")  z_size_reco(" << z_size_reco<<") " << std::endl;
  }

  bool is_match = (dphi_2 < (phi_spread_2*m_cluster_nphiwidths_2) 
               && (dz_2   < (z_spread_2  *m_cluster_nzwidths_2)));
  if (!is_match) return { false, 0. };
  else           return { true,  dphi_2/phi_spread_2 + dz_2/z_spread_2 };
}

float TruthRecoTrackMatching::sigma_CompCluster(PossibleMatch& match) {
  auto index_true = match[PM_idtrue];
  auto index_reco = match[PM_idreco];

  auto keys_true = truekey_arr55(index_true);
  auto keys_reco = recokey_arr55vec(index_reco);

  // output metric: mss_phierror + mss_zerror (Mean-Sum-Square)
  int match_cnt {0};
  float mss_sum {0.};
  for (int i{0};i<55;++i) {
    if ( keys_true[i] != 0 ) {
      for (auto& key_reco : keys_reco[i]) {
        std::pair<bool, double> comp_arr = compare_clusters(keys_true[i],key_reco);
        if (comp_arr.first) {
          ++match_cnt;
          mss_sum += comp_arr.second;
          break;
        }
      }
    } 
  }
  if (match_cnt > 0) {
    return mss_sum / match_cnt;
  } else {
    std::cout << PHWHERE << " error in calling sigma_CompCluster: no matched clusters " << std::endl;
    return -1.;
  }
}

inline bool TruthRecoTrackMatching::skip_match(
    PossibleMatch& match,
    std::set<int>& true_matched,
    std::set<int>& reco_matched) 
{
  return 
      (true_matched.count(match[PM_idtrue]) != 0)
  ||  (reco_matched.count(match[PM_idreco]) != 0);
}

inline bool TruthRecoTrackMatching::match_pass_cuts(PossibleMatch& match) {
  return match[PM_nmatch] >= m_nmincluster_match
     &&  static_cast<float>(match[PM_nmatch]) / match[PM_ntrue] >= m_nmincluster_ratio  ;
}
