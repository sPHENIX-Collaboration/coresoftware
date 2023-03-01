#include "TrkrTruthTrackContainerv1.h"
#include "TrkrTruthTrackv1.h"

#include <algorithm>
#include <g4main/PHG4Particle.h>
#include <g4main/PHG4TruthInfoContainer.h>
#include <g4main/PHG4VtxPoint.h>
#include <iostream>
#include <phool/phool.h>
#include <set>

using std::endl;
using std::cout;
using std::set;

void TrkrTruthTrackContainerv1::Reset()
{
  for (auto entry : m_data)
  {
    delete entry.second;
  }
  m_data.clear();
}

void TrkrTruthTrackContainerv1::addTruthTrack(TrkrTruthTrack* track)
{
  if (!track) return;
  const auto id { track->getTrackid() };
  if (hasTrackid(id)) {
    cout << PHWHERE << "Warning, replacing existing track-id("<<id<<")" << endl;
  }
  m_data[id] = track;
}

TrkrTruthTrack* TrkrTruthTrackContainerv1::getTruthTrack(unsigned int trackid)
{
  if (!hasTrackid(trackid)) {
    cout << PHWHERE << " Asking for TrkrTruthTrack " << trackid 
      << " which is not present. Returning empty track." << endl;
    TrkrTruthTrack* track = new TrkrTruthTrackv1();
    return track;
  }
  return m_data[trackid];
}

TrkrTruthTrack* TrkrTruthTrackContainerv1::getTruthTrack (unsigned int id, PHG4TruthInfoContainer* truth_info) { 
  // return the track if already in m_data, otherwise make it and return the newly made track
  if (hasTrackid(id)) return m_data[id];
      PHG4Particle* particle = /*(PHG4Particlev3*)*/ truth_info->GetParticle(id);
      if (particle == nullptr) {
        cout << PHWHERE << " Note: embedded track from PHG4TruthInfoContainer, id( " 
          << id <<" )without an associated PHG4Particle" << endl;
        auto current_track = new TrkrTruthTrackv1();
        current_track->setTrackid(id);
        m_data[id] = current_track;
        return current_track;
      }
      int vtxid = particle->get_vtx_id();
      PHG4VtxPoint* vtx = truth_info->GetVtx(vtxid);
      auto current_track = new TrkrTruthTrackv1(id, particle, vtx) ;
      m_data[id] = current_track;
      return current_track;
}


/* void TrkrTruthTrackContainerv1::addAllTruthTracks(PHG4TruthInfoContainer* truth_info, bool only_if_empty) */
/* { */
/*   bool is_empty = (m_data.size()==0); */
/*   if (only_if_empty && !is_empty) { */
/*     cout << PHWHERE << "Warning, ingnoring call TrkrTruthTrackContainverv1::addAllTruthTracks, because it already has some tracks" << endl; */
/*     cout << " -> therefore ignoring call" << endl; */
/*     cout << " If you really wish to do this, add `true` as the final argument to addAllTruthTracks" << endl; */
/*     return; */
/*   } */

/*   auto iter_truth_tracks = truth_info->GetEmbeddedTrkIds(); */
/*   if (iter_truth_tracks.first == iter_truth_tracks.second) return; */

/*   if (is_empty) { */
/*     // just fill in all the data */
/*     for (auto it = iter_truth_tracks.first; it != iter_truth_tracks.second; ++it) { */
/*       auto trkid = it->first; */
/*       if (!truth_info->isEmbeded(trkid)) continue; */
///*       PHG4Particle* particle = /*(PHG4Particlev3*)*/ truth_info->GetParticle(trkid); */
/*       if (particle == nullptr) { */
/*         cout << PHWHERE << " Note: embedded track from PHG4TruthInfoContainer, id( " << trkid <<" )without an associated PHG4Particle" << endl; */
/*         auto current_track = new TrkrTruthTrackv1(); */
/*         current_track->setTrackid(trkid); */
/*         m_data.push_back(current_track); */
/*         continue; */
/*       } */
/*       int vtxid = particle->get_vtx_id(); */
/*       PHG4VtxPoint* vtx = truth_info->GetVtx(vtxid); */
/*       auto current_track = new TrkrTruthTrackv1(trkid, particle, vtx) ; */
/*       m_data.push_back(current_track); */
/*     } */
/*     return; */
/*   } */ 

/*   // @ this point, is_empty == false && only_if_empty == false */
/*   // -> only add new tracks that do not already exist */
/*   std::set<unsigned int> existing_ids; */
/*   for (auto& track : getTruthTracks()) { */
/*     existing_ids.insert(track->getTrackid()); */
/*   } */

/*   int cnt_added {0}; */
/*   int cnt_total {0}; */

/*   for (auto it = iter_truth_tracks.first; it != iter_truth_tracks.second; ++it) { */
/*     ++cnt_total; */
/*     auto trkid = it->first; */
/*     if (!truth_info->isEmbeded(trkid) ||  (existing_ids.count(trkid) != 0)) continue; */
/*     ++cnt_added; */
///*     PHG4Particle* particle = /*(PHG4Particlev3*)*/ truth_info->GetParticle(trkid); */
/*     if (particle == nullptr) { */
/*       cout << PHWHERE << " Note: embedded track from PHG4TruthInfoContainer, id( " << trkid <<" )without an associated PHG4Particle" << endl; */
/*       auto current_track = new TrkrTruthTrackv1(); */
/*       current_track->setTrackid(trkid); */
/*       m_data.push_back(current_track); */
/*       continue; */
/*     } */
/*     int vtxid = particle->get_vtx_id(); */
/*     PHG4VtxPoint* vtx = truth_info->GetVtx(vtxid); */
/*     auto current_track = new TrkrTruthTrackv1(trkid, particle, vtx) ; */
/*     m_data.push_back(current_track); */
/*   } */
/*   std::sort(m_data.begin(), m_data.end(), TrkrTruthTrack::Comp()); */
/*   cout << " Added " << cnt_added << " tracks in addAllTruthTracks, out of possible " << cnt_total << "." << endl; */
/* } */


TrkrTruthTrackContainer::ConstRange TrkrTruthTrackContainerv1::getTruthTrackRange() const
{
  return std::make_pair(m_data.begin(), m_data.end());
}

bool TrkrTruthTrackContainerv1::hasTrackid(unsigned int id) const
{
  return (m_data.find(id) != m_data.end());
}

TrkrTruthTrackContainer::Map& TrkrTruthTrackContainerv1::getMap()
{
  return m_data;
}


void TrkrTruthTrackContainerv1::identify(std::ostream& os) const
{
  os << " TrkrTruthTrackContainer data.  Containter " << (int)m_data.size() << " tracks" << std::endl;
  int cnt = 0;
  for (auto& entry : m_data)
  {
    os << " Track(" << cnt << "): " << std::endl;
    entry.second->identify(os);
    ++cnt;
  }
}
