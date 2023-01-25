#include "EmbRecoMatchContainerv1.h"
#include "EmbRecoMatch.h"
#include "EmbRecoMatchv1.h"

#include <algorithm>
#include <iomanip>

// WARNING to user:
// methods depend on the internal data being sorted appropriately
// -> BE SURE TO CALL *SORT*

void EmbRecoMatchContainerv1::sort() {
  std::sort(m_data.begin(),              m_data.end(), EmbRecoMatch::Comp());
  std::sort(m_RecoToTruth.begin(),       m_RecoToTruth.end());
  std::sort(m_idsTruthUnmatched.begin(), m_idsTruthUnmatched.end());
}

void EmbRecoMatchContainerv1::addMatch(EmbRecoMatch* match) {
  m_data.push_back(match);
  auto id_true = match->idTruthTrack();
  auto id_reco = match->idRecoTrack();

  m_RecoToTruth.push_back({id_reco, id_true}); // vector of which to go to

  if (m_nTruthPerReco.find(id_reco) == m_nTruthPerReco.end()) {
    m_nTruthPerReco[id_reco] = 1;
  } else {
    m_nTruthPerReco[id_reco] += 1;
  }
  if (m_nRecoPerTruth.find(id_true) == m_nRecoPerTruth.end()) {
    m_nRecoPerTruth[id_true] = 1;
  } else {
    m_nRecoPerTruth[id_true] += 1;
  }
}

void EmbRecoMatchContainerv1::checkfill_idsTruthUnmatched(unsigned short id_true) {
  if (m_nRecoPerTruth.find(id_true) == m_nRecoPerTruth.end()) {
    m_idsTruthUnmatched.push_back(id_true);
  }
}

void EmbRecoMatchContainerv1::Reset() {
  for (auto &m : m_data) delete m;
  m_data              .clear();
  m_RecoToTruth       .clear();
  m_idsTruthUnmatched .clear();
  m_nTruthPerReco     .clear();
  m_nRecoPerTruth     .clear();
}

std::vector<unsigned short>  EmbRecoMatchContainerv1::ids_TruthMatched() const {
  std::vector<unsigned short> vec;
  for (auto& id : m_data) vec.push_back(id->idTruthTrack());
  return vec;
}

std::vector<unsigned short>  EmbRecoMatchContainerv1::ids_RecoMatched() const {
  std::vector<unsigned short> vec;
  for (auto& id : m_RecoToTruth) vec.push_back(id.first);
  return vec;
}

EmbRecoMatch* EmbRecoMatchContainerv1::getMatchTruth(unsigned short idtruth, unsigned short offset) {
  auto iter = std::lower_bound(m_data.begin(), m_data.end(), idtruth, EmbRecoMatch::Comp());
  iter += offset;

  if (iter >= m_data.end() || (*iter)->idTruthTrack() != idtruth) {
    std::cout << "Error: asking for match (offset by " << offset <<") for truth track id " << idtruth 
      << " which is not present. Returning null pointer." << nullptr;
    return nullptr;
  } 
  return *iter;
}

EmbRecoMatch* EmbRecoMatchContainerv1::getMatchReco(unsigned short idreco, unsigned short offset) {
  auto iter = std::lower_bound(m_RecoToTruth.begin(), m_RecoToTruth.end(), idreco, CompShortToPair());
  iter += offset;
  if (iter >= m_RecoToTruth.end() || iter->first != idreco) {
    std::cout << "Error: asking for match (offset by " << offset <<") for reco track id " << idreco
      << " which is not present. Returning null pointer." << nullptr;
    return nullptr;
  } 
  return getMatchTruth(iter->second, offset);
}

bool EmbRecoMatchContainerv1::hasTruthMatch(unsigned short idtruth) {
  return std::binary_search(m_data.begin(), m_data.end(), idtruth, EmbRecoMatch::Comp());
}

bool EmbRecoMatchContainerv1::hasRecoMatch(unsigned short idreco) {
  return std::binary_search(m_RecoToTruth.begin(), m_RecoToTruth.end(), idreco, CompShortToPair());
}

void EmbRecoMatchContainerv1::identify(std::ostream& os) const {
  os << " EmbRecoMatchContainerv1 data. N(matched emb. tracks) = " 
     << nMatches() << ",  N(not matched) " << nTruthUnmatched() << std::endl;
  os << " Matched track data: " << std::endl;
  //    id-Emb  id-Reco  id-Seed  id-SvtxSeed  nClus-Emb nClus-Reco meanZClus meanPhiClus
  os << std::setw(6) << "id-Emb " 
     << std::setw(7) << "id-Reco " 
     << std::setw(7) << "id-Seed "
     << std::setw(11) << "id-SvtxSeed "
     << std::setw(9)  << "nClus-Emb "
     << std::setw(10) << "nClus-Reco "
     << std::setw(12) << "nClus_Matched " << std::endl;
  for (auto& _m : m_data) {
    auto m = static_cast<EmbRecoMatchv1*> (_m);
  os << std::setw(6)  << m->idTruthTrack()      << " " //id-Emb" 
     << std::setw(7)  << m->idRecoTrack()       << " " //"id-Reco"
     << std::setw(7)  << m->idTpcTrackSeed()    << " " //"id-Seed"
     << std::setw(11) << m->idSvtxTrackSeed()   << " " //id-SvtxSeed"
     << std::setw(9)  << m->nClustersTruth()    << " " //nClus-Emb"
     << std::setw(10) << m->nClustersReco()     << " " //"nClus-Reco"
     << std::setw(12) << m->nClustersMatched()  << std::endl; // //"nClus-Reco"
     /* << std::setw(9)  << m->meanClusterZDiff() << " " //"meanZclus" */
     /* << std::setw(11) << m->meanClusterPhiDiff() << std::endl; //"meanPhiclus" << std::endl; */

  /* os << Form(" %7i %7i %8s %11s %11s %10s %10s %10s", */
    /* m->idTruthTrack(), m->idRecoTrack(), m->idTrackSeed(), m->idSvtxTrackSeed(), */
    /* m->nClustersTruth(), m->nClustersReco(), m->meanClusterZDiff(), m->meanClusterPhiDiff()) << std::endl; */
  }
  os << " IDs of embedded tracks that were not reconstructed: " << std::endl;
  for (const auto n : m_idsTruthUnmatched) os << n << " ";
  os << std::endl;
}
