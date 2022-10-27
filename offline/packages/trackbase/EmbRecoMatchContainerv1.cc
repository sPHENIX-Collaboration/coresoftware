#include "EmbRecoMatchContainerv1.h"
#include "EmbRecoMatch.h"
#include "EmbRecoMatchv1.h"

#include <algorithm>
#include <iomanip>

void EmbRecoMatchContainerv1::Reset() {
  m_data.clear();
  m_ids_Unmatched.clear();
}

std::vector<unsigned short>  EmbRecoMatchContainerv1::ids_Matched() const {
  std::vector<unsigned short> vec;
  for (auto& id : m_data) vec.push_back(id->idTruthTrack());
  return vec;
}

EmbRecoMatch* EmbRecoMatchContainerv1::getMatch(unsigned short idEmb) {
  auto iter = std::lower_bound(m_data.begin(), m_data.end(), idEmb, EmbRecoMatch::Comp());
  if (iter == m_data.end()) {
    std::cout << "Asking for EmbTrackMatch for Embedded Track id " << idEmb 
              << " which is not present. Returning empty match." << std::endl;
    EmbRecoMatch* track = new EmbRecoMatchv1();
    return track;
  }
  return *iter;
}

bool EmbRecoMatchContainerv1::hasMatch(unsigned short idEmb) {
  return std::binary_search(m_data.begin(), m_data.end(), idEmb, EmbRecoMatch::Comp());
}

void EmbRecoMatchContainerv1::identify(std::ostream& os) const {
  os << " EmbRecoMatchContainerv1 data. N(matched emb. tracks) = " 
     << nMatches() << ",  N(not matched) " << nUnmatched() << std::endl;
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
     << std::setw(7)  << m->idTrackSeed()       << " " //"id-Seed"
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
  for (const auto n : m_ids_Unmatched) os << n << " ";
  os << std::endl;
}
