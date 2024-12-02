#ifndef TRACKBASE_EMBRECOMATCHV1_H
#define TRACKBASE_EMBRECOMATCHV1_H

#include "EmbRecoMatch.h"

#include <array>
#include <iostream>
#include <limits>
#include <vector>

/* class VtxPoint; */
class EmbRecoMatchv1 : public EmbRecoMatch
{
 public:
  ~EmbRecoMatchv1() override = default;

  unsigned short idTruthTrack() const override { return m_idTruthTrack; };
  unsigned short nClustersTruth() const override { return m_nClustersTruth; };
  unsigned short nMatches() const override { return m_nMatches; };

 private:
  static constexpr int MATCH_idReco = 0;
  static constexpr int MATCH_idTpcSeed = 1;
  static constexpr int MATCH_idSvtxSeed = 2;
  static constexpr int MATCH_nClustersReco = 3;
  static constexpr int MATCH_nClustersMatched = 4;

 public:
  unsigned short idRecoTrack(unsigned short i = 0) const override { return m_matches[i][MATCH_idReco]; };
  unsigned short idTpcTrackSeed(unsigned short i = 0) const override { return m_matches[i][MATCH_idTpcSeed]; };
  unsigned short idSvtxTrackSeed(unsigned short i = 0) const override { return m_matches[i][MATCH_idSvtxSeed]; };

  unsigned short nClustersReco(unsigned short i = 0) const override { return m_matches[i][MATCH_nClustersReco]; };
  unsigned short nClustersMatched(unsigned short i = 0) const override { return m_matches[i][MATCH_nClustersMatched]; };

  /* float meanClusterZDiff()   const override { return m_meanClusterZDiff; }; */
  /* float meanClusterPhiDiff() const override { return m_meanClusterPhiDiff; }; */

  EmbRecoMatchv1(
      unsigned short id_truth = std::numeric_limits<unsigned short>::max(), unsigned short id_reco = std::numeric_limits<unsigned short>::max(), unsigned short nclustruth = std::numeric_limits<unsigned short>::max(), unsigned short nclusreco = std::numeric_limits<unsigned short>::max(), unsigned short nclusmatched = std::numeric_limits<unsigned short>::max(), unsigned short id_tpctrackseed = std::numeric_limits<unsigned short>::max(), unsigned short id_svtxtrackseed = std::numeric_limits<unsigned short>::max())
    : m_idTruthTrack{id_truth}
    , m_nClustersTruth{nclustruth}
    , m_nMatches{1}
    , m_matches{{id_reco, id_tpctrackseed, id_svtxtrackseed, nclusreco, nclusmatched}} {};

  unsigned short add_match(
      unsigned short id_reco = std::numeric_limits<unsigned short>::max(), unsigned short nclusreco = std::numeric_limits<unsigned short>::max(), unsigned short nclusmatched = std::numeric_limits<unsigned short>::max(), unsigned short id_tpctrackseed = std::numeric_limits<unsigned short>::max(), unsigned short id_svtxtrackseed = std::numeric_limits<unsigned short>::max());

 private:
  unsigned short m_idTruthTrack{std::numeric_limits<unsigned short>::max()};
  unsigned short m_nClustersTruth{std::numeric_limits<unsigned short>::max()};
  unsigned short m_nMatches{std::numeric_limits<unsigned short>::max()};

  std::vector<std::array<unsigned short, 5>> m_matches;
  // each match has metrics defined by five unsigned ints: idReco, idTpcSee, idSvtxSeed, nClustersReco, nClusterMatched

 public:
  // PHObject virtual overloads
  void identify(std::ostream& os = std::cout) const override
  {
    os << "EmbRecoMatchv1 base class" << std::endl;
  };

 protected:
  ClassDefOverride(EmbRecoMatchv1, 1)
};

#endif  // G4TPC_TruthTrack_h
