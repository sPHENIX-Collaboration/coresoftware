#ifndef TRACKBASE_EMBRECOMATCHCONTAINERV1_H
#define TRACKBASE_EMBRECOMATCHCONTAINERV1_H

#include "EmbRecoMatchContainer.h"
#include <vector>

class EmbRecoMatch;

class EmbRecoMatchContainerv1 : public EmbRecoMatchContainer 
{
  public:
  void Reset() override;

  //! add a match set
  unsigned short nMatches()   const override { return m_data.size();          };
  unsigned short nTruthUnmatched() const override { return m_idsTruthUnmatched.size(); };
  std::vector<unsigned short>& ids_TruthUnmatched() override { return m_idsTruthUnmatched; };
  std::vector<unsigned short>  ids_TruthMatched() const override;
  std::vector<unsigned short>  ids_RecoMatched()  const override;

  std::map<unsigned short, unsigned short>& map_nRecoPerTruth() override { return m_nRecoPerTruth; };
  std::map<unsigned short, unsigned short>& map_nTruthPerReco() override { return m_nTruthPerReco; };

  ConstRange getMatchedRange() const override { return  std::make_pair(m_data.begin(), m_data.end()); };
  Vector&    getMatches()            override { return m_data; };

  void addMatch(EmbRecoMatch* match)       override;
  bool hasTruthMatch(unsigned short idEmb) override;
  bool hasRecoMatch(unsigned short idReco) override;
  EmbRecoMatch* getMatchTruth(unsigned short idEmb, unsigned short offset=0) override;
  EmbRecoMatch* getMatchReco (unsigned short idEmb, unsigned short offset=0) override;

  // PHObject virtual overload
  void identify(std::ostream& os = std::cout) const override;

  struct CompShortToPair {
    bool operator() (const unsigned short lhs, const std::pair<unsigned short, unsigned short> rhs) const
    {return lhs < rhs.first;}
    bool operator() (const std::pair<unsigned short, unsigned short> lhs, const unsigned short rhs) const
    {return lhs.first < rhs;}
  };

  void sort() override;
  void checkfill_idsTruthUnmatched(unsigned short) override; // add tracks one-by-one and add to m_idsTruthUnmatched if not matched

  private:
  Vector m_data; // using Vector = std::vector<EmbRecoMatch*>;
  std::vector<std::pair<unsigned short, unsigned short>> m_RecoToTruth;
  std::vector<unsigned short> m_idsTruthUnmatched;
  Map_nMatches m_nRecoPerTruth;
  Map_nMatches m_nTruthPerReco;
  ClassDefOverride(EmbRecoMatchContainerv1, 1)
};

#endif // TRACKBASE_EMBRECOMATCHCONTAINERV1_H

