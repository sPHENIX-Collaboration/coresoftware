#ifndef TRACKBASE_EMBRECOMATCHCONTAINER_H
#define TRACKBASE_EMBRECOMATCHCONTAINER_H

/* #include "TrkrDefs.h" */

#include <phool/PHObject.h>
#include <vector>
#include <climits>
#include <map>

class EmbRecoMatch;

class EmbRecoMatchContainer : public PHObject
{
  public:
  using Vector        = std::vector<EmbRecoMatch*>;
  using ConstIterator = Vector::const_iterator;
  using ConstRange    = std::pair<ConstIterator, ConstIterator>;
  using Map_nMatches  = std::map<unsigned short, unsigned short>;
  /* using Iterator      = Vector::iterator; */
  /* using Range         = std::pair<Iterator, Iterator>; */

  //! reset method
  void Reset() override {}

  //! add a match set
  virtual unsigned short nMatches()   const { return USHRT_MAX; };
  virtual unsigned short nTruthUnmatched() const { return USHRT_MAX; };

  virtual std::vector<unsigned short>& ids_TruthUnmatched(); // id's of the TrkrTruthTrack's that are not matched
  virtual std::vector<unsigned short>  ids_TruthMatched() const { return {}; }; // id's of the TrkrTruthTrack's that are matched
  virtual std::vector<unsigned short>  ids_RecoMatched()  const { return {}; }; // id's of the TrkrTruthTrack's that are matched

  virtual std::map<unsigned short, unsigned short>& map_nRecoPerTruth();
  virtual std::map<unsigned short, unsigned short>& map_nTruthPerReco();

  virtual ConstRange getMatchedRange() const;
  virtual Vector&    getMatches()           ;

  virtual void addMatch(EmbRecoMatch*) {};
  virtual bool hasTruthMatch (unsigned short /*Truth Track Id*/) { return false; };
  virtual bool hasRecoMatch  (unsigned short /*Reco  Track Id*/) { return false; };
  virtual EmbRecoMatch* getMatchTruth (unsigned short /*idEmb*/,  unsigned short /*offset*/=0) { return nullptr; };
  virtual EmbRecoMatch* getMatchReco  (unsigned short /*idReco*/, unsigned short /*offset*/=0) { return nullptr; };

  virtual void checkfill_idsTruthUnmatched(unsigned short){}; // add tracks one-by-one and add to m_idsTruthUnmatched if not matched

  virtual void sort() {}; // only to be used when filling the contianer. Sorts the internal vectors for the sake of using them.

  // PHObject virtual overload
  void identify(std::ostream& os = std::cout) const override
  {
    os << "EmbREcoMatchContainer base class" << std::endl;
  };

 protected:
  EmbRecoMatchContainer() = default;
  ClassDefOverride(EmbRecoMatchContainer, 1)
};

#endif // TRACKBASE_EMBRECOMATCHCONTAINER_H
