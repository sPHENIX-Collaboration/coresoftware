#ifndef TRACKBASE_EMBRECOMATCHCONTAINER_H
#define TRACKBASE_EMBRECOMATCHCONTAINER_H

/* #include "TrkrDefs.h" */

#include <phool/PHObject.h>
#include <vector>

class EmbRecoMatch;

class EmbRecoMatchContainer : public PHObject
{
  public:
  using Vector        = std::vector<EmbRecoMatch*>;
  using ConstIterator = Vector::const_iterator;
  using ConstRange    = std::pair<ConstIterator, ConstIterator>;
  /* using Iterator      = Vector::iterator; */
  /* using Range         = std::pair<Iterator, Iterator>; */

  //! reset method
  void Reset() override {}

  //! add a match set
  virtual unsigned short nMatches() const  = 0;
  virtual unsigned short nUnmatched() const = 0;
  virtual std::vector<unsigned short>& ids_Unmatched() = 0; // id's of the TrkrTruthTrack's that are not matched
  virtual std::vector<unsigned short> ids_Matched() const = 0;   // id's of the TrkrTruthTrack's that are matched
  /* virtual Range           getMatches() = 0; */
  virtual ConstRange getMatchedRange() const = 0;
  virtual Vector&    getMatches()     = 0;

  virtual void addMatch(EmbRecoMatch*) = 0;
  virtual bool hasMatch(unsigned short idEmb) = 0;
  virtual EmbRecoMatch* getMatch(unsigned short idEmb) = 0;

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
