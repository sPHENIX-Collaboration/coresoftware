#ifndef TRACKBASE_EMBRECOMATCHCONTAINERV1_H
#define TRACKBASE_EMBRECOMATCHCONTAINERV1_H

/* #include "TrkrDefs.h" */
#include "EmbRecoMatchContainer.h"

/* #include <phool/PHObject.h> */
#include <vector>
/* #include <map> */

class EmbRecoMatch;

class EmbRecoMatchContainerv1 : public EmbRecoMatchContainer 
{
  public:
  void Reset() override;

  //! add a match set
  unsigned short          nMatches()       const override { return m_data.size();          };
  unsigned short          nUnmatched()     const override { return m_ids_Unmatched.size(); };
  std::vector<unsigned short>& ids_Unmatched()  override { return m_ids_Unmatched; };
  std::vector<unsigned short>  ids_Matched() const   override;
  ConstRange              getMatchedRange() const override { return std::make_pair(m_data.begin(), m_data.end()); };
  Vector&                 getMatches()      override { return m_data; };

  void addMatch(EmbRecoMatch* match)        override { m_data.push_back(match); };
  bool hasMatch(unsigned short idEmb) override;
  EmbRecoMatch* getMatch(unsigned short idEmb) override;

  // PHObject virtual overload
  void identify(std::ostream& os = std::cout) const override;

  private:
  Vector m_data;
  std::vector<unsigned short> m_ids_Unmatched;
  ClassDefOverride(EmbRecoMatchContainerv1, 1)
};

#endif // TRACKBASE_EMBRECOMATCHCONTAINERV1_H

