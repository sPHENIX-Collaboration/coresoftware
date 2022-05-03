#ifndef TRACKRECO_ASSOCINFOCONTAINER_H
#define TRACKRECO_ASSOCINFOCONTAINER_H

#include <phool/PHObject.h>

#include <trackbase/TrkrDefs.h>

#include <iostream>              // for cout, ostream
#include <map>
#include <utility>               // for pair
#include <vector>                // for vector

class AssocInfoContainer : public PHObject
{
 public:
  typedef std::multimap<TrkrDefs::cluskey, unsigned int> ClusterTrackMap;

  ~AssocInfoContainer() override{}

  void identify(std::ostream& os = std::cout) const override;

  virtual void SetClusterTrackAssoc(const TrkrDefs::cluskey& /*cluster_id*/, const unsigned int& /*track_id*/) {return;}

  virtual std::vector<unsigned int> GetTracksFromCluster(const TrkrDefs::cluskey& /*cluster_id*/) const { std::vector<unsigned int> emptyvec; return emptyvec;}

 protected:
  AssocInfoContainer(){};

 private:

  ClassDefOverride(AssocInfoContainer, 1)
};

#endif
