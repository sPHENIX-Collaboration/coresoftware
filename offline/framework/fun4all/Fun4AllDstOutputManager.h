#ifndef FUN4ALLDSTOUTPUTMANAGER_H__
#define FUN4ALLDSTOUTPUTMANAGER_H__

#include "Fun4AllOutputManager.h"

#include <set>
#include <string>

class PHNodeIOManager;
class PHCompositeNode;

class Fun4AllDstOutputManager : public Fun4AllOutputManager
{
 public:
  Fun4AllDstOutputManager(const std::string &myname = "DSTOUT", const std::string &filename = "dstout.root");
  virtual ~Fun4AllDstOutputManager();

  int AddNode(const std::string &nodename);
  int StripNode(const std::string &nodename);
  int StripRunNode(const std::string &nodename);
  int outfileopen(const std::string &fname);

  void Print(const std::string &what = "ALL") const;

  int Write(PHCompositeNode *startNode);
  int WriteNode(PHCompositeNode *thisNode);

 protected:
  std::set<std::string> savenodes;
  std::set<std::string> stripnodes;
  std::set<std::string> striprunnodes;
  PHNodeIOManager *dstOut;
};

#endif /* FUN4ALLDSTOUTPUTMANAGER_H__ */
