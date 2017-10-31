#ifndef __FUN4ALLDSTOUTPUTMANAGER_H__
#define __FUN4ALLDSTOUTPUTMANAGER_H__


#include "Fun4AllOutputManager.h"
#include <string>
#include <vector>

class PHNodeIOManager;
class PHCompositeNode;

class Fun4AllDstOutputManager: public Fun4AllOutputManager
{
 public:

  Fun4AllDstOutputManager(const std::string &myname = "DSTOUT" , const std::string &filename = "dstout.root");
  virtual ~Fun4AllDstOutputManager();

  int AddNode(const std::string  &nodename);
  int StripNode(const std::string  &nodename);
  int outfileopen(const std::string &fname);
  int RemoveNode(const std::string &nodename);

  void Print(const std::string &what = "ALL") const;

  int Write(PHCompositeNode *startNode);
  int WriteNode(PHCompositeNode *thisNode, const PHTreeType);

 protected:
  std::vector <std::string> savenodes;
  std::vector <std::string> stripnodes;
  PHNodeIOManager *dstOut;
};

#endif /* __FUN4ALLDSTOUTPUTMANAGER_H__ */
