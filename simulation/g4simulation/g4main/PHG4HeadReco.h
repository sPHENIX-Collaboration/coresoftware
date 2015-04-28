#ifndef PHG4HEADRECO_H__
#define PHG4HEADRECO_H__

#include <fun4all/SubsysReco.h>

class PHG4HeadReco: public SubsysReco
{
 public:
  PHG4HeadReco(const std::string &name="PHG4HeadReco");
  virtual ~PHG4HeadReco(){}
  int Init(PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode);

 protected:

  int evtseq;

};

#endif
