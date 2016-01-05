#ifndef PHG4INEVENTREADBACK_H
#define PHG4INEVENTREADBACK_H

#include <fun4all/SubsysReco.h>

class PHG4InEventReadBack: public SubsysReco
{
 public:
  PHG4InEventReadBack(const std::string &name = "PHG4InEventReadBack");
  virtual ~PHG4InEventReadBack() {}
  int InitRun(PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode);
  int End(PHCompositeNode *topNode);

 protected:
  
};

#endif
