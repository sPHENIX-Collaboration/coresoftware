#ifndef __PHG4HITREADBACK__H__
#define __PHG4HITREADBACK__H__

#include <fun4all/SubsysReco.h>

#include <string>                // for string

class PHCompositeNode;

class PHG4HitReadBack : public SubsysReco
{
 public:
  PHG4HitReadBack(const std::string &name="PHG4HITREADBACK");
  int InitRun(PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode);
 protected:
};


#endif
