#ifndef HEPMCNODEREADER_H__
#define HEPMCNODEREADER_H__

#include <fun4all/SubsysReco.h>

#include <string>

class PHCompositeNode;

class HepMCNodeReader : public SubsysReco
{
 public:
  HepMCNodeReader(const std::string &name = "HEPMCREADER");
  virtual ~HepMCNodeReader() {}

  int Init(PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode);

  void Embed(const int i=1) {_embed_flag = i;}

private:
  int _embed_flag;
  
};

#endif

