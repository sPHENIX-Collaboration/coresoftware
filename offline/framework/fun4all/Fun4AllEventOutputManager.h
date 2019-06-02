// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef FUN4ALL_FUN4ALLEVENTOUTPUTMANAGER_H
#define FUN4ALL_FUN4ALLEVENTOUTPUTMANAGER_H

#include "Fun4AllOutputManager.h"

#include <string>

class Fun4AllEventOutStream;
class PHCompositeNode;

class Fun4AllEventOutputManager : public Fun4AllOutputManager
{
 public:
  Fun4AllEventOutputManager(const std::string &myname = "EVENTOUT", const std::string &filename = "eventout.prdf", const unsigned int sizeInMB = 0, const int offset = 0, const int increment = 1);
  virtual ~Fun4AllEventOutputManager();

  int outfileopen(const std::string & /*fname*/) { return 0; }

  void Print(const std::string &what = "ALL") const;

  int Write(PHCompositeNode *startNode);

  int AddPacket(const int ipkt);
  int DropPacket(const int ipkt);
  int AddPacketRange(const int ipktmin, const int ipktmax);
  int DropPacketRange(const int ipktmin, const int ipktmax);
  void SetOutfileName(const std::string &fname);

 protected:
  std::string outfilerule;
  Fun4AllEventOutStream *outstream;
};

#endif /* FUN4ALL_FUN4ALLEVENTOUTPUTMANAGER_H */
