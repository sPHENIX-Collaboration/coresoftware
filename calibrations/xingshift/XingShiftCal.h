// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef XINGSHIFTCAL_H
#define XINGSHIFTCAL_H

#include <fun4all/SubsysReco.h>

#include <string>

class PHCompositeNode;
class Event;
class Packet;

class XingShiftCal : public SubsysReco
{
 public:

  XingShiftCal(const std::string &name = "XingShiftCal");

  ~XingShiftCal() override;


  int Init(PHCompositeNode *topNode) override;
  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;
  int ResetEvent(PHCompositeNode *topNode) override;
  int EndRun(const int runnumber) override;
  int End(PHCompositeNode *topNode) override;
  int Reset(PHCompositeNode * /*topNode*/) override;

  void Print(const std::string &what = "ALL") const override;

 private:
  Packet *p;
  int nevt = 0;
  int threshold = 5000;
  int evtcap = 50000;
  int done = 0;

  static const int NTRIG = 16;
  static const int NBUNCHES = 120;

  const int packet_BEGRUN = 900;
  const int packet_GL1 = 14001;

  bool success = 0;

  //default xingshift
  int xingshift = 5;

  long long scalercounts[NTRIG][NBUNCHES];

 public:
  int Calibrate(const int final=0);
  int CalculateCrossingShift(int& xingshift, long long counts[NTRIG][NBUNCHES], bool& success);


};

#endif // XINGSHIFTCAL_H
