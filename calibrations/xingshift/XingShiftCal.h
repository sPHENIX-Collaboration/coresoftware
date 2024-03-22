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
  static const int NTRIG = 16;
  static const int NBUNCHES = 120;

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

  int Calibrate(const int final = 0);
  int CalculateCrossingShift(int &xingshift, uint64_t counts[NTRIG][NBUNCHES], bool &success);

 private:
  Packet *p = nullptr;
  Packet *pBlueSpin = nullptr;
  Packet *pYellSpin = nullptr;
  int nevt = 0;
  int threshold = 1000;
  int evtcap = 50000;
  int done = 0;

  const int packet_BLUESPIN = 14902;
  const int packet_YELLSPIN = 14903;
  const int packet_GL1 = 14001;

  bool success = 0;

  // default xingshift
  int xingshift = 5;

  int blueSpinPattern[NBUNCHES] = {0};
  int yellSpinPattern[NBUNCHES] = {0};

  uint64_t scalercounts[NTRIG][NBUNCHES]{};
};

#endif  // XINGSHIFTCAL_H
