// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef XINGSHIFT_XINGSHIFTCAL_H
#define XINGSHIFT_XINGSHIFTCAL_H

#include <fun4all/SubsysReco.h>

#include <cstdint>
#include <string>

class PHCompositeNode;
class Packet;

class XingShiftCal : public SubsysReco
{
  static const int NTRIG = 16;
  static const int NBUNCHES = 120;

 public:
  XingShiftCal(const std::string &name = "XingShiftCal", const int poverwriteSpinEntry = 0);

  ~XingShiftCal() override;

  int Init(PHCompositeNode *topNode) override;
  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;
  int ResetEvent(PHCompositeNode *topNode) override;
  // int EndRun(const int runnumber) override;
  int End(PHCompositeNode *topNode) override;
  int Reset(PHCompositeNode * /*topNode*/) override;

  void Print(const std::string &what = "ALL") const override;

  int Calibrate(const int final = 0);
  int CalculateCrossingShift(int &xingshift, uint64_t counts[NTRIG][NBUNCHES], bool &success);
  int WriteToCDB(const std::string &fname);
  int CommitToSpinDB();
  std::string SQLArrayConstF(float x, int n);

 private:
  Packet *p{nullptr};

  //  Packet *pBlueSpin {nullptr};
  //  Packet *pYellSpin {nullptr};
  Packet *pBluePol{nullptr};
  Packet *pYellPol{nullptr};
  //  Packet *pBlueAsym {nullptr};
  //  Packet *pYellAsym {nullptr};
  Packet *pBlueIntPattern{nullptr};
  Packet *pYellIntPattern{nullptr};
  Packet *pBluePolPattern{nullptr};
  Packet *pYellPolPattern{nullptr};
  Packet *pBlueFillNumber{nullptr};
  Packet *pYellFillNumber{nullptr};

  int nevt{0};
  int threshold{1000};
  int evtcap{50000};
  int done{0};

  //  const int packet_BLUESPIN  {14902};
  //  const int packet_YELLSPIN  {14903};
  const int packet_BLUEPOL{14905};
  //  const int packet_BLUEASYM  {14906};
  const int packet_YELLPOL{14907};
  //  const int packet_YELLASYM  {14908};
  const int packet_BLUEINTPATTERN{14910};
  const int packet_BLUEPOLPATTERN{14911};
  const int packet_YELLINTPATTERN{14912};
  const int packet_YELLPOLPATTERN{14913};
  const int packet_BLUEFILLNUMBER{14915};
  const int packet_YELLFILLNUMBER{14916};
  const int packet_GL1{14001};

  int runnumber{0};
  bool success{0};
  int commitSuccessCDB{0};
  int commitSuccessSpinDB{0};
  int overwriteSpinEntry{0};

  // default xingshift
  int xingshift{0};

  int blueSpinPattern[NBUNCHES]{0};
  int yellSpinPattern[NBUNCHES]{0};
  int blueFillPattern[NBUNCHES]{0};
  int yellFillPattern[NBUNCHES]{0};

  float polBlue{0};
  float polBlueErr{0};
  float polYellow{0};
  float polYellowErr{0};

  int fillnumberBlue{0};
  int fillnumberYellow{0};

  uint64_t scalercounts[NTRIG][NBUNCHES]{};
};

#endif  // XINGSHIFT_XINGSHIFTCAL_H
