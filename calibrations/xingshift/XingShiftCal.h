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

  /** Called during initialization.
      Typically this is where you can book histograms, and e.g.
      register them to Fun4AllServer (so they can be output to file
      using Fun4AllServer::dumpHistos() method).
   */
  int Init(PHCompositeNode *topNode) override;


  /** Called for first event when run number is known.
      Typically this is where you may want to fetch data from
      database, because you know the run number. A place
      to book histograms which have to know the run number.
   */
  int InitRun(PHCompositeNode *topNode) override;


  /** Called for each event.
      This is where you do the real work.
   */
  int process_event(PHCompositeNode *topNode) override;


  /// Clean up internals after each event.
  int ResetEvent(PHCompositeNode *topNode) override;


  /// Called at the end of each run.
  int EndRun(const int runnumber) override;


  /// Called at the end of all processing.
  int End(PHCompositeNode *topNode) override;


  /// Reset
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

  const int packet_BEGINRUN = 900;
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
