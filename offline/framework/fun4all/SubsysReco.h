// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef FUN4ALL_SUBSYSRECO_H
#define FUN4ALL_SUBSYSRECO_H

#include "Fun4AllBase.h"

#include <string>

class PHCompositeNode;

/** Base class for all reconstruction and analysis modules to be
 *  used under the Fun4All framework.
 *
 *  If you write a reconstruction/analysis module, you must derive
 *  from this base class and you have to implement this class methods.
 *  None of these are strictly required as far as C++ is concerned, but as
 *  far as your job is concerned, at least process_event(), to do the
 *  job, and InitRun(), to initialize, should be implemented.  
 *  
 */

class SubsysReco : public Fun4AllBase
{
 public:
  /** dtor. 
      Does nothing as this is a base class only.
  */
  ~SubsysReco() override {}

  /// Called at the end of all processing.
  virtual int End(PHCompositeNode * /*topNode*/) { return 0; }

  /// Called at the end of each run.
  virtual int EndRun(const int /*runnumber*/) { return 0; }

  /** Called during initialization.
      Typically this is where you can book histograms, and e.g.
      register them to Fun4AllServer (so they can be output to file
      using Fun4AllServer::dumpHistos() method).
   */
  virtual int Init(PHCompositeNode * /*topNode*/) { return 0; }

  /** Called for first event when run number is known.
      Typically this is where you may want to fetch data from
      database, because you know the run number.
   */
  virtual int InitRun(PHCompositeNode * /*topNode*/) { return 0; }

  /** Called for each event.
      This is where you do the real work.
  */
  virtual int process_event(PHCompositeNode * /*topNode*/) { return 0; }

  /// Reset.
  virtual int Reset(PHCompositeNode * /*topNode*/) { return 0; }

  /// Clean up after each event.
  virtual int ResetEvent(PHCompositeNode * /*topNode*/) { return 0; }

  void Print(const std::string & /*what*/ = "ALL") const override {}

 protected:
  /** ctor.
      @param name is the reference used inside the Fun4AllServer
  */
  SubsysReco(const std::string &name = "NONAME")
    : Fun4AllBase(name)
  {
  }
};

#endif
