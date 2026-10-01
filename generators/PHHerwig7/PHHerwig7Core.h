// PHHerwig7Core.h -- drive Herwig 7 in-process through its public API.
//
// No Fun4All and no HepMC dependency here, so the same code is exercised by
// the PHHerwig7 SubsysReco on RCF (HepMC2) and by the standalone smoke test
// on a machine with only Herwig (hwgen_test.cc, HepMC3).
//
// What Herwig's own `Herwig run file.run -N n -s seed` does, and what this
// reproduces, is (Herwig-7.3.0/API/HerwigAPI.cc, HerwigGenericRun):
//   PersistentIStream(runfile) >> eg;   eg->setSeed(seed);   eg->initialize();
//   loop: eg->shoot();                  eg->finalize();
// prepareRun() covers everything up to and including initialize().
#ifndef PHHERWIG7CORE_H
#define PHHERWIG7CORE_H

#include <Herwig/API/HerwigAPI.h>
#include <Herwig/API/HerwigUI.h>

#include <ThePEG/EventRecord/Event.h>
#include <ThePEG/Repository/EventGenerator.h>
#include <ThePEG/Repository/Repository.h>
#include <ThePEG/Repository/BaseRepository.h>
#include <ThePEG/Interface/InterfaceBase.h>
#include <ThePEG/Config/Unitsystem.h>

#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

// The command-line options Herwig's CLI would carry, as a plain struct.
class PHHerwigUI : public Herwig::HerwigUI
{
 public:
  PHHerwigUI() = default;
  ~PHHerwigUI() override = default;

  Herwig::RunMode::Mode runMode() const override { return Herwig::RunMode::RUN; }
  std::string repository() const override { return m_Repository; }
  std::string inputfile() const override { return m_RunFile; }
  std::string setupfile() const override { return m_SetupFile; }
  bool resume() const override { return false; }
  bool tics() const override { return false; }
  std::string tag() const override { return m_Tag; }
  std::string integrationList() const override { return std::string(); }
  const std::vector<std::string> &prependReadDirectories() const override { return m_Prepend; }
  const std::vector<std::string> &appendReadDirectories() const override { return m_Append; }
  long N() const override { return m_N; }
  int seed() const override { return m_Seed; }
  int jobs() const override { return 1; }
  unsigned int jobSize() const override { return 0; }
  unsigned int maxJobs() const override { return 0; }
  void quitWithHelp() const override { throw std::runtime_error("PHHerwigUI: Herwig asked to quit (see messages above)"); }
  void quit() const override { throw std::runtime_error("PHHerwigUI: Herwig asked to quit (see messages above)"); }
  std::ostream &outStream() const override { return std::cout; }
  std::ostream &errStream() const override { return std::cerr; }
  std::istream &inStream() const override { return std::cin; }

  std::string m_Repository = "HerwigDefaults.rpo";
  std::string m_RunFile;
  std::string m_SetupFile;
  std::string m_Tag;
  std::vector<std::string> m_Prepend;
  std::vector<std::string> m_Append;
  long m_N = 0;
  int m_Seed = 0;
};

class PHHerwig7Core
{
 public:
  PHHerwig7Core() = default;
  ~PHHerwig7Core()
  {
    if (m_Generator && !m_Finished)
    {
      finish();
    }
  }
  PHHerwig7Core(const PHHerwig7Core &) = delete;
  PHHerwig7Core &operator=(const PHHerwig7Core &) = delete;

  PHHerwigUI &ui() { return m_UI; }

  // Load the .run file, seed, initialize.  Throws std::runtime_error on failure.
  void init()
  {
    if (m_UI.m_RunFile.empty())
    {
      throw std::runtime_error("PHHerwig7Core::init - no run file set");
    }
    ThePEG::Repository::exitOnError() = 1;  // as the Herwig CLI does
    m_Generator = Herwig::API::prepareRun(m_UI);
    if (!m_Generator)
    {
      throw std::runtime_error("PHHerwig7Core::init - prepareRun returned no EventGenerator");
    }
    // The .run file carries the card's NumberOfEvents and shoot() returns
    // nothing once that many have been generated.  Fun4All owns the event
    // count here, so lift the limit (Herwig's own -N does the same).
    // N(long) is protected; go through the interface as a setup file would.
    const ThePEG::InterfaceBase *ifc = ThePEG::BaseRepository::FindInterface(m_Generator, "NumberOfEvents");
    if (ifc)
    {
      ifc->exec(*m_Generator, "set", std::to_string(m_MaxEvents));
    }
    else
    {
      std::cerr << "PHHerwig7Core::init - cannot find NumberOfEvents interface; generator will stop after "
                << m_Generator->N() << " events" << std::endl;
    }
    m_Finished = false;
  }

  // generator-side cap on events shot; default -1 = unlimited (Fun4All decides)
  void set_max_events(long n) { m_MaxEvents = n; }

  // One event.  Returns a null pointer only if Herwig gives up on the event.
  ThePEG::EventPtr next()
  {
    ThePEG::EventPtr ev = m_Generator->shoot();
    if (ev)
    {
      ++m_NGenerated;
      m_SumWeights += ev->weight();
    }
    return ev;
  }

  void finish()
  {
    if (m_Generator && !m_Finished)
    {
      m_Generator->finalize();
      m_Finished = true;
    }
  }

  ThePEG::EventGenerator *generator() { return m_Generator ? &(*m_Generator) : nullptr; }

  // integrated cross section of the running sampler, in pb
  double xsec_pb() const
  {
    if (!m_Generator) return 0.0;
    const double x = m_Generator->integratedXSec() / ThePEG::Units::picobarn;
    return x;
  }
  double xsec_err_pb() const
  {
    if (!m_Generator) return 0.0;
    const double x = m_Generator->integratedXSecErr() / ThePEG::Units::picobarn;
    return x;
  }
  long n_generated() const { return m_NGenerated; }
  double sum_weights() const { return m_SumWeights; }

 private:
  PHHerwigUI m_UI;
  ThePEG::EGPtr m_Generator;
  bool m_Finished = true;
  long m_NGenerated = 0;
  double m_SumWeights = 0.0;
  long m_MaxEvents = -1;  // ThePEG: negative = unlimited
};

#endif
