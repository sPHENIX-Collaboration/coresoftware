// PHHerwig7Core.h
// This class aims at removing the need for a .run card, so .in card is the only input needed to run in Fun4All 
#ifndef PHHERWIG7CORE_H
#define PHHERWIG7CORE_H

#include <ThePEG/Config/Unitsystem.h>
#include <ThePEG/EventRecord/Event.h>
#include <ThePEG/Interface/InterfaceBase.h>
#include <ThePEG/Repository/BaseRepository.h>
#include <ThePEG/Repository/EventGenerator.h>
#include <ThePEG/Repository/Repository.h>

#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

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

  // Read the config (and setup) file, seed, initialize. Throws std::runtime_error on failure.
  void init()
  {
    if (m_ConfigFile.empty())
    {
      throw std::runtime_error("PHHerwig7Core::init - no config file set");
    }
    const char *offline_main = getenv("OFFLINE_MAIN");
    if (!offline_main)
    {
      throw std::runtime_error("PHHerwig7Core::init - OFFLINE_MAIN not set");
    }
    const std::string herwig_share = std::string(offline_main) + "/share/Herwig";

    // Equivalent to "Herwig read [card].in": load HerwigDefaults.rpo, set read path, read the card, then update()
    ThePEG::Repository::exitOnError() = 1; 
    std::string msg = ThePEG::Repository::load(herwig_share + "/HerwigDefaults.rpo");
    if (!msg.empty())
    {
      throw std::runtime_error("PHHerwig7Core::init - loading HerwigDefaults.rpo: " + msg);
    }
    // load() resets the read path, so the directories go in after it
    ThePEG::Repository::prependReadDir(m_ReadDirs);
    ThePEG::Repository::appendReadDir(herwig_share);
    for (const std::string &file : {m_ConfigFile, m_SetupFile})
    {
      if (file.empty())
      {
        continue;
      }
      // (motivated from Coderabbit review) safeguard against "saverun" or "run" in the .in card, which would create a .run file on disk
      std::ifstream card(file);
      std::string line;
      while (std::getline(card, line))
      {
        std::istringstream words(line);
        std::string command;
        words >> command;
        if (command == "saverun" || command == "run")
        {
          throw std::runtime_error("PHHerwig7Core::init - " + file + " contains " + command + "; PHHerwig7 creates the run in memory, remove it");
        }
      }
      msg = ThePEG::Repository::read(file, std::cout);
      if (!msg.empty())
      {
        throw std::runtime_error("PHHerwig7Core::init - reading " + file + ": " + msg);
      }
    }
    ThePEG::Repository::update();

    const std::string run_name = std::filesystem::path(m_ConfigFile).stem().string() + "-S" + std::to_string(m_Seed);
    // This replaces the Herwig "saverun", so the content of the .run card is only kept in memory instead of being produced as a separate file
    m_Generator = ThePEG::Repository::makeRun(ThePEG::Repository::GetObject<ThePEG::EGPtr>("/Herwig/Generators/EventGenerator"), run_name);
    m_Generator->setSeed(m_Seed);
    m_Generator->initialize();

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

  // event loop. Returns a null pointer only if Herwig gives up on the event
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

  // finalize
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

  std::string m_ConfigFile;
  std::string m_SetupFile;
  std::vector<std::string> m_ReadDirs;
  int m_Seed = 0;

 private:
  ThePEG::EGPtr m_Generator;
  bool m_Finished = true;
  long m_NGenerated = 0;
  double m_SumWeights = 0.0;
  long m_MaxEvents = -1;  // ThePEG: negative = unlimited
};

#endif
