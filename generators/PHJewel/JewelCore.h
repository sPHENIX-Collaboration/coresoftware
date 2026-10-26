// JewelCore.h -- drive a JEWEL 2.6.0 core library in-process.
//
// No Fun4All and no HepMC dependency, so the same code runs inside the
// PHJewel SubsysReco (RCF, HepMC2) and in the standalone jwgen_test (any
// machine with JEWEL's toolchain).
//
// The core library (libjewel_<medium>.so) is opened at run time with
// dlopen(RTLD_LOCAL) and exports only jw_* functions, so its private PYTHIA
// 6.4.25 and common blocks are invisible to the rest of the process.  JEWEL
// keeps its state in common blocks: one JewelCore per process.
//
// What PROGRAM JEWEL does, and what this reproduces:
//   init();  for each isospin channel: initpythia(b1,b2); genevent() x N
// For a fixed ISOCHANNEL (PP, PN, NP, NN) the sequence of calls -- and so the
// events -- are identical to the standalone binary with NJOB = seed/1000.
// For ISOCHANNEL XX the standalone program draws all channel counts up front
// from NEVENT, which a trigger-driven run does not know.  Here the channel is
// drawn per block of `block` JEWEL-accepted events with probabilities
// Z^2, Z(A-Z), (A-Z)Z, (A-Z)^2 over A^2, so the channel mix is exact in
// expectation; PYTHIA is re-initialised only when the channel changes.
#ifndef PHJEWEL_JEWELCORE_H
#define PHJEWEL_JEWELCORE_H

#include <dlfcn.h>

#include <cstdio>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

class JewelCore
{
 public:
  static constexpr int kAbiVersion = 1;

  struct EventInfo
  {
    double weight = 0;      // PARI(10)
    double xsec_mb = 0;     // PARI(1), running estimate of the current channel
    double b = -1;          // impact parameter [fm] (simple); LBRICK (brick); -1 (vac)
    double centrality = -1; // getcentrality(): percent (simple); TBRICK (brick); -1 (vac)
    int isub = 0;           // MSTI(1)
    int id1 = 0, id2 = 0;   // MSTI(15), MSTI(16)
    double x1 = 0, x2 = 0;  // PARI(33), PARI(34)
    double pthat = 0;       // PARI(17)
    double Q2 = 0;          // PARI(22)
    long ngood = 0;         // JEWEL's event number (HepMC E line, EVINFOFILE line)
    int nscatcen = 0;       // thermal partons stored for this event
    int channel = 0;        // 1 pp, 2 pn, 3 np, 4 nn
  };

  struct Counts
  {
    double ngood = 0, ndisc = 0, nstrange = 0, nevdone = 0, sumw = 0, wdisc = 0;
    double ntotspliti = 0, noverspliti = 0, ntotpdf = 0, noverpdf = 0;
    double ntotxsec = 0, noverxsec = 0, ntotsuda = 0, noversuda = 0;
  };

  struct Config  // what JEWEL itself read from the card
  {
    int mass = 0, nproton = 0, iso = -1;
    std::string process;
    bool weighted = false, keeprecoils = false, writescatcen = false, writedummies = false;
    bool dosubtraction = false, hadro = false, shorthepmc = false, isrscat = false, mpi = false;
    int nf = 0;
    double sqrts = 0, ptmin = 0, ptmax = 0, etamax = 0, q0 = 0, lqcd = 0, lps = 0;
  };

  JewelCore() = default;
  ~JewelCore()
  {
    // never dlclose: gfortran registers atexit handlers inside the library
  }
  JewelCore(const JewelCore &) = delete;
  JewelCore &operator=(const JewelCore &) = delete;

  // ---- setup ----------------------------------------------------------------
  // lib: "libjewel_simple.so" (found via LD_LIBRARY_PATH) or a full path
  void load(const std::string &lib)
  {
    // RTLD_DEEPBIND: the core resolves its own undefined symbols (pdfset_,
    // structm_ from LHAPDF 6) in its dependencies first. Without it they bind
    // to whatever is already global -- in a G4 job that is LHAPDF 5.9.1,
    // loaded by libPythia6 for the G4 decayer, and JEWEL stops in PDFSET.
    m_Handle = dlopen(lib.c_str(), RTLD_NOW | RTLD_LOCAL | RTLD_DEEPBIND);
    if (!m_Handle)
    {
      throw std::runtime_error("JewelCore: dlopen(" + lib + ") failed: " + dlerror());
    }
    sym(f_abi, "jw_abi_version");
    sym(f_err, "jw_last_error");
    sym(f_init, "jw_init");
    sym(f_chan, "jw_set_channel");
    sym(f_next, "jw_next");
    sym(f_text, "jw_event_text");
    sym(f_info, "jw_info");
    sym(f_counts, "jw_counts");
    sym(f_config, "jw_config");
    sym(f_finish, "jw_finish");
    if (f_abi() != kAbiVersion)
    {
      throw std::runtime_error("JewelCore: " + lib + " has ABI " + std::to_string(f_abi()) +
                               ", this code expects " + std::to_string(kAbiVersion) + "; rebuild the package");
    }
    m_Lib = lib;
  }

  // card   : canonical JEWEL card (JewelCard::writeCanonical)
  // hepmc  : "" or a path: JEWEL's own Fortran writer also writes the events
  // evinfo : "" or a path: EVINFOFILE
  void init(const std::string &card, const std::string &log, int seed, const std::string &hepmc = "",
            const std::string &evinfo = "")
  {
    if (!m_Handle) throw std::runtime_error("JewelCore::init before load");
    if (f_init(card.c_str(), log.c_str(), hepmc.c_str(), evinfo.c_str(), seed) != 0)
    {
      throw std::runtime_error(std::string("JewelCore: jw_init: ") + f_err());
    }
    m_Seed = seed;
    m_Config = readConfig();
    m_Rng.seed(static_cast<unsigned long long>(seed) * 2654435761ULL + 12345ULL);
    const double A = m_Config.mass, Z = m_Config.nproton, N = A - Z;
    m_Prob = {Z * Z / (A * A), Z * N / (A * A), N * Z / (A * A), N * N / (A * A)};
    m_Init = true;
  }

  void set_isospin_block(long n) { m_Block = n > 0 ? n : 1; }
  // give up when JEWEL discards this many passes in a row (0 = never)
  void set_max_discards(long n) { m_MaxDiscards = n; }

  // ---- events ---------------------------------------------------------------
  // Generate until JEWEL writes an event.  Returns its HepMC2 lines.
  const char *next()
  {
    if (!m_Init) throw std::runtime_error("JewelCore::next before init");
    if (m_Channel == 0 || (m_Config.iso == 0 && m_InBlock >= m_Block)) chooseChannel();
    long discards = 0;
    while (true)
    {
      ++m_Pass;
      if (f_next(static_cast<int>(m_Pass)) == 1) break;
      ++m_NDiscardedPasses;
      if (m_MaxDiscards > 0 && ++discards >= m_MaxDiscards)
      {
        throw std::runtime_error("JewelCore: " + std::to_string(discards) + " consecutive passes discarded by JEWEL");
      }
    }
    m_Last = readInfo();
    ++m_InBlock;
    ++m_SessN;
    ++m_NAccepted;
    m_SumW += m_Last.weight;
    m_ChanN[m_Channel - 1] += 1;
    m_ChanW[m_Channel - 1] += m_Last.weight;
    m_SessXsec = m_Last.xsec_mb;
    return f_text();
  }

  const EventInfo &info() const { return m_Last; }
  // full HepMC2 IO_GenEvent stream for the last event (header + event + footer)
  std::string event_stream() const
  {
    return std::string("HepMC::Version 2.06.05\nHepMC::IO_GenEvent-START_EVENT_LISTING\n") + f_text() +
           "HepMC::IO_GenEvent-END_EVENT_LISTING\n";
  }

  void finish()
  {
    if (m_Init && !m_Finished)
    {
      commitSession();
      f_finish();
      m_Finished = true;
    }
  }

  // ---- bookkeeping ----------------------------------------------------------
  long n_accepted() const { return m_NAccepted; }       // events JEWEL wrote
  long n_discarded_passes() const { return m_NDiscardedPasses; }
  double sum_weights() const { return m_SumW; }         // over events JEWEL wrote
  // cross section [mb] of the sample: per channel PYTHIA's estimate at the end
  // of each PYTHIA initialisation, averaged with the number of events made
  // under it.  For a fixed ISOCHANNEL this is PARI(1) at the last event.
  double xsec_mb() const
  {
    const double sn = m_XsecN + m_SessXsec * m_SessN;
    const double n = m_XsecCount + m_SessN;
    return n > 0 ? sn / n : 0.0;
  }
  long channel_events(int ich) const { return m_ChanN[ich - 1]; }
  double channel_sumw(int ich) const { return m_ChanW[ich - 1]; }
  int channel_switches() const { return m_Switches; }

  Counts counts() const
  {
    double d[16] = {0};
    f_counts(d);
    Counts c;
    c.ngood = d[0];
    c.ndisc = d[1];
    c.nstrange = d[2];
    c.nevdone = d[3];
    c.sumw = d[4];
    c.wdisc = d[5];
    c.ntotspliti = d[6];
    c.noverspliti = d[7];
    c.ntotpdf = d[8];
    c.noverpdf = d[9];
    c.ntotxsec = d[10];
    c.noverxsec = d[11];
    c.ntotsuda = d[12];
    c.noversuda = d[13];
    return c;
  }
  const Config &config() const { return m_Config; }
  const std::string &library() const { return m_Lib; }

 private:
  template <class F>
  void sym(F &f, const char *name)
  {
    f = reinterpret_cast<F>(dlsym(m_Handle, name));
    if (!f) throw std::runtime_error(std::string("JewelCore: ") + m_Lib + " lacks " + name);
  }

  Config readConfig() const
  {
    int ic[16] = {0};
    double dc[8] = {0};
    f_config(ic, dc);
    Config c;
    c.mass = ic[0];
    c.nproton = ic[1];
    c.iso = ic[2];
    for (int s = 24; s >= 0; s -= 8) c.process += static_cast<char>((ic[3] >> s) & 0xff);
    c.weighted = ic[4];
    c.keeprecoils = ic[5];
    c.writescatcen = ic[6];
    c.writedummies = ic[7];
    c.dosubtraction = ic[8];
    c.hadro = ic[9];
    c.shorthepmc = ic[10];
    c.isrscat = ic[11];
    c.mpi = ic[12];
    c.nf = ic[13];
    c.sqrts = dc[0];
    c.ptmin = dc[1];
    c.ptmax = dc[2];
    c.etamax = dc[3];
    c.q0 = dc[4];
    c.lqcd = dc[5];
    c.lps = dc[6];
    return c;
  }

  EventInfo readInfo() const
  {
    double d[16] = {0};
    f_info(d);
    EventInfo e;
    e.weight = d[0];
    e.xsec_mb = d[1];
    e.b = d[2];
    e.centrality = d[3];
    e.isub = static_cast<int>(d[4]);
    e.id1 = static_cast<int>(d[5]);
    e.id2 = static_cast<int>(d[6]);
    e.x1 = d[7];
    e.x2 = d[8];
    e.pthat = d[9];
    e.Q2 = d[10];
    e.ngood = static_cast<long>(d[11]);
    e.nscatcen = static_cast<int>(d[12]);
    e.channel = static_cast<int>(d[13]);
    return e;
  }

  // one PYTHIA initialisation = one cross-section estimate
  void commitSession()
  {
    m_XsecN += m_SessXsec * m_SessN;
    m_XsecCount += m_SessN;
    m_SessN = 0;
  }

  void chooseChannel()
  {
    int ich = m_Config.iso;  // 1..4 fixed
    if (ich == 0)
    {
      std::uniform_real_distribution<double> u(0.0, 1.0);
      const double r = u(m_Rng);
      double acc = 0;
      ich = 4;
      for (int i = 0; i < 4; ++i)
      {
        acc += m_Prob[i];
        if (r < acc)
        {
          ich = i + 1;
          break;
        }
      }
    }
    if (ich < 1 || ich > 4) throw std::runtime_error("JewelCore: JEWEL reports an unknown ISOCHANNEL");
    m_InBlock = 0;
    if (ich == m_Channel) return;  // same channel drawn again: PYTHIA keeps running
    commitSession();
    if (f_chan(ich) != 0) throw std::runtime_error(std::string("JewelCore: ") + f_err());
    m_Channel = ich;
    ++m_Switches;
  }

  void *m_Handle = nullptr;
  std::string m_Lib;
  int (*f_abi)() = nullptr;
  const char *(*f_err)() = nullptr;
  int (*f_init)(const char *, const char *, const char *, const char *, int) = nullptr;
  int (*f_chan)(int) = nullptr;
  int (*f_next)(int) = nullptr;
  const char *(*f_text)() = nullptr;
  void (*f_info)(double *) = nullptr;
  void (*f_counts)(double *) = nullptr;
  void (*f_config)(int *, double *) = nullptr;
  void (*f_finish)() = nullptr;

  bool m_Init = false;
  bool m_Finished = false;
  int m_Seed = 0;
  Config m_Config;
  std::mt19937_64 m_Rng;
  std::vector<double> m_Prob;
  long m_Block = 1000;
  long m_MaxDiscards = 100000;

  int m_Channel = 0;
  int m_Switches = 0;
  long m_InBlock = 0;  // events since the channel was last drawn
  long m_SessN = 0;    // events under the current PYTHIA initialisation
  double m_SessXsec = 0;
  long m_Pass = 0;
  long m_NAccepted = 0;
  long m_NDiscardedPasses = 0;
  double m_SumW = 0;
  double m_XsecN = 0;
  double m_XsecCount = 0;
  long m_ChanN[4] = {0, 0, 0, 0};
  double m_ChanW[4] = {0, 0, 0, 0};
  EventInfo m_Last;
};

#endif
