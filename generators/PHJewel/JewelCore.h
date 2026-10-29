// JewelCore.h -- drive a JEWEL 2.6.0 core library in-process.
//
// No Fun4All and no HepMC dependency.  Only PHJewel uses it; the header is
// not installed.
//
// The core library (libjewel_<medium>.so) is opened at run time with
// dlopen(RTLD_LOCAL | RTLD_DEEPBIND) and exports only jw_* functions, so its private PYTHIA
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

#include <random>
#include <string>
#include <vector>

class JewelCore
{
 public:
  static constexpr int kAbiVersion = 1;

  struct EventInfo
  {
    double weight = 0;       // PARI(10)
    double xsec_mb = 0;      // PARI(1), running estimate of the current channel
    double b = -1;           // impact parameter [fm] (simple); LBRICK (brick); -1 (vac)
    double centrality = -1;  // getcentrality(): percent (simple); TBRICK (brick); -1 (vac)
    int isub = 0;            // MSTI(1)
    int id1 = 0, id2 = 0;    // MSTI(15), MSTI(16)
    double x1 = 0, x2 = 0;   // PARI(33), PARI(34)
    double pthat = 0;        // PARI(17)
    double Q2 = 0;           // PARI(22)
    long ngood = 0;          // JEWEL's event number (HepMC E line, EVINFOFILE line)
    int nscatcen = 0;        // thermal partons stored for this event
    int channel = 0;         // 1 pp, 2 pn, 3 np, 4 nn
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
  ~JewelCore() = default;  // never dlclose: gfortran registers atexit handlers inside the library
  JewelCore(const JewelCore &) = delete;
  JewelCore &operator=(const JewelCore &) = delete;
  JewelCore(JewelCore &&) = delete;
  JewelCore &operator=(JewelCore &&) = delete;

  // ---- setup ----------------------------------------------------------------
  // lib: "libjewel_simple.so" (found via LD_LIBRARY_PATH) or a full path
  void load(const std::string &lib);

  // card   : canonical JEWEL card (JewelCard::writeCanonical)
  // hepmc  : "" or a path: JEWEL's own Fortran writer also writes the events
  // evinfo : "" or a path: EVINFOFILE
  void init(const std::string &card, const std::string &log, int seed, const std::string &hepmc = "",
            const std::string &evinfo = "");

  void set_isospin_block(long n) { m_Block = n > 0 ? n : 1; }
  // give up when JEWEL discards this many passes in a row (0 = never)
  void set_max_discards(long n) { m_MaxDiscards = n; }

  // ---- events ---------------------------------------------------------------
  // Generate until JEWEL writes an event.  Returns its HepMC2 lines.
  const char *next();

  const EventInfo &info() const { return m_Last; }
  // full HepMC2 IO_GenEvent stream for the last event (header + event + footer)
  std::string event_stream() const;

  void finish();

  // ---- bookkeeping ----------------------------------------------------------
  long n_accepted() const { return m_NAccepted; }  // events JEWEL wrote
  long n_discarded_passes() const { return m_NDiscardedPasses; }
  double sum_weights() const { return m_SumW; }  // over events JEWEL wrote
  // cross section [mb] of the sample: per channel PYTHIA's estimate at the end
  // of each PYTHIA initialisation, averaged with the number of events made
  // under it.  For a fixed ISOCHANNEL this is PARI(1) at the last event.
  double xsec_mb() const;
  long channel_events(int ich) const { return m_ChanN[ich - 1]; }
  double channel_sumw(int ich) const { return m_ChanW[ich - 1]; }
  int channel_switches() const { return m_Switches; }

  Counts counts() const;
  const Config &config() const { return m_Config; }
  const std::string &library() const { return m_Lib; }

 private:
  template <class F>
  void sym(F &f, const char *name);

  Config readConfig() const;

  EventInfo readInfo() const;

  // one PYTHIA initialisation = one cross-section estimate
  void commitSession();

  void chooseChannel();

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
