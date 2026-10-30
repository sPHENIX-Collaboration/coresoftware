// jewel_capi.cc -- the C interface of a JEWEL core library (libjewel_<medium>).
//
// Built into the same shared library as the JEWEL Fortran and its private
// copy of PYTHIA 6.4.25.  Only the jw_* symbols below are exported (libtool
// -export-symbols-regex '^jw_' / a macOS exported-symbols list), so the
// library's PYTHIA, common blocks and Fortran routines can neither clash with
// nor be interposed by any other PYTHIA 6 in the process.
//
// No HepMC, no Fun4All: an event leaves here as HepMC2 IO_GenEvent text, the
// same lines JEWEL writes to its HEPMCFILE but with full double precision.
// The caller parses it with HepMC's own reader.
//
// JEWEL keeps all its state in Fortran common blocks, so a process can hold
// exactly one initialised JEWEL per loaded core library.  jw_init refuses a
// second call.

#include <cstdio>
#include <cstring>
#include <string>

extern "C"
{
  // ---- Fortran side (jewel-api.f) -----------------------------------------
  void jwfinit_();
  void jwfchan_(int *ich);
  void jwfnext_(int *iev);
  void jwfinfo_(double *d);
  void jwfcounts_(double *d);
  void jwfconfig_(int *ic, double *dc);
  void jwffinish_();

  // COMMON /JWAPI/ APIMODE, APIWRITE, APISEED   (LOGICAL = 4-byte int)
  extern struct
  {
    int apimode;
    int apiwrite;
    int apiseed;
  } jwapi_;
  // COMMON /JWAPIC/ APICARD, APILOG, APIHEPMC, APIEVINFO   (CHARACTER*512)
  extern struct
  {
    char card[512];
    char log[512];
    char hepmc[512];
    char evinfo[512];
  } jwapic_;
}

namespace
{
  std::string g_event;    // lines of the event being written / last written
  bool g_open = false;    // E line seen, HMCEND not yet
  bool g_ready = false;   // a complete event is in g_event
  bool g_init = false;
  bool g_finished = false;
  std::string g_error;

  void put_fortran_string(char *dst, const char *src)
  {
    const size_t n = src ? std::strlen(src) : 0;
    std::memset(dst, ' ', 512);
    std::memcpy(dst, src, n < 512 ? n : 512);
  }

  void add_int(int v)
  {
    char b[32];
    std::snprintf(b, sizeof(b), " %d", v);
    g_event += b;
  }
  void add_dbl(double v)
  {
    char b[40];
    std::snprintf(b, sizeof(b), " %.17g", v);
    g_event += b;
  }
}  // namespace

// ---- event buffer, called by the HMC* helpers in jewel-lib.f ---------------
extern "C"
{
  void jwbufe_(int *i1, int *i2, double *e1, double *e2, double *e3, int *i3, int *i4, int *i5,
               int *i6, int *i7, int *i8, int *i9, double *e4)
  {
    g_event.clear();
    g_open = true;
    g_ready = false;
    g_event += "E";
    add_int(*i1);
    add_int(*i2);
    add_dbl(*e1);
    add_dbl(*e2);
    add_dbl(*e3);
    add_int(*i3);
    add_int(*i4);
    add_int(*i5);
    add_int(*i6);
    add_int(*i7);
    add_int(*i8);
    add_int(*i9);
    add_dbl(*e4);
    g_event += "\n";
  }
  void jwbufn_() { g_event += "N 1 \"0\"\n"; }
  void jwbufu_() { g_event += "U GEV MM\n"; }
  void jwbufc_(double *e1, double *e2)
  {
    g_event += "C";
    add_dbl(*e1);
    add_dbl(*e2);
    g_event += "\n";
  }
  void jwbufh_(int *i1, int *i2, int *i3, int *i4, int *i5, int *i6, int *i7, int *i8, int *i9,
               double *e1, double *e2, double *e3, double *e4)
  {
    g_event += "H";
    add_int(*i1);
    add_int(*i2);
    add_int(*i3);
    add_int(*i4);
    add_int(*i5);
    add_int(*i6);
    add_int(*i7);
    add_int(*i8);
    add_int(*i9);
    add_dbl(*e1);
    add_dbl(*e2);
    add_dbl(*e3);
    add_dbl(*e4);
    g_event += "\n";
  }
  void jwbuff_(int *i1, int *i2, double *e1, double *e2, double *e3, double *e4, double *e5,
               int *i3, int *i4)
  {
    g_event += "F";
    add_int(*i1);
    add_int(*i2);
    add_dbl(*e1);
    add_dbl(*e2);
    add_dbl(*e3);
    add_dbl(*e4);
    add_dbl(*e5);
    add_int(*i3);
    add_int(*i4);
    g_event += "\n";
  }
  void jwbufv_(int *i1, int *i2, int *i3, int *i4, int *i5, int *i6, int *i7, int *i8, int *i9)
  {
    g_event += "V";
    add_int(*i1);
    add_int(*i2);
    add_int(*i3);
    add_int(*i4);
    add_int(*i5);
    add_int(*i6);
    add_int(*i7);
    add_int(*i8);
    add_int(*i9);
    g_event += "\n";
  }
  void jwbufp_(int *i1, int *i2, double *e1, double *e2, double *e3, double *e4, double *e5,
               int *i3, int *i4, int *i5, int *i6, int *i7)
  {
    g_event += "P";
    add_int(*i1);
    add_int(*i2);
    add_dbl(*e1);
    add_dbl(*e2);
    add_dbl(*e3);
    add_dbl(*e4);
    add_dbl(*e5);
    add_int(*i3);
    add_int(*i4);
    add_int(*i5);
    add_int(*i6);
    add_int(*i7);
    g_event += "\n";
  }
  void jwbufend_()
  {
    if (g_open)
    {
      g_open = false;
      g_ready = true;
    }
  }
}

// ---- exported interface ----------------------------------------------------
extern "C"
{
  // bump when any signature or array layout below changes
  int jw_abi_version() { return 1; }

  const char *jw_last_error() { return g_error.c_str(); }

  // card    : the (canonical) JEWEL parameter file
  // logfile : JEWEL's LOGFILE (PYTHIA writes there too)
  // hepmc   : "" = no file; else JEWEL's own Fortran HepMC writer also writes
  //           here (validation: must be identical to a standalone run)
  // evinfo  : "" = none; else EVINFOFILE, one line per accepted event
  // seed    : PYTHIA MRPY(1), 0 <= seed < 900000000
  int jw_init(const char *card, const char *logfile, const char *hepmc, const char *evinfo, int seed)
  {
    if (g_init)
    {
      g_error = "jw_init called twice: JEWEL keeps its state in common blocks and cannot be re-initialised";
      return 1;
    }
    if (seed < 0 || seed >= 900000000)
    {
      g_error = "seed outside the PYTHIA 6 range 0 <= seed < 900000000";
      return 2;
    }
    const char *args[4] = {card, logfile, hepmc, evinfo};
    for (const char *a : args)
    {
      if (a && std::strlen(a) >= 512)
      {
        g_error = std::string("path longer than 511 characters: ") + a;
        return 3;
      }
    }
    put_fortran_string(jwapic_.card, card);
    put_fortran_string(jwapic_.log, logfile);
    put_fortran_string(jwapic_.hepmc, hepmc);
    put_fortran_string(jwapic_.evinfo, evinfo);
    jwapi_.apimode = 1;
    jwapi_.apiwrite = (hepmc && hepmc[0]) ? 1 : 0;
    jwapi_.apiseed = seed;
    jwfinit_();
    g_init = true;
    return 0;
  }

  // 1 pp, 2 pn, 3 np, 4 nn: re-initialises PYTHIA for that beam pair
  int jw_set_channel(int ich)
  {
    if (!g_init || g_finished || ich < 1 || ich > 4)
    {
      g_error = "jw_set_channel: not initialised, finished, or channel outside 1..4";
      return 1;
    }
    jwfchan_(&ich);
    return 0;
  }

  // One pass of JEWEL's event loop.  Returns 1 if an event was written,
  // 0 if JEWEL discarded this one (its own DISCARD / PYTHIA error paths).
  int jw_next(int iev)
  {
    g_ready = false;
    g_open = false;
    jwfnext_(&iev);
    return g_ready ? 1 : 0;
  }

  // the last event's HepMC2 lines (E ... P), without file header/footer
  const char *jw_event_text() { return g_event.c_str(); }

  // d[16]: weight, xsec_mb, b, centrality, isub, id1, id2, x1, x2, pthat, Q2,
  //        NGOOD, nscatcen, channel, 0, 0
  void jw_info(double *d) { jwfinfo_(d); }

  // d[16]: NGOOD NDISC NSTRANGE nevdone sumw wdisc ntotspliti noverspliti
  //        ntotpdf noverpdf ntotxsec noverxsec ntotsuda noversuda 0 0
  void jw_counts(double *d) { jwfcounts_(d); }

  // ic[16]: mass nproton iso(0 XX 1 PP 2 PN 3 NP 4 NN) process(packed)
  //         weighted keeprecoils writescatcen writedummies dosubtraction
  //         hadro shorthepmc isrscat mpi nf 0 0
  // dc[8] : sqrts ptmin ptmax etamax q0 lambdaqcd lambdaqcdps 0
  void jw_config(int *ic, double *dc) { jwfconfig_(ic, dc); }

  void jw_finish()
  {
    if (g_init && !g_finished)
    {
      jwffinish_();
      g_finished = true;
    }
  }
}
