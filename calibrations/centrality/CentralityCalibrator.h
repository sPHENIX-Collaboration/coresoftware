// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef CENTRALITYCALIB_CENTRALITYCALIBRATOR_H
#define CENTRALITYCALIB_CENTRALITYCALIBRATOR_H

#include <array>
#include <string>

// Au+Au MBD centrality calibration for the three bunch classes (main, -5 ns and +5 ns bunches) that
// MinimumBiasClassifier and CentralityReco use, from the trees of MbdCentralityTree. Three steps:
//
//  1. FitT0Offset: the run's main-bunch MBD t0 (Gaussian core fit), so that t0 - offset puts the main bunch at 0.
//     The rows of every run, after T0TableHeader(), form the t0 offset table the other two steps read.
//  2. Calibrate: per class, the vertex scales (charge sum in each z bin balanced to the 0-5 cm bin) and the run scale
//     (the reference run's balanced mean charge sum over this run's). Calibrating the reference run itself (no
//     reference file) also gives the NBD+Glauber fit of its charge sum, whose trigger efficiency sets the number of
//     divisions, and the divisions: that many equal-population bins of its balanced charge sum.
//  3. WritePayloads: the run's Centrality, CentralityScale and CentralityVertexScale CDB files, with the classes as
//     column suffixes ("", "_m5", "_p5"). Divisions are always the reference run's; a run without its own calibration
//     gets the reference run's vertex scales and scale 1.
//
// Event selection: GL1 scaled trigger bit; >= 2 hit tubes per side (q > 0.4, |t_hit - t0 offset| < 25 ns), whose charge
// is summed; ZDC south and north > 40; |z_vtx| <= 60 cm. The class is decided by the event t0 - offset: main
// [-s_m/2, s_p/2], -5 ns [-3 s_m/2, -s_m/2), +5 ns (s_p/2, 3 s_p/2], with the bunch spacings s_m = s_p = 5 ns.
class CentralityCalibrator
{
 public:
  static constexpr int NCLASSES = 3;
  //! class directory names in the calibration files and payload column suffixes, in MinimumBiasInfo::BUNCH order
  static const std::array<std::string, NCLASSES> class_names;
  static const std::array<std::string, NCLASSES> class_suffixes;

  //! what WritePayloads used for a run
  enum PayloadSource
  {
    OWN = 0,
    REFERENCE_RUN = 1,      //!< the run is the reference run
    REFERENCE_DEFAULT = 2,  //!< no usable calibration of its own: reference vertex scales, scale 1
  };

  // event selection (defaults: Run 25 Au+Au; Run 24 Au+Au uses trigger bit 10)
  void SetTriggerBit(int bit) { m_trigger_bit = bit; }
  void SetVertexCut(double zmax) { m_vertex_cut = zmax; }
  void SetZdcCut(double emin) { m_zdc_cut = emin; }
  //! a tube is hit if q > qmin (float, as in MinimumBiasClassifier and CentralityReco) and |t_hit - t0 offset| < tmax
  void SetHitCuts(float qmin, double tmax, int nhits)
  {
    m_hit_charge_cut = qmin;
    m_hit_time_cut = tmax;
    m_nhit_cut = nhits;
  }
  void SetBunchSpacing(double spacing_m, double spacing_p)
  {
    m_spacing_m = spacing_m;
    m_spacing_p = spacing_p;
  }
  //! number of centrality divisions over the min-bias sample: the trigger efficiency of the NBD+Glauber fit in %
  void SetNDivisions(int n) { m_ndivisions = n; }
  //! Glauber Npart distribution (histogram "hNpart") for the NBD+Glauber fit of the reference run
  void SetGlauberFile(const std::string &file) { m_glauber_file = file; }
  //! balanced charge sums above this are left out of the divisions
  void SetMaxChargeSum(float qmax) { m_max_charge_sum = qmax; }

  // t0 offset fit: its selection has no t0 offset yet, a wider vertex cut and stops after n_max events
  void SetT0FitVertexCut(double zmax) { m_t0fit_vertex_cut = zmax; }
  void SetT0FitMaxEvents(long n) { m_t0fit_max_events = n; }
  //! the main-bunch fit is seeded at the tallest bin within +-window ns (default 3; wider when a whole run is shifted)
  void SetT0FitSeedWindow(double window) { m_t0fit_seed_window = window; }

  //! step 1: writes <outdir>/t0offset_<run>.csv (one table row) and <outdir>/t0offset_<run>.root (the t0 histogram)
  int FitT0Offset(int run, const std::string &treefile, const std::string &outdir) const;
  static std::string T0TableHeader() { return "run,nsel,offset,offset_err,sigma,chi2ndf,m5_rel,m5_rel_err,p5_rel,p5_rel_err,status"; }
  //! offset and status of one run in a t0 offset table; false if the run is not in it
  static bool ReadT0Offset(const std::string &t0table, int run, double &offset, int &status);

  //! step 2: reffile = the reference run's Calibrate output, or empty to calibrate the reference run itself (which
  //! needs the Glauber file). A run with an unusable t0 offset (status other than 0 or 3) is skipped: no output, and
  //! not an error unless it is the reference run.
  int Calibrate(int run, const std::string &treefile, const std::string &t0table, const std::string &reffile, const std::string &outfile) const;

  //! step 3: calibfile = the run's Calibrate output, or empty for none; returns a PayloadSource, or -1 on failure
  static int WritePayloads(int run, const std::string &calibfile, const std::string &reffile, const std::string &t0table, const std::string &outdir);

 private:
  int ClassOf(double t0) const;

  int m_trigger_bit{14};
  double m_vertex_cut{60.};
  double m_zdc_cut{40.};
  float m_hit_charge_cut{0.4F};
  double m_hit_time_cut{25.};
  int m_nhit_cut{2};
  double m_spacing_m{5.};
  double m_spacing_p{5.};
  int m_ndivisions{92};
  float m_max_charge_sum{2100};
  std::string m_glauber_file;

  double m_t0fit_vertex_cut{150.};
  long m_t0fit_max_events{500000};
  double m_t0fit_seed_window{3.};
};

#endif
