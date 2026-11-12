#ifndef RHOUECALIB_RHOUECALIBPAYLOAD_H
#define RHOUECALIB_RHOUECALIBPAYLOAD_H

#include <string>

class TH1;

/// \class RhoUECalibPayload
///
/// \brief Turns the merged (hadd) output of the rho UE calibration makers into the knee
/// table and the CDBHistos payloads read by DetermineTowerRho and TowerRhoToBackground.
///
///   FitMbdKnees            MbdChargeHistMaker output -> knee table ("run knee knee_err sigma chi2ndf status events t0 t0_rms")
///   MakeKneePayloads       knee table                -> per-run knee payloads (MbdChargeKnee, set_mbd_knee_cdb_tag)
///   MakeTowerShapePayload  TowerShapeMapMaker output -> s-map payload (DetermineTowerRho::set_tower_shape_path)
///   MakeRhoUEPayload       RhoUECalibMaker output    -> UE payload (TowerRhoToBackground::set_eta_calib_path)
///
/// Payload contents (CDBHistos, format_version = 2):
///   zvtx_bins, mbdq_bins  TH1D, class bin edges in z [cm] and Q'; in the UE payload the
///                         content is the event-weighted bin centre (interpolation node)
///   mbdq_knee             TH1D over the run number, content = Q_knee
///   s-map: s_<L>          TH3F (ieta 24, iphi 64, class = izbin * n_Q' + iQ')
///   UE:    w_<L>, alpha_<L>, beta_<L>  TH3F (ieta, z, Q'); UE = w x or alpha x + beta
/// with <L> = CEMC (retowered), HCALIN, HCALOUT. All functions return 0 on success.
class RhoUECalibPayload
{
 public:
  // fit the central edge of each run's MBD charge spectrum ("q_<run>"):
  // f(Q) = (a + b (Q - Q_knee)) erfc((Q - Q_knee) / (sqrt2 sigma)) / 2, Poisson likelihood in
  // [0.55, 1.35] x the 99.5% quantile of Q; and the main-bunch t0 centre of each run ("t0_<run>":
  // the mode, then the mean within +- t0_window iterated five times) with its rms
  static int FitMbdKnees(const std::string &sums_file, const std::string &table, double t0_window = 2.5);

  // one CDBTTree knee payload per run of the knee table, <out_dir>/mbdq_knee_<run>.root, read by
  // DetermineTowerRho and TowerRhoToBackground (set_mbd_knee_cdb_tag / set_mbd_knee_path). Single
  // values: mbdq_knee (the knee in use), mbdq_knee_fit, mbdq_knee_sigma, events (floats), run,
  // status (ints). Runs whose fit failed, with fewer than min_events events, or whose edge width
  // differs from the median width by more than max_sigma_dev (relative) get the median knee of
  // the table and status 1 (0 = the run's own fit). <out_dir>/mbdq_knee_iov.txt lists the
  // payloads for the upload, each valid for its run.
  static int MakeKneePayloads(const std::string &table, const std::string &out_dir, double min_events = 1e5,
                              double max_sigma_dev = 0.5);

  // s = <E_T,tower> / <E_T> over the towers of the layer with >= min_entries entries, per class;
  // s = 1 for other towers, s <= 0.05, and whole (layer, class) maps whose two halves differ by
  // more than max_half_rms (rms over towers of (s_even - s_odd) / s: the map is noise there)
  static int MakeTowerShapePayload(const std::string &sums_file, const std::string &payload,
                                   double max_half_rms = 0.25, double min_entries = 100);

  // w = sy / sx and the weighted least-squares y = alpha x + beta per (layer, ieta, z, Q') bin
  // with >= min_events events; other bins take the nearest fitted Q' bin at the same z, then the
  // nearest z
  static int MakeRhoUEPayload(const std::string &sums_file, const std::string &payload, double min_events = 50);

 private:
  struct Knee
  {
    double knee{-1};
    double knee_err{-1};
    double sigma{-1};
    double chi2ndf{-1};
    int status{-1};
  };
  static Knee FitKnee(TH1 *hq);
  // knee per run from the makers' knee_<run> histograms, as a TH1D over the run number
  static TH1 *KneeHisto(const std::string &sums_file);
};

#endif
