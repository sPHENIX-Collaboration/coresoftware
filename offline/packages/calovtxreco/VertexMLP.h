#ifndef VERTEXMLP_H
#define VERTEXMLP_H

// Loads and evaluates the calo-vertex-z MLP exported by
// export_nn_weights.py + export_nn_root.C (calovertex/vertex_mlp_weights.root).
//
// Usage in a Fun4All SubsysReco module:
//   in InitRun():        VertexMLP::Load("vertex_mlp_weights.root");
//   in process_event():  double z = VertexMLP::PredictVertexZ(features);
//
// Retraining the model requires no recompilation here: rerun
// export_nn_weights.py + export_nn_root.C to regenerate the .root file,
// point Load() at it (same path or a new one), done.
//
// Input feature order (features[0..20] must be filled in exactly this
// order -- it's the order the model was trained on). Features 0-3 are
// binary "this calo had a valid shower" flags (1 = at least one constituent
// above threshold, 0 = none -- fill the corresponding zmean/zsig/zskew with
// 0 in that case, matching how the training CSV imputes them):
//   0  emcal_lead (flag)      7  emcal_lead_energy     14 emcal_sublead_zskew
//   1  ohcal_lead (flag)      8  ohcal_lead_zmean      15 emcal_sublead_energy
//   2  emcal_sublead (flag)   9  ohcal_lead_zsig       16 ohcal_sublead_zmean
//   3  ohcal_sublead (flag)   10 ohcal_lead_zskew      17 ohcal_sublead_zsig
//   4  emcal_lead_zmean       11 ohcal_lead_energy     18 ohcal_sublead_zskew
//   5  emcal_lead_zsig        12 emcal_sublead_zmean   19 ohcal_sublead_energy
//   6  emcal_lead_zskew       13 emcal_sublead_zsig    20 exj
// (also stored, in this order, in the "feature_names" TObjArray inside the
// .root file -- worth checking against if this list and that file ever
// disagree).

#include <TFile.h>
#include <TMatrixD.h>
#include <TVectorD.h>

#include <algorithm>
#include <array>
#include <iostream>
#include <string>

namespace VertexMLP {

constexpr int kNFeatures = 21;

inline TMatrixD gW1, gW2, gW3, gW4;
inline TVectorD gB1, gB2, gB3, gB4;
inline TVectorD gMedian, gIqr, gLo, gHi;
inline bool gLoaded = false;

namespace detail {

inline bool GetMatrix(TFile *f, const char *name, TMatrixD &out) {
  TMatrixD *m = static_cast<TMatrixD *>(f->Get(name));
  if (!m) {
    std::cerr << "VertexMLP::Load: missing TMatrixD \"" << name << "\"" << std::endl;
    return false;
  }
  out.ResizeTo(*m);
  out = *m;
  return true;
}

inline bool GetVector(TFile *f, const char *name, TVectorD &out) {
  TVectorD *v = static_cast<TVectorD *>(f->Get(name));
  if (!v) {
    std::cerr << "VertexMLP::Load: missing TVectorD \"" << name << "\"" << std::endl;
    return false;
  }
  out.ResizeTo(*v);
  out = *v;
  return true;
}

// out = ReLU(W*in + b); set relu=false for the final (output) layer.
inline TVectorD ApplyLayer(const TMatrixD &W, const TVectorD &b, const TVectorD &in, bool relu) {
  TVectorD out = W * in + b;
  if (relu) {
    for (int i = 0; i < out.GetNrows(); i++) out[i] = std::max(0.0, out[i]);
  }
  return out;
}

} // namespace detail

// Call once, from InitRun(). Returns false (and leaves the model unloaded,
// so PredictVertexZ() will refuse to run) if the file or any expected
// object inside it is missing.
inline bool Load(const std::string &filename = "vertex_mlp_weights.root") {
  gLoaded = false;

  TFile *f = TFile::Open(filename.c_str(), "READ");
  if (!f || f->IsZombie()) {
    std::cerr << "VertexMLP::Load: cannot open " << filename << std::endl;
    return false;
  }

  bool ok = detail::GetMatrix(f, "W1", gW1) && detail::GetMatrix(f, "W2", gW2) &&
            detail::GetMatrix(f, "W3", gW3) && detail::GetMatrix(f, "W4", gW4) &&
            detail::GetVector(f, "b1", gB1) && detail::GetVector(f, "b2", gB2) &&
            detail::GetVector(f, "b3", gB3) && detail::GetVector(f, "b4", gB4) &&
            detail::GetVector(f, "median", gMedian) && detail::GetVector(f, "iqr", gIqr) &&
            detail::GetVector(f, "lo", gLo) && detail::GetVector(f, "hi", gHi);

  f->Close();
  delete f;

  if (!ok || gMedian.GetNrows() != kNFeatures) {
    std::cerr << "VertexMLP::Load: failed to load a complete model from " << filename << std::endl;
    return false;
  }

  gLoaded = true;
  std::cout << "VertexMLP::Load: loaded " << filename << std::endl;
  return true;
}

// Predicts the vertex z [cm] given the 21 input features (see the ordering
// table above). Returns 0 and prints an error if Load() hasn't succeeded.
inline double PredictVertexZ(const std::array<double, kNFeatures> &features) {
  if (!gLoaded) {
    std::cerr << "VertexMLP::PredictVertexZ: model not loaded -- call "
                 "VertexMLP::Load() in InitRun() first"
              << std::endl;
    return 0.0;
  }

  // Same preprocessing as training: percentile clip, then median/IQR scale.
  TVectorD x(kNFeatures);
  for (int i = 0; i < kNFeatures; i++) {
    double v = std::clamp(features[i], gLo[i], gHi[i]);
    x[i] = (v - gMedian[i]) / gIqr[i];
  }

  TVectorD h1 = detail::ApplyLayer(gW1, gB1, x, true);
  TVectorD h2 = detail::ApplyLayer(gW2, gB2, h1, true);
  TVectorD h3 = detail::ApplyLayer(gW3, gB3, h2, true);
  TVectorD out = detail::ApplyLayer(gW4, gB4, h3, false);

  return out[0];
}

} // namespace VertexMLP

#endif
