#include "LaserClusterTruthMatcher.h"

// TODO: fix these include paths / class names to match the real headers --
// written from memory of the macro and the conventions described, not from
// the actual LaserClusterContainer / LaserClusterv5 / TpcDistortionCorrection
// headers, since those weren't available while drafting this.
#include <trackbase/LaserClusterContainer.h>
#include <trackbase/LaserClusterv5.h>
#include <trackbase/TpcDefs.h>
#include <trackbase/TrkrDefs.h>

#include <cdbobjects/CDBTTree.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>


#include <TCanvas.h>
#include <TColor.h>
#include <TGraph.h>
#include <TH1D.h>
#include <TLine.h>
#include <TMath.h>
#include <TStyle.h>

#include <algorithm>
#include <cmath>
#include <format>
#include <iostream>
#include <numeric>
#include <set>

// ============================================================================
//  row/stripe finding + reco-to-truth row matching -- unchanged from the
//  standalone macro (assignClustersToTruth.C), just kept file-local here.
//  See that macro for the extensive comments on WHY each piece works the way
//  it does; only the driving code around it (module vs. macro main) changed.
// ============================================================================

namespace
{
  static constexpr double kPetalDefault = TMath::Pi() / 9.0;  // 20 degrees

  struct Cluster
  {
    double R;
    double phi;
    double adc;
    TrkrDefs::cluskey key;  // index into the module's per-side working vector, NOT a tree entry
  };

  struct Result
  {
    std::vector<double> grid;
    std::vector<double> density;
    std::vector<double> peakR;
    std::vector<double> bound;
  };

  double layerPitch(double R)
  {
    if (R < 40.0) return 0.625;
    if (R < 60.0) return 1.25;
    return 1.125;
  }

  double clusSigma(const Cluster &c, double sigma_ref, double adc_ref)
  {
    double s = sigma_ref * sqrt(adc_ref / std::max(c.adc, 1.0));
    s = std::max(s, 0.5 * layerPitch(c.R));
    s = std::min(s, 1.0);
    return s;
  }

  Result findRows(const std::vector<Cluster> &cl,
                          double sigma_ref = 0.15, double promfrac = 0.05,
                          double rlo = 28.0, double rhi = 78.0, int ngrid = 4000)
  {
    Result out;
    if (cl.empty()) return out;

    std::vector<double> adcs;
    adcs.reserve(cl.size());
    for (const auto &c : cl) adcs.push_back(c.adc);
    std::nth_element(adcs.begin(), adcs.begin() + adcs.size() / 2, adcs.end());
    const double adc_ref = std::max(adcs[adcs.size() / 2], 1.0);

    const double dr = (rhi - rlo) / ngrid;
    out.grid.resize(ngrid);
    out.density.assign(ngrid, 0.0);
    for (int b = 0; b < ngrid; ++b) out.grid[b] = rlo + (b + 0.5) * dr;

    for (const auto &c : cl)
    {
      if (c.R < rlo || c.R > rhi) continue;
      const double s = clusSigma(c, sigma_ref, adc_ref);
      const int lo = std::max(0, (int)((c.R - 4 * s - rlo) / dr));
      const int hi = std::min(ngrid - 1, (int)((c.R + 4 * s - rlo) / dr));
      for (int b = lo; b <= hi; ++b)
      {
        const double d = out.grid[b] - c.R;
        out.density[b] += exp(-0.5 * d * d / (s * s)) / s;
      }
    }

    std::vector<int> pk;
    for (int i = 1; i < ngrid - 1; ++i)
    {
      if (out.density[i] <= out.density[i - 1] || out.density[i] < out.density[i + 1]) continue;
      double minL = out.density[i];
      for (int j = i - 1; j >= 0; --j)
      {
        if (out.density[j] > out.density[i]) break;
        minL = std::min(minL, out.density[j]);
      }
      double minR = out.density[i];
      for (int j = i + 1; j < ngrid; ++j)
      {
        if (out.density[j] > out.density[i]) break;
        minR = std::min(minR, out.density[j]);
      }
      if (out.density[i] - std::max(minL, minR) > promfrac * out.density[i]) pk.push_back(i);
    }

    for (int i : pk) out.peakR.push_back(out.grid[i]);

    for (size_t ir = 0; ir + 1 < pk.size(); ++ir)
    {
      int imin = pk[ir];
      for (int b = pk[ir]; b <= pk[ir + 1]; ++b)
        if (out.density[b] < out.density[imin]) imin = b;
      out.bound.push_back(out.grid[imin]);
    }
    return out;
  }

  int rowOf(double R, const std::vector<double> &bound)
  {
    return std::upper_bound(bound.begin(), bound.end(), R) - bound.begin();
  }

  std::vector<double> refine(const std::vector<Cluster> &cl, const Result &res)
  {
    const size_t n = res.peakR.size();
    std::vector<double> sw(n, 0.0), swr(n, 0.0);
    for (const auto &c : cl)
    {
      if (c.R < res.grid.front() || c.R > res.grid.back()) continue;
      const int ir = rowOf(c.R, res.bound);
      sw[ir] += c.adc;
      swr[ir] += c.adc * c.R;
    }
    std::vector<double> R(n, 0.0);
    for (size_t i = 0; i < n; ++i) R[i] = (sw[i] > 0.0) ? swr[i] / sw[i] : res.peakR[i];
    return R;
  }

  // Folds a phi value to be between lo and lo+petal
  double foldPetal(double phi, double lo, double petal)
  {
    double p = fmod(phi - lo, petal);
    if (p < 0.0) p += petal;
    return lo + p;
  }

  double circDiff(double a, double b, double petal)
  {
    double d = a - b;
    while (d > 0.5 * petal) d -= petal;
    while (d < -0.5 * petal) d += petal;
    return d;
  }

  struct Stripes
  {
    std::vector<double> grid;
    std::vector<double> density;
    std::vector<double> peakPhi;
  };

  std::map<int, TruthRowPattern> buildTruthRowPattern(CDBTTree &cdbttree, int side)
  {
    std::map<int, TruthRowPattern> truthRows;
    for(int row=0; row<32; row++)
    {
      std::vector<std::pair<int,double>> entries;
      for (int iphi = 0; iphi < 12; iphi++)
      {
        unsigned int truthIndex = (side ? 18 : 0)*10000 + (row*100) + iphi;
        double phiVal = cdbttree.GetDoubleValue(truthIndex, "truthPhi");
        if (std::isnan(phiVal)) continue;
        entries.push_back({iphi, phiVal});
      }
      std::sort(entries.begin(), entries.end(), [](auto &a, auto &b){ return a.second < b.second; });

      TruthRowPattern truthRow;
      for (auto &[iphi, phiVal] : entries)
      {
        truthRow.stripePhi.push_back(phiVal);
        truthRow.stripeIPhi.push_back(iphi);
      }
      truthRows[row] = truthRow;
    }
    return truthRows;
  }

  // tuning knobs -- identical meaning to the macro's local consts
  constexpr double sigR = 0.10;
  constexpr double sigGap = 0.15;
  constexpr double stripeEdgeFloorFrac = 0.25;
  constexpr double stripePromFrac = 0.05;
  constexpr double lamSearchGate = 0.04;
  constexpr int lamRowHalfWindow = 3;
  constexpr double lamMatchTol = 0.010;
  constexpr double lamMatchFrac = 0.7;
  constexpr int lamRecoverySearchRows = 8;
  constexpr double lamRecoveryMatchTol = 0.005;
  constexpr double lamRecoveryFloorFrac = 0.10;
  const double phiPetalLo[2] = {-TMath::Pi() / 18.0, 0.0};

  //! Find stripe peaks in one row's folded-phi KDE.
  //! Seeds from the tallest intetior peak and walks outward in steps of ~recoDphi
  //! to find the rest, then extends past the walk's ends by scanning for any
  //! remaining local peaks (isLocalPeak) above a height floor set by the
  //! already-found peaks (medPeakHeight * edgeFloorFrac).
  Stripes findStripes(const std::vector<Cluster> &clusters, double lo,
                      double sigma_ref = 0.004,   // rad, median-ADC cluster
                      double sigma_min = 0.0015,  // rad, ~half a pad pitch
                      int ngrid = 720, double petal = kPetalDefault,
                      double recoDphi = kPetalDefault / 11,
                      double edgeFloorFrac = 0.5)
  {
    Stripes out;
    if(clusters.empty()) return out;

    std::vector<double> adcs;
    adcs.reserve(clusters.size());
    for(const auto &c : clusters) { adcs.push_back(c.adc); }
    std::nth_element(adcs.begin(), adcs.begin() + adcs.size() / 2, adcs.end());
    double adc_ref = std::max(adcs[adcs.size() / 2], 1.0);

    //! set up grid and density for Gaussian KDE
    double dp = petal / ngrid;
    out.grid.resize(ngrid);
    out.density.assign(ngrid, 0.0);
    for(int b=0; b<ngrid; b++)
    {
      out.grid[b] = lo + (b+0.5) * dp;
    }

    for(const auto &c : clusters)
    {
      const double p = foldPetal(c.phi, lo, petal);
      //! Set width of Gaussian for KDE.
      //! Higher ADC clusters have lower width
      //! Bounded by sigma_min and petal/8
      double s = sigma_ref * sqrt(adc_ref / std::max(c.adc, 1.0));
      s = std::max(s, sigma_min);
      s = std::min(s, petal / 8.0);

      //! c0 = grid point closest to cluster phi
      const int c0 = (int) ((p - lo) / dp);
      //! half -> convert 4 standard deviations into a number of bins
      const int half = (int) ceil(4.0 * s / dp);
      for(int k=-half; k<=half; k++)
      {
        const int b = (((c0 + k) % ngrid) + ngrid) % ngrid; // double modulo to wrap around the petal boundary
        const double d = circDiff(out.grid[b], p, petal); // signed shortest distance between p and grid point on circle (dealing with wrapping)
        out.density[b] += exp(-0.5 * d * d / (s * s)) / s;
      }
    }

    //! easy accessors for grid coordinates and density values
    auto dens = [&](int i){ return out.density[((i % ngrid) + ngrid) % ngrid]; };

    //! seed: tallest peak in the central half of the row (away from the edges)
    int maxCand = -1;
    double maxCont = 0.0;
    for(int i=(int)floor(ngrid/4.0); i<(int)ceil(3.0*ngrid/4.0); i++)
    {
      double cont = dens(i);
      if(cont > maxCont)
      {
        maxCont = cont;
        maxCand = i;
      }
    }

    std::vector<int> cand;
    cand.push_back(maxCand);

    //! walk outward in steps of ~recoDphi, searching a window around each step
    int searchStep = (int) round(recoDphi / dp);
    int minSearch = (int) floor(0.8 * recoDphi / dp);
    int maxSearch = (int) ceil(1.2 * recoDphi / dp);
    int searchWidth = std::max(maxSearch - searchStep, searchStep - minSearch);

    int nLoop = 0;
    //! Starting from last peak in phi, search for new peaks at larger phi
    //! by moving one recoDphi step and searching in a window around that
    while(true)
    {
      int prevIdx = cand.back();
      int j0 = prevIdx + searchStep - searchWidth;
      int j1 = std::min(prevIdx + searchStep + searchWidth, ngrid-1);
      if(j0 > j1) break;

      int localCand = -1;
      double localMax = 0.0;
      for(int j=j0; j<=j1; j++)
      {
        double cont = dens(j);
        if(cont > localMax)
        {
          localMax = cont;
          localCand = j;
        }
      }
      if(localCand < 0 || nLoop > 20) break;
      cand.push_back(localCand);
      nLoop++;
    }

    nLoop = 0;
    //! Starting from first peak in phi, search for new peaks at smaller phi
    //! by moving one recoDphi step and searching in a window around that
    while(true)
    {
      int prevIdx = cand.front();
      int j0 = std::max(prevIdx - searchStep - searchWidth, 0);
      int j1 = prevIdx - searchStep + searchWidth;
      if(j0 > j1) break;

      int localCand = -1;
      double localMax = 0.0;
      for(int j=j0; j<=j1; j++)
      {
        double cont = dens(j);
        if(cont > localMax)
        {
          localMax = cont;
          localCand = j;
        }
      }
      if(localCand < 0 || nLoop > 20) break;
      cand.insert(cand.begin(), localCand);
      nLoop++;
    }

    std::sort(cand.begin(), cand.end());

    //! Check first and last peaks to ensure they don't continue rising on either side (peak split by petal wrapping)
    const int checkSamples = 8;  // samples checked on each side to confirm a local peak
    auto isLocalPeak = [&](int idx)
    {
      for (int k = 1; k <= checkSamples; ++k)
      {
        if (dens(idx - k) > dens(idx - k + 1))
        {
          return false;
        }
        if (dens(idx + k) > dens(idx + k - 1))
        {
          return false;
        }
      }
      return true;
    };

    //! Iteratively remove first or last peak if it is rising toward wrap boundary
    while (cand.size() > 1 && !isLocalPeak(cand.front())) cand.erase(cand.begin());
    while (cand.size() > 1 && !isLocalPeak(cand.back())) cand.pop_back();

    double medPeakHeight = 0.0;
    if(!cand.empty())
    {
      std::vector<double> candH;
      for(const auto &c : cand)
      {
        candH.push_back(dens(c));
      }
      std::nth_element(candH.begin(), candH.begin() + candH.size() / 2, candH.end());
      medPeakHeight = candH[candH.size() / 2];
    }

    //! Find the single tallest peak in range [i0, i1]
    auto scanForPeak = [&](int i0, int i1) -> int   // inclusive range, i0 <= i1
    {
      int best = -1;
      double bestH = 0.0;
      for(int i=i0; i<=i1; i++)
      {
        if(!isLocalPeak(i)) continue;
        if(dens(i) > bestH)
        {
          bestH = dens(i);
          best = i;
        }
      }
      return best;
    };

    //! Recursively scan [i0,i1] for peaks, not just the tallest one.
    //! Start with highest peak in range, then split range around that peak and exclude it
    //! to recursively identify all peaks above threshold
    //! Used to identify lamination sitting next to a taller ordinary stripe in the same
    //! region that would would otherwise be missed
    auto scanRangeForPeaks = [&](int i0, int i1, bool prepend)
    {
      std::vector<std::pair<int,int>> ranges = { {i0,i1} };
      std::vector<int> found;

      while(!ranges.empty())
      {
        auto [a, b] = ranges.back();
        ranges.pop_back();
        if(a > b) continue;

        int best = scanForPeak(a, b);
        if(best < 0 || dens(best) <= edgeFloorFrac * medPeakHeight) continue;

        found.push_back(best);
        ranges.push_back({a, best - 1});
        ranges.push_back({best + 1, b});
      }

      std::sort(found.begin(), found.end());
      if(prepend) cand.insert(cand.begin(), found.begin(), found.end());
      else cand.insert(cand.end(), found.begin(), found.end());
    };

    //! Search between first grid point and first peak and the last peak and last grid point
    //! for any additional peaks that may have been missed (searching for laminations).
    if(!cand.empty())
    {
      if(cand.front() > 0) scanRangeForPeaks(0, cand.front() - 1, true);
      if(cand.back() < ngrid - 1) scanRangeForPeaks(cand.back() + 1, ngrid - 1, false);
    }


    for(int i : cand)
    {
      out.peakPhi.push_back(out.grid[i]);
    }
    std::sort(out.peakPhi.begin(), out.peakPhi.end());
    return out;
  } // end findStripes

  //! Find the peak within searchGate of "edge" whose position is most stable across a window
  //! of neighbouring rows -- matches some peak in at least matchFrac of the rows within matchTol.
  //! A lamination sits at a near-fixed phi, so it's stable across many rows; a real stripe's
  //! position shifts for each row, so a single neighbour can spuriously match but a wide
  //! window won't.
  int findLaminationNear(const std::vector<std::vector<double>> &peaksPerRow, int rowIdx, double edge, double searchGate,
                          int rowHalfWindow, double matchTol, double matchFrac, double petal)
  {
    const auto &peaks = peaksPerRow[rowIdx];
    const int nrow = (int) peaksPerRow.size();

    int bestIdx = -1;
    int bestMatches = -1;
    for(size_t idx=0; idx<peaks.size(); idx++)
    {
      if(fabs(circDiff(peaks[idx], edge, petal)) > searchGate) continue;

      int nMatched = 0, nChecked = 0;
      for(int j=std::max(0, rowIdx - rowHalfWindow);
          j<=std::min(nrow - 1, rowIdx + rowHalfWindow);
          j++)
      {
        if(j == rowIdx || peaksPerRow[j].empty()) continue;
        nChecked++;
        double best = 1e9;
        for(double q : peaksPerRow[j])
        {
          best = std::min(best, fabs(circDiff(peaks[idx], q, petal)));
        }
        if(best < matchTol) nMatched++;
      }
      if(nChecked == 0) continue;
      if((double) nMatched / nChecked >= matchFrac && nMatched > bestMatches)
      {
        bestMatches = nMatched;
        bestIdx = (int) idx;
      }
    }
    return bestIdx;
  }       
  
  //! Snapshot of what findLaminationNear already confirmed for one edge, used
  //! as prediction sources for the recovery pass below. Captured BEFORE
  //! recovery runs and never updated during it, so recovery for one row can
  //! never depend on another row's recovered (rather than independently
  //! confirmed) result -- keeps recovery non-cascading, so uncertainty can't
  //! compound across rows.
  struct LamEdgeSnapshot
  {
    std::vector<bool> confirmed;
    std::vector<double> pos;  // valid only where confirmed is true
  };

  LamEdgeSnapshot snapshotLamEdge(const std::vector<std::vector<double>> &peaksPerRow, const std::vector<int> &edgeIdx)
  {
    LamEdgeSnapshot snap;
    int nrow = (int) peaksPerRow.size();
    snap.confirmed.assign(nrow, false);
    snap.pos.assign(nrow, 0.0);
    for(int i=0; i<nrow; i++)
    {
      if(edgeIdx[i] < 0) continue;
      snap.confirmed[i] = true;
      snap.pos[i] = peaksPerRow[i][edgeIdx[i]];
    }
    return snap;
  }

  //! Predict row i's lamination position for one edge from the nearest ORIGINALLY-confirmed
  //! neighbouring rows (within maxSearchRows on each side): linear interpolation if confirmed
  //! rows exist on both sides, extrapolation (nearest value) if only one side has data.
  std::pair<double, bool> predictLamPos(const LamEdgeSnapshot &snap, int row, int maxSearchRows, double petal)
  {
    int nrow = (int) snap.confirmed.size();
    int below = -1, above = -1;
    for(int j=row-1; j>=std::max(0, row-maxSearchRows); j--)
    {
      if(snap.confirmed[j])
      {
        below = j;
        break;
      }
    }
    for(int j=row+1; j<=std::min(nrow-1, row+maxSearchRows); j++)
    {
      if(snap.confirmed[j])
      {
        above = j;
        break;
      }
    }

    if(below >= 0 && above >= 0)
    {
      double frac = double(row - below) / (double)(above - below);
      return { snap.pos[below] + frac * circDiff(snap.pos[above], snap.pos[below], petal), true};
    }
    if(below >= 0) return {snap.pos[below], true};
    if(above >= 0) return {snap.pos[above], true};
    return {0.0, false};
  }

  //! median density at this row's oen existing peak positions -- used as a 
  //! row-local scale for the recovery floor below.
  double medianPeakDensity(const Stripes &stripeRow, const std::vector<double> &peaks)
  {
    if(peaks.empty() || stripeRow.grid.size() < 2) return 0.0;
    const double dp = stripeRow.grid[1] - stripeRow.grid[0];
    std::vector<double> h;
    for(double p : peaks)
    {
      int idx = (int) std::round((p - stripeRow.grid[0]) / dp);
      idx = std::max(0, std::min((int) stripeRow.grid.size() - 1, idx));
      h.push_back(stripeRow.density[idx]);
    }
    std::nth_element(h.begin(), h.begin() + h.size() / 2, h.end());
    return h[h.size() / 2];
  }

  //! Recover a lamination for one row/edge that fineLaminationNear missed, using a position predicted from
  //! independently-confirmed neighbouring rows.
  //! (a) checks whether an EXISTING peak is already close to the prediction -- recovers a peak that WAS found
  //! but failed the cross-row stability test.
  //! (b) if no existing peak is close, scans the row's own raw LDE near the prediction for the tallest bin and,
  //! if it clears absFloor, inserts it as a NEW peak -- recovers a lamination that was fully abosorbed into
  //! a taller neighbour and never became a distinct peak
  bool recoverLamination(std::vector<double> &peaks, std::vector<bool> &isLam, const Stripes &stripeRow,
                          double predicted, double matchTol, double petal, double absFloor)
  {
    for(size_t k=0; k<peaks.size(); k++)
    {
      if(fabs(circDiff(peaks[k], predicted, petal)) < matchTol)
      {
        isLam[k] = true;
        return true;
      }
    }

    if(stripeRow.grid.size() < 2) return false;
    const double dp = stripeRow.grid[1] - stripeRow.grid[0];
    const int half = std::max(1, (int) (matchTol / dp));
    const int ngrid = (int) stripeRow.grid.size();

    int c0 = -1;
    double bestD = 1e300;
    for(int b=0; b<ngrid; b++)
    {
      double d = fabs(circDiff(stripeRow.grid[b], predicted, petal));
      if(d < bestD)
      {
        bestD = d;
        c0 = b;
      }
    }
    if(c0 < 0) return false;

    int bestBin = -1;
    double bestH = 0.0;
    for(int k=-half; k<=half; k++)
    {
      int b = ((c0 + k) % ngrid + ngrid) % ngrid;
      if(stripeRow.density[b] > bestH)
      {
        bestH = stripeRow.density[b];
        bestBin = b;
      }
    }
    if(bestBin < 0 || bestH < absFloor) return false;

    double newPeakPos = stripeRow.grid[bestBin];
    size_t insertAt = std::lower_bound(peaks.begin(), peaks.end(), newPeakPos) - peaks.begin();
    peaks.insert(peaks.begin() + insertAt, newPeakPos);
    isLam.insert(isLam.begin() + insertAt, true);
    return true;
  }
  
  struct Cost
  {
    double nMismatch;
    double gapChi2;
    double rChi2;
    bool operator<(const Cost &o) const
    {
      if (nMismatch != o.nMismatch) return nMismatch < o.nMismatch;
      if (gapChi2 != o.gapChi2) return gapChi2 < o.gapChi2;
      return rChi2 < o.rChi2;
    }
    Cost operator+(const Cost &o) const
    {
      return {nMismatch + o.nMismatch, gapChi2 + o.gapChi2, rChi2 + o.rChi2};
    }
  };

  const Cost INF = {1e300, 1e300, 1e300};

  struct MatchResult
  {
    std::vector<int> truthOf;
    std::vector<int> skipped;
    Cost dp{1e300, 1e300, 1e300};
  };

  Cost pairCost(double recoR, int recoN, const TruthRow &t, double sigmaR)
  {
    Cost c{0.0, 0.0, 0.0};
    if(recoN > 0 && t.nstripes > 0)
    {
      double d = recoN - t.nstripes;
      c.nMismatch = d * d;
    }
    double dR = recoR - t.R;
    c.rChi2 = dR * dR / (sigmaR * sigmaR);
    return c;
  }

  //! Cost of the reco row-to-row R gap disagreeing with the corresponding truth gap. Uses actual truth R centers
  //! so a skip spanning a module boundary is judged against the true, larger gap rather than needing
  //! special-case handling.
  inline double gapCost(double recoR_prev, double recoR_curr, double truthR_prev, double truthR_curr, double sigmaGap)
  {
    double recoGap = recoR_curr - recoR_prev;
    double truthGap = truthR_curr - truthR_prev;
    double d = recoGap - truthGap;
    return d * d / (sigmaGap * sigmaGap);
  }

  //! dp[i][j]: best Cost matching reco rows 0..i, with row i matched to truth row j exactly (needed so gapCost knows
  //! the previous match's truth row)
  void matchDP(const std::vector<double> &recoR, const std::vector<int> &recoN, const std::vector<TruthRow> &truth, double sigmaR, double sigmaGap, MatchResult &best, int verbosity)
  {
    const int nr = recoR.size();
    const int nt = truth.size();
    best = MatchResult();
    if(nr == 0 || nt == 0 || nr > nt)
    {
      if(verbosity) std::cout << "rowmatch: " << nr << " reco rows vs " << nt << " truth rows -- cannot align" << std::endl;
      return;
    }

    std::vector<std::vector<Cost>> dp(nr, std::vector<Cost>(nt, INF));
    std::vector<std::vector<int>> prevJ(nr, std::vector<int>(nt, -1));

    for(int j=0; j<nt; j++)
    {
      dp[0][j] = pairCost(recoR[0], recoN[0], truth[j], sigmaR);  // no gap cost for the first match
    }

    for(int i=1; i<nr; i++)
    {
      for(int j=i; j<nt; j++)
      {
        Cost matchP = pairCost(recoR[i], recoN[i], truth[j], sigmaR);
        for(int jp = i - 1; jp < j; jp++)
        {
          if(!(dp[i-1][jp] < INF)) continue;
          double gc = gapCost(recoR[i-1], recoR[i], truth[jp].R, truth[j].R, sigmaGap);
          Cost cand = dp[i-1][jp] + Cost{matchP.nMismatch, gc, matchP.rChi2};
          if(cand < dp[i][j])
          {
            dp[i][j] = cand;
            prevJ[i][j] = jp;
          }
        }
      }
    }

    int bestJ = -1;
    Cost bestCost = INF;
    for(int j=nr-1; j<nt; j++)
    {
      if(dp[nr-1][j] < bestCost)
      {
        bestCost = dp[nr-1][j];
        bestJ = j;
      }
    }

    if(bestJ < 0) return;

    std::vector<int> truthOf(nr);
    int i = nr - 1, j = bestJ;
    while(i>=0)
    {
      truthOf[i] = j;
      int pj = (i > 0 ? prevJ[i][j] : -1);
      i--;
      j = pj;
    }

    std::vector<bool> used(nt, false);
    for(int t : truthOf) used[t] = true;
    std::vector<int> skipped;
    for(int t=0; t<nt; t++)
    {
      if(!used[t]) skipped.push_back(t);
    }

    best = {truthOf, skipped, bestCost};
    return;
  }

  //! Rebase a row's real (non-lamination) peaks so the lamination sits at `lo`, wrap into [lo, lo+petal), then match against
  //! truthRowPattern's pattern if the count and spacing agree. Returns one entry per input peak (lamination and unmatched
  //! entries get -1); origIdx tracks which output slot each repased peak came from.
  std::vector<int> alignPhiPeaksToTruth(const std::vector<double> &peaks, const std::vector<bool> &isLam, double dPhi, double lo, double petal, const TruthRowPattern &truthRowPattern)
  {
    std::vector<double> peaksToUse;
    std::vector<int> origIdx;
    double lamPhi = 0.0;
    bool haveLam = false;
    for(size_t i=0; i<isLam.size(); i++)
    {
      if(isLam[i])
      {
        lamPhi = peaks[i];
        haveLam = true;
        continue;
      }
      peaksToUse.push_back(peaks[i]);
      origIdx.push_back((int) i);
    }

    std::vector<int> iphiAssignment(peaks.size(), -1);
    if(!haveLam || peaksToUse.empty()) return iphiAssignment;

    for(auto &p : peaksToUse)
    {
      p = p - lamPhi + lo;
      while(p < lo) p += petal;
      while(p >= lo + petal) p -= petal;
    }

    std::vector<std::size_t> phiSortIdx(peaksToUse.size());
    std::iota(phiSortIdx.begin(), phiSortIdx.end(), 0);

    std::sort(phiSortIdx.begin(), phiSortIdx.end(), [&](std::size_t a, std::size_t b){ return peaksToUse[a] < peaksToUse[b]; });

    if(peaksToUse.size() == truthRowPattern.stripePhi.size())
    {
      bool allGood = true;
      for(size_t p = 0; p+1 < peaksToUse.size(); p++)
      {
        double gap = peaksToUse[phiSortIdx[p+1]] - peaksToUse[phiSortIdx[p]];
        if(gap > 1.2 * dPhi || gap < 0.8 * dPhi)
        {
          allGood = false;
          break;
        }
      }
      if(allGood)
      {
        for(size_t p = 0; p < peaksToUse.size(); p++)
        {
          iphiAssignment[origIdx[phiSortIdx[p]]] = truthRowPattern.stripeIPhi[p];
        }
      }
    }
    return iphiAssignment;
  }
  
  struct PetalCand
  {
    double dPhi;
    int petal;    
  };

  int nearestPeak(double foldedPhi, const std::vector<double> &peaks, double petal)
  {
    int best = 0;
    double bestD = 1e300;
    for(size_t k = 0; k<peaks.size(); k++)
    {
      double d = fabs(circDiff(foldedPhi, peaks[k], petal));
      if(d < bestD)
      {
        bestD = d;
        best = (int) k;
      }
    }
    return best;
  }

  unsigned int encodeTruthIndex(int side, int petalIdx, int row, int stripeIdx)
  {
    int pp = (side == 0) ? petalIdx : petalIdx + 18;
    return (unsigned int) (pp*10000 + row*100 + stripeIdx);
  }

  std::map<TrkrDefs::cluskey, std::pair<int,bool>> matchSide(const std::vector<Cluster>& clusters, CDBTTree &cdbttree, const std::map<int, TruthRowPattern> &truthRowPatterns, const std::vector<TruthRow> &truth, int side, std::string QABase, int verbosity)
  {
    const char *sname = side ? "North" : "South";
    const double lo = phiPetalLo[side];
    const double petal = kPetalDefault;

    std::map<TrkrDefs::cluskey, std::pair<int,bool>> keyMatch;
    for(const auto &c : clusters)
    {
      keyMatch[c.key] = {-999,false};
    }

    if (clusters.empty())
    {
      std::cout << "LaserClusterTruthMatcher - matchSide: no clusters found for side " << sname << std::endl;
      return keyMatch;
    }

    // ---- rows in R ----
    Result rowResult = findRows(clusters, 0.15, 0.05);
    std::vector<double> rowR = refine(clusters, rowResult);
    const int nrow = (int)rowR.size();
    if (nrow == 0)
    {
      std::cout << "LaserClusterTruthMatcher: found 0 rows for side " << sname << std::endl;
      return keyMatch;
    }

    // ---- assign clusters to rows ----
    std::vector<std::vector<Cluster>> rowClus(nrow);
    for (const auto &c : clusters)
    {
      rowClus[rowOf(c.R, rowResult.bound)].push_back(c);
    }

    if(!QABase.empty())
    {
      auto palette = TColor::GetPalette();
      const Int_t nColors = palette.GetSize();

      std::vector<TH1D *> hRADCrow(nrow);
      TCanvas *c1 = new TCanvas();
      gStyle->SetOptStat(0);
      TH1D *hRADC = new TH1D(std::format("hRADC_{}", sname).c_str(),
                        std::format("R of {} ADC-Weighted Aggregated Clusters;R [cm]",
                                    sname).c_str(),
                        400, 28, 78);      
      for(int i=0; i<nrow; i++)
      {
        hRADCrow[i] = new TH1D(std::format("hRADC_{}_row{}", sname, i).c_str(),
                        std::format("R of {} ADC-Weighted Aggregated Clusters;R [cm]",
                                    sname).c_str(),
                        400, 28, 78);
        const Int_t code = palette.At(i * (nColors / std::max(nrow, 1)));
        hRADCrow[i]->SetLineColor(code);
        hRADCrow[i]->SetFillColor(code);

        for(const auto& c : rowClus[i])
        {
          hRADCrow[i]->Fill(c.R, c.adc);
        }

        hRADC->Add(hRADCrow[i]);
      }

      hRADC->Draw();
      for(int i=0; i<nrow; i++) hRADCrow[i]->Draw("HIST SAME");


      for (double b : rowResult.bound)
      {
        TLine *l = new TLine(b, 0, b, hRADC->GetMaximum());
        l->SetLineColor(kRed);
        l->Draw("same");
      }
      c1->SaveAs(std::format("{}_RPeaks_{}.pdf",QABase, sname).c_str());
    }

    // ---- pass 1: stripe spacing per row from the neighbour-gap spectrum ----

    std::vector<TH1D *> hDPhi(nrow);
    std::vector<double> recoDphi(nrow, -1.0);

    for(int i=0; i<nrow; i++)
    {
      hDPhi[i] = new TH1D(std::format("hDPhi_{}_row{}",sname,i).c_str(),
                          std::format("#Delta#phi, {} row {};#Delta#phi [rad]",sname, i).c_str(),
                          400, 0, petal);

      std::sort(rowClus[i].begin(), rowClus[i].end(),
                [](const Cluster &a, const Cluster &b){ return a.phi < b.phi; });
      
      for(int j=0; j<(int)rowClus[i].size(); j++)
      {
        const double dPhi = (j == (int) rowClus[i].size() - 1)
          ? rowClus[i][0].phi + 2*TMath::Pi() - rowClus[i][j].phi // wrap around in 2pi
          : rowClus[i][j+1].phi - rowClus[i][j].phi;
        hDPhi[i]->Fill(dPhi);
      }
      if(hDPhi[i]->GetEntries() > 20)
      {
        recoDphi[i] = hDPhi[i]->GetBinCenter(hDPhi[i]->GetMaximumBin());
      }
    }

    if(!QABase.empty())
    {
      TCanvas *c1 = new TCanvas();
      c1->SaveAs(std::format("{}_dPhi_{}.pdf[",QABase, sname).c_str());
      for(int i=0; i<nrow; i++)
      {
        c1->Clear();
        hDPhi[i]->Draw();
        c1->SaveAs(std::format("{}_dPhi_{}.pdf",QABase, sname).c_str());
      }
      c1->SaveAs(std::format("{}_dPhi_{}.pdf]",QABase, sname).c_str());
    }

    // ---- pass 2a: stripe peak positions per row, stepping by dPhi ----
    std::vector<Stripes> stripes(nrow);
    std::vector<std::vector<double>> peaksPerRow(nrow);

    for(int i=0; i<nrow; i++)
    {
      stripes[i] = findStripes(rowClus[i], lo, 0.004, 0.0015, 720,
                                petal, recoDphi[i], stripeEdgeFloorFrac);
      peaksPerRow[i] = stripes[i].peakPhi;                            
    }

    // ---- pass 2b: classify each row's first/last peaks as lamination or stripe ----
    std::vector<std::vector<bool>> isLamRow(nrow);
    std::vector<int> recoN(nrow, 0);
    std::vector<int> lamFrontIdxPerRow(nrow, -1), lamBackIdxPerRow(nrow, -1);

    for(int i=0; i<nrow; i++)
    {
      isLamRow[i] = std::vector<bool>(peaksPerRow[i].size(), false);

      int lamFrontIdx = findLaminationNear(peaksPerRow, i, lo, lamSearchGate, lamRowHalfWindow,
                                            lamMatchTol, lamMatchFrac, petal);
      int lamBackIdx = findLaminationNear(peaksPerRow, i, lo+petal, lamSearchGate, lamRowHalfWindow,
                                            lamMatchTol, lamMatchFrac, petal);  
                                            
      if (lamFrontIdx >= 0) isLamRow[i][lamFrontIdx] = true;
      if (lamBackIdx >= 0 && lamBackIdx != lamFrontIdx) isLamRow[i][lamBackIdx] = true;
      
      lamFrontIdxPerRow[i] = lamFrontIdx;
      lamBackIdxPerRow[i]  = lamBackIdx;

      recoN[i] = (int) std::count(isLamRow[i].begin(), isLamRow[i].end(), false);

      if(verbosity > 1)
      {
        std::cout << "  row " << std::setw(2) << i
                << "   R = "      << std::setw(8) << rowR[i]
                << "   clusters " << std::setw(5) << (int) rowClus[i].size()
                << "   spacing "  << std::setw(10) << recoDphi[i]
                << "   peaks "    << std::setw(3) << (int) peaksPerRow[i].size()
                << "   stripes "  << std::setw(3) << recoN[i]
                << "   lam "      << ((int) peaksPerRow[i].size() - recoN[i])
                << std::endl;
      }
    }

    // ---- pass 2c: recover laminations findLaminationNear missed, using
    // positions predicted from independently-confirmed neighbouring rows ----
    auto snapFront = snapshotLamEdge(peaksPerRow, lamFrontIdxPerRow);
    auto snapBack = snapshotLamEdge(peaksPerRow, lamBackIdxPerRow);

    for(int i=0; i<nrow; i++)
    {
      if(lamFrontIdxPerRow[i] < 0)
      {
        auto [pred, have] = predictLamPos(snapFront, i, lamRecoverySearchRows, petal);
        if(have)
        {
          double floor = lamRecoveryFloorFrac * medianPeakDensity(stripes[i], peaksPerRow[i]);
          recoverLamination(peaksPerRow[i], isLamRow[i], stripes[i], pred,
                                lamRecoveryMatchTol, petal, floor);
          if(verbosity)
          {
            std::cout << "  row " << i << " front lamination recovered near " << pred << std::endl;
          }                                
        }
      }
      if(lamBackIdxPerRow[i] < 0)
      {
        auto [pred, have] = predictLamPos(snapBack, i, lamRecoverySearchRows, petal);
        if(have)
        {
          double floor = lamRecoveryFloorFrac * medianPeakDensity(stripes[i], peaksPerRow[i]);
          recoverLamination(peaksPerRow[i], isLamRow[i], stripes[i], pred,
                                lamRecoveryMatchTol, petal, floor);
          if(verbosity)
          {
            std::cout << "  row " << i << " back lamination recovered near " << pred << std::endl;
          }                                
        }
      }
      recoN[i] = std::count(isLamRow[i].begin(), isLamRow[i].end(), false);
    }

    if(!QABase.empty())
    {
      TCanvas *c1 = new TCanvas();
      gStyle->SetOptStat(0);
      c1->SaveAs(std::format("{}_phiPeaks_{}.pdf[",QABase, sname).c_str());
      for(int i=0; i<nrow; i++)
      {
        TGraph *gr = new TGraph(stripes[i].grid.size(), &stripes[i].grid[0], &stripes[i].density[0]);
        gr->SetTitle(std::format("{} row {} (R = {:.2f}, spacing = {:.4f}, stripes = {});"
                              "folded #phi [rad];density",
                              sname, i, rowR[i], recoDphi[i], recoN[i]).c_str());
        gr->GetXaxis()->SetRangeUser(lo, lo+petal);
        gr->SetMarkerStyle(20);
        gr->SetMarkerSize(0.5);
        gr->Draw("ALP");
        for (size_t ip = 0; ip < peaksPerRow[i].size(); ip++)
        {
          const double p = peaksPerRow[i][ip];
          TLine *l = new TLine(p, 0, p, gr->GetHistogram()->GetMaximum());
          l->SetLineColor(isLamRow[i][ip] ? kRed : kBlue);
          l->Draw("same");
        }
        c1->SaveAs(std::format("{}_phiPeaks_{}.pdf",QABase, sname).c_str());
      }
      c1->SaveAs(std::format("{}_phiPeaks_{}.pdf]",QABase, sname).c_str());
    }

    // ---- step 3: match to truth R pattern using number of stripes in each row (uses gaps and dR as backups) ----
    if(truth.empty())
    {
      if(verbosity) std::cout << "no truth table for " << sname << ", skipping matching" << std::endl;
      return keyMatch;
    }

    MatchResult best;
    matchDP(rowR, recoN, truth, sigR, sigGap, best, verbosity);
    if(best.truthOf.empty())
    {
      if(verbosity) std::cout << "no good matching found for " << sname << ", returning without matches" << std::endl;
      return keyMatch;
    }

    if(verbosity > 1)
    {
      std::cout << "--- " << sname << " truth matching ---" << std::endl;
      std::cout << "  best nMismatch  " << best.dp.nMismatch << "  (primary)" << std::endl;
      std::cout << "  best gapChi2    " << best.dp.gapChi2   << "  (spacing tiebreak)" << std::endl;
      std::cout << "  best rChi2      " << best.dp.rChi2     << "  (last resort)" << std::endl;

      std::cout << "  unmatched truth rows:";
      if (best.skipped.empty()) std::cout << " none";
      for (int t : best.skipped)
      {
        std::cout << " " << t << " (R=" << truth[t].R << ")";
      }
      std::cout << std::endl;
    }

    if(verbosity > 2)
    {
      for (int i = 0; i < nrow; i++)
      {
        const TruthRow &t = truth[best.truthOf[i]];
        const double dR = rowR[i] - t.R;
        std::cout << "  reco " << std::setw(2) << i
                  << " -> truth " << std::setw(2) << best.truthOf[i]
                  << "   dR = " << std::setw(10) << dR;
        if (recoDphi[i] > 0.0 && t.dphi > 0.0)
        {
          std::cout << "   dphi ratio = " << recoDphi[i] / t.dphi;
        }
        if (recoN[i] > 0 && t.nstripes > 0)
        {
          std::cout << "   stripes " << recoN[i] << " vs " << t.nstripes;
        }
        std::cout << std::endl;
      }
    }
    
    // ---- step 4: truth iphi assignment for each peak in each row ----
    std::vector<std::vector<int>> iphiAssignments;
    int nRowsFullyAligned = 0, nRowsWithRealPeaks = 0;
    for(int i=0; i<nrow; i++)
    {
      std::vector<int> iphiAssignment = alignPhiPeaksToTruth(peaksPerRow[i], isLamRow[i], recoDphi[i],
                                                              lo, petal, truthRowPatterns.at(best.truthOf[i]));

      bool haveRealPeak = false, allAssigned = true;
      for(size_t k=0; k<isLamRow[i].size(); k++)
      {
        if(isLamRow[i][k]) continue;
        haveRealPeak = true;
        if(iphiAssignment[k] < 0) allAssigned = false;
      }          
      if(haveRealPeak)
      {
        nRowsWithRealPeaks++;
        if(allAssigned) nRowsFullyAligned++;
      }
      
      if(verbosity > 1)
      {
        std::cout << "   row " << std::setw(3) << i << " " << std::setw(8) << "peak" << std::setw(12) << "assignment" << std::endl;
        for (int j = 0; j < (int) peaksPerRow[i].size(); j++)
        {
          std::cout << "           " << std::setw(8) << peaksPerRow[i][j] << std::setw(12) << iphiAssignment[j] << std::endl;
        }
      }
      iphiAssignments.push_back(iphiAssignment);
    }

    if(verbosity) std::cout << "  phi alignment: " << nRowsFullyAligned << "/" << nRowsWithRealPeaks << " rows fully assigned" << std::endl;

    // ---- step 5a: per cluster truth assignment ----
    //! For every cluster, find its (row, iphi); lamination never get a phi-slot and never enert
    //! the petal search at all. For everything else, rank ALL 18 petal candidated by phi distance
    //! (not just the closest)
    std::vector<int> clusterRow(clusters.size(), -1);
    std::vector<int> clusterIphi(clusters.size(), -1);
    std::vector<bool> clusterIsLam(clusters.size(), false);
    std::vector<std::vector<PetalCand>> clusterPetalCands(clusters.size());

    for(size_t ci = 0; ci < clusters.size(); ci++)
    {
      const auto &c = clusters[ci];
      const int row = rowOf(c.R, rowResult.bound);
      clusterRow[ci] = row;

      int peak = nearestPeak(foldPetal(c.phi, lo, petal), peaksPerRow[row], petal);
      if(isLamRow[row][peak])
      {
        clusterIsLam[ci] = true;
        clusterIphi[ci] = -1;
        continue; // lamination cluster -- never attempt a truth match
      }

      const int iphi = iphiAssignments[row][peak];
      clusterIphi[ci] = iphi;
      if(iphi < 0) continue;

      std::vector<PetalCand> cands;
      cands.reserve(18);
      for(int p=0; p<18; p++)
      {
        unsigned int currentIndex = encodeTruthIndex(side, p, best.truthOf[row], iphi);
        double candTruthPhi = cdbttree.GetDoubleValue(currentIndex, "truthPhi");
        double dPhi = std::abs(circDiff(c.phi, candTruthPhi, 2.0*TMath::Pi()));
        cands.push_back({dPhi, p});
      }
      std::sort(cands.begin(), cands.end(), [](const PetalCand &a, const PetalCand &b){ return a.dPhi < b.dPhi; });
      clusterPetalCands[ci] = std::move(cands);
    }

    std::map<std::pair<int,int>, std::vector<int>> groups;
    for(size_t ci = 0; ci < clusters.size(); ci++)
    {
      if(clusterIphi[ci] < 0) continue;
      groups[{clusterRow[ci], clusterIphi[ci]}].push_back((int) ci);
    }

    std::vector<int> resolvedPetal(clusters.size(), -1);
    int nNoPetal = 0, nLostToConflict = 0;
    for(auto &[key, members] : groups)
    {
      std::sort(members.begin(), members.end(), [&](int a, int b)
      {
        return clusters[a].adc > clusters[b].adc; // highest ADC claims first
      });

      std::set<int> takenPetals;
      for(int ci : members)
      {
        bool gotOne = false;
        for(const auto &c : clusterPetalCands[ci])
        {
          if(c.dPhi > petal / 2.0) break;
          if(takenPetals.count(c.petal)) continue;
          resolvedPetal[ci] = c.petal;
          takenPetals.insert(c.petal);
          gotOne = true;
          break;
        }
        if(!gotOne)
        {
          bool haveValidCandidate = !clusterPetalCands[ci].empty() && clusterPetalCands[ci].front().dPhi <= petal / 2.0;
          if(haveValidCandidate) nLostToConflict++;
          else nNoPetal++;
        }
      }
    }

    for(int ci = 0; ci<(int)clusters.size(); ci++)
    {
      if (clusterIsLam[ci])
      {
        keyMatch[clusters[ci].key] = {-999, true};
      }
      else if (clusterIphi[ci] < 0)
      {
        keyMatch[clusters[ci].key] = {-999, false};
      }
      else if (resolvedPetal[ci] >= 0)
      {
        unsigned int matchedIndex = encodeTruthIndex(side, resolvedPetal[ci], best.truthOf[clusterRow[ci]], clusterIphi[ci]);
        keyMatch[clusters[ci].key] = {matchedIndex, false};
      }
    }
    return keyMatch;
  }// end matchSide


}  // namespace

LaserClusterTruthMatcher::LaserClusterTruthMatcher(const std::string &name)
  : SubsysReco(name)
{
}

int LaserClusterTruthMatcher::getNodes(PHCompositeNode *topNode)
{
  m_laserClusterContainer = findNode::getClass<LaserClusterContainer>(topNode, m_laserClusterNodeName);
  if (!m_laserClusterContainer)
  {
    std::cout << PHWHERE << " no LaserClusterContainer named " << m_laserClusterNodeName
              << " on the node tree" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  m_laserClusterHelper.set_useZ(false);
  m_laserClusterHelper.set_useDouble(true);
  m_laserClusterHelper.set_useGlobal(m_useGlobal);
  if(m_usePHGarfieldDistortions)
  {
    m_laserClusterHelper.set_garfield_cmvoltage(m_garfield_cmvoltage);
    m_laserClusterHelper.set_garfield_zerofield(m_garfield_zerofield);
    m_laserClusterHelper.set_garfield_keffside0(m_garfield_keffside0);
    m_laserClusterHelper.set_garfield_keffside1(m_garfield_keffside1);
    m_laserClusterHelper.set_garfield_stepns(m_garfield_stepns);
  }
  m_laserClusterHelper.loadNodes(topNode);

  return Fun4AllReturnCodes::EVENT_OK;
}

int LaserClusterTruthMatcher::InitRun(PHCompositeNode *topNode)
{

  // truth pattern is reloaded every run since it's cheap; distortions are
  // what actually differ run to run, and those come from m_dcc above
  delete m_cdbttree;
  m_cdbttree = new CDBTTree(m_truthFile);
  m_cdbttree->LoadCalibrations();

  for(int side=0; side<2; side++)
  {
    m_truthRowPatterns[side] = buildTruthRowPattern(*m_cdbttree, side);
    for(int row=0; row<(int)m_truthRowPatterns[side].size(); row++)
    {
      TruthRow newRow;
      unsigned int truthIndex = (side ? 18 : 0)*10000 + (row*100) + 0;
      newRow.R = m_cdbttree->GetDoubleValue(truthIndex, "truthR");
      newRow.dphi = m_truthRowPatterns[side][row].stripePhi[1] - m_truthRowPatterns[side][row].stripePhi[0];
      newRow.nstripes = m_truthRowPatterns[side][row].stripePhi.size();

      m_truthRows[side].push_back(newRow);
    }
  }

  return getNodes(topNode);
}

int LaserClusterTruthMatcher::process_event(PHCompositeNode * /*topNode*/)
{
  // the aggregated clusters represent one run's worth of laser data
  // collapsed into a single pseudo-event, so this whole-side pipeline only
  // needs to run once per side per call here -- there is no meaningful
  // "next event" for this data

  std::vector<Cluster> clusters[2];
  
  auto clusrange = m_laserClusterContainer->getClusters();
  for (auto cmitr = clusrange.first; cmitr != clusrange.second; ++cmitr)
  {
    const auto &[cmkey, cmclus_orig] = *cmitr;
    LaserCluster *cmclus = cmclus_orig;

    int side = TpcDefs::getSide(cmkey);

    Acts::Vector3 pos;
    if(m_usePHGarfieldDistortions) pos = m_laserClusterHelper.getClusterCentroidWithPHGarfield(cmclus);
    else pos = m_laserClusterHelper.getClusterCentroid(cmclus);
    if(pos.hasNaN())
    {
      continue;
    }

    double R = sqrt(pos[0]*pos[0] + pos[1]*pos[1]);
    double phi = atan2(pos[1], pos[0]);

    clusters[side].push_back({R, phi, cmclus->getAdcDouble(), cmkey});
  }

  for(int side=0; side<2; side++)
  {
    std::map<TrkrDefs::cluskey, std::pair<int,bool>> matched = matchSide(clusters[side], *m_cdbttree, m_truthRowPatterns[side], m_truthRows[side], side, m_QABase, Verbosity());
    for(auto &[key, matchedIndex] : matched)
    {
      LaserCluster *clus = m_laserClusterContainer->findCluster(key);
      if(clus)
      {
        clus->setTruthIndex(matchedIndex.first);
        clus->setIsLamination(matchedIndex.second);
      }
    }
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int LaserClusterTruthMatcher::End(PHCompositeNode * /*topNode*/)
{
  delete m_cdbttree;
  m_cdbttree = nullptr;
  return Fun4AllReturnCodes::EVENT_OK;
}