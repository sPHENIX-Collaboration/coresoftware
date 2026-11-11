#ifndef QA_TPCCALIBRATIONQA_INTERNAL_H
#define QA_TPCCALIBRATIONQA_INTERNAL_H

#include <TCanvas.h>
#include <TDirectory.h>
#include <TF1.h>
#include <TFile.h>
#include <TFitResult.h>
#include <TFitResultPtr.h>
#include <TGraphErrors.h>
#include <TH1D.h>
#include <TH2.h>
#include <TH2D.h>
#include <TH3.h>
#include <TKey.h>
#include <TLegend.h>
#include <TLine.h>
#include <TNamed.h>
#include <TPad.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TSystem.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace TpcCalibrationQADetail
{
  using Hist = std::unique_ptr<TH1>;
  using Hists = std::map<std::string, Hist>;
  // ROOT installations may still use C++17; Jenkins checks also support C++20.
  template <class Container, class Key>
  bool hasKey(const Container &container, const Key &key)
  {
#if __cplusplus >= 202002L
    return container.contains(key);
#else
    return container.find(key) != container.end();
#endif
  }
  struct Sample
  {
    std::string label;
    std::string list;
    std::string driftList;
    double rawTimeShift{std::numeric_limits<double>::quiet_NaN()};
  };
  struct Curve
  {
    std::string title;
    std::string xTitle;
    std::string yTitle;
    std::string comparison{"difference"};
    std::vector<double> x;
    std::vector<double> ex;
    std::vector<double> y;
    std::vector<double> ey;
    std::vector<double> edges;
  };
  using Curves = std::map<std::string, Curve>;
  Curve makeCurve(std::string title, std::string xTitle, std::string yTitle,
                  std::string comparison = "difference");
  struct Config
  {
    double ptMin, ptMax;
    int minEntries, phiGroup, etaGroup;
    std::string prefix;
    bool saveSlices;
    double ratioFloor;
    bool phiQpt;
    double rawTimeShift;
    bool fitBackground{false};
  };
  struct Fit
  {
    int status{1}, rootStatus{-1}, covariance{-1};
    double background{0};
    int ndf{0};
    double n{0}, amplitude{0}, mean{0}, meanError{0}, sigma{0}, sigmaError{0}, low{0}, high{0},
        chi2ndf{0};
  };

  struct Location
  {
    std::string folder, name;
  };
  extern std::map<std::string, std::string> outputIndex;
  struct DriftInput
  {
    Hist time[2], rawTime[2], z[2];
    double appliedShift[2]{0., 0.};
    Hist conditions;
    double referenceV[2]{NAN, NAN}, referenceT0[2]{NAN, NAN};
    std::string referenceSource[2];
    double phgMean[2]{0, 0}, weightedV[2]{0, 0}, weightedADC[2]{0, 0};
    double calibrations[2]{0, 0};
    double pad{102.605}, cmHalf{.28}, clock{56.8};
    std::string source[2];
  };
  struct Edge
  {
    bool valid{false};
    double x10{0}, x50{0}, x90{0}, left{0}, right{0};
  };
  bool sameAxis(const TAxis *a, const TAxis *b);
  bool sameBins(const TH1 *a, const TH1 *b);
  std::vector<std::string> wanted();
  Hists load(const Sample &sample, const Config &cfg, std::ostream &availability);
  void point(Curve &c, double x, double ex, double y, double ey);
  std::unique_ptr<TGraphErrors> graph(const Curve &c, const std::string &name);
  std::vector<double> edges(const TAxis *axis);
  TDirectory *folder(TDirectory *root, const std::string &path);
  bool ends(const std::string &s, const std::string &suffix);
  std::string plain(std::string s);
  Location location(std::string key, const std::string &comparison = "difference");
  Location sliceLocation(const std::string &key);
  std::string writeCurve(const Curve &c, const std::string &name, TDirectory *dir, bool flat = false);
  void slices(TH2 *h, const std::string &key, const std::string &title, Curves &curves,
              const Sample &sample, const Config &cfg, TDirectory *dir, std::ostream &csv);
  bool samePoint(const Curve &a, size_t i, const Curve &b, size_t j);
  Curve combine(const Curve &a, const Curve &b, const std::string &op, double ratioFloor = .02);
  void chargeCurves(TH1 *pos, TH1 *neg, const std::string &key, const std::string &title,
                    Curves &curves, double frames);
  std::pair<int, int> ptRange(const TAxis *axis, double low, double high);
  void rebin(TH2 *h, int group);
  void rebin1(TH1 *h, int group);
  void phiQoverPt(TH3 *pos, TH3 *neg, const std::string &side, const Sample &sample,
                  const Config &cfg, TDirectory *dir, Curves &curves, std::ostream &csv);
  void plot(const Curves &curves, const Curves &reference, const Sample &sample,
            const std::string &refLabel, const std::string &pdf, TDirectory *dir, const Config &cfg);
  void plotDcaDistributions(TDirectory *dir, const std::string &pdf);
  double median(std::vector<double> x);
  double windowMedian(const std::vector<double> &x, const std::vector<double> &y, double lo, double hi);
  Edge edge(const std::vector<double> &x, const std::vector<double> &raw, double ll, double lh, double rl,
            double rh, double sl, double sh, double expected, int direction, bool clean = true);
  Hist clone1(TH1 *h, const std::string &name);
  Hist readTime(TFile *file, int side, std::string &source);
  Hist readZ(TFile *file, int side);
  std::pair<std::vector<double>, std::vector<double>> xy(const TH1 *h, int side, bool fold);
  std::pair<std::vector<double>, std::vector<double>> physicalZ(const TH1 *h, int side, double pad);
}  // namespace TpcCalibrationQADetail
#endif
