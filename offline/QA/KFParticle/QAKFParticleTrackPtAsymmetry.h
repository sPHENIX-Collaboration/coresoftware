#ifndef QAKFPARTICLETRACKPTASYMMETRY_H
#define QAKFPARTICLETRACKPTASYMMETRY_H

#include <string>
#include <vector>

class KFParticle;
class SvtxTrackMap;
class Fun4AllHistoManager;
class TH1;
class TH2;

class QAKFParticleTrackPtAsymmetry
{
  public:
    QAKFParticleTrackPtAsymmetry(const std::string &histo_prefix,                                                                      //
                                 double min_m = 0.4,                                                                                   //
                                 double max_m = 0.6,                                                                                   //
                                 const std::vector<double> &mother_eta_bins = std::vector<double>{-2.0, -1.0, -0.5, 0, 0.5, 1.0, 2.0}, //
                                 const std::vector<double> &mother_phi_bins = std::vector<double>{-3.2, -2.5, -2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.2});

    void bookHistograms(Fun4AllHistoManager *hm);

    void analyzeTrackPtAsymmetry(SvtxTrackMap *trackMap, KFParticle *mother);

    void setVerbosity(int verb) { m_verbosity = verb; }

  protected:
  private:
    int m_verbosity{0};
    std::string m_prefix;
    double m_min_mass{0.4}; // default to Kshort mass window
    double m_max_mass{0.6}; // default to Kshort mass window
    // mother eta and phi bins
    std::vector<double> m_mother_eta_bins{-2.0, -1.0, -0.5, 0, 0.5, 1.0, 2.0};                                    // default bins
    std::vector<double> m_mother_phi_bins{-3.2, -2.5, -2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.2}; // default bins

    // basic histograms
    TH1 *h_trackPtAsymmetry{nullptr};
    TH2 *h2_trackPtAsymmetry_vs_mass{nullptr};
    // differential histograms in (eta, phi) bins
    std::vector<std::vector<TH2 *>> h2_trackPtAsymmetry_eta_phi_bins;
};

#endif // QAKFPARTICLETRACKPTASYMMETRY_H