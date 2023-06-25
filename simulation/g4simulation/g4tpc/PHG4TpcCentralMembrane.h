#ifndef PHG4TPCCENTRALMEMBRANE_H
#define PHG4TPCCENTRALMEMBRANE_H

#include <fun4all/SubsysReco.h>
#include <phparameter/PHParameterInterface.h>

#include <array>
#include <cmath>
#include <string>  // for string
#include <vector>

class PHCompositeNode;
class PHG4Hit;

// all distances in mm, all angles in rad
// class that generates stripes and dummy hit coordinates
// stripes have width of one mm, length of one pad width, and are centered in middle of sector gaps

class PHG4TpcCentralMembrane : public SubsysReco, public PHParameterInterface
{
 public:
  /// constructor
  PHG4TpcCentralMembrane(const std::string& name = "PHG4TpcCentralMembrane");

  /// destructor
  ~PHG4TpcCentralMembrane() override;

  /// run initialization
  int InitRun(PHCompositeNode*) override;

  /// per event processing
  int process_event(PHCompositeNode*) override;

  /// default parameters
  void SetDefaultParameters() override;

  /// detector name
  void Detector(const std::string& d)
  {
    detector = d;
    hitnodename = "G4HIT_" + d;
  }

  /// check if coords are in a stripe
  int getSearchResult(double xcheck, double ycheck) const;

  int getStripeID(double xcheck, double ycheck) const;

  /// adjust central membrane hits delay with respect to trigger time
  void setCentralMembraneDelay(int ns) { m_centralMembraneDelay = ns; };

  /// set modulo for events in which to generate CM hits
  void setCentralMembraneEventModulo(int mod) { m_eventModulo = mod; };

 private:
  /// detector name
  std::string detector = "TPC";

  /// g4hitnode name
  std::string hitnodename = "G4HIT_TPC";
  std::vector<PHG4Hit*> PHG4Hits;

  int m_eventModulo = 10;
  int m_eventNum = 0;

  static constexpr double mm = 1.0;
  static constexpr double cm = 10.0;

  /// inner radius of CM
  static constexpr double begin_CM = 221.4019814 * mm;

  /// outer radius of CM
  static constexpr double end_CM = 759.2138 * mm;

  static constexpr int nRadii = 8;
  static constexpr int nStripes_R1 = 6;
  static constexpr int nStripes_R2 = 8;
  static constexpr int nStripes_R3 = 12;

  static constexpr int nPads_R1 = 6 * 16;
  static constexpr int nPads_R2 = 8 * 16;
  static constexpr int nPads_R3 = 12 * 16;

  //
  static constexpr double padfrac_R1 = 0.5 * 5.59106385 * mm;
  static constexpr double padfrac_R2 = 0.5 * 10.13836283 * mm;
  static constexpr double padfrac_R3 = 0.5 * 10.90189537 * mm;

  /// radius of arc on end of a stripe
  static constexpr double arc_r = 0.5 * mm;

  /// stripe radii
  static constexpr std::array<double, nRadii> R1_e = {{227.0902789 * mm, 238.4100043 * mm, 249.7297296 * mm, 261.049455 * mm, 272.3691804 * mm, 283.6889058 * mm, 295.0086312 * mm, 306.3283566 * mm}};
  static constexpr std::array<double, nRadii> R1 = {{317.648082 * mm, 328.9678074 * mm, 340.2875328 * mm, 351.6072582 * mm, 362.9269836 * mm, 374.246709 * mm, 385.5664344 * mm, 396.8861597 * mm}};
  static constexpr std::array<double, nRadii> R2 = {{421.705532 * mm, 442.119258 * mm, 462.532984 * mm, 482.9467608 * mm, 503.36069 * mm, 523.774416 * mm, 544.188015 * mm, 564.601868 * mm}};
  static constexpr std::array<double, nRadii> R3 = {{594.6048725 * mm, 616.545823 * mm, 638.4867738 * mm, 660.4277246 * mm, 682.3686754 * mm, 704.3096262 * mm, 726.250577 * mm, 748.1915277 * mm}};

  double str_width_R1_e[nStripes_R1][nRadii] = {};
  double str_width_R1[nStripes_R1][nRadii] = {};
  double str_width_R2[nStripes_R2][nRadii] = {};
  double str_width_R3[nStripes_R3][nRadii] = {};

  static constexpr std::array<double, nRadii> widthmod_R1_e = {{1.493, 1.398, 1.334, 1.284, 1.243, 1.208, 1.178, 1.152}};
  static constexpr std::array<double, nRadii> widthmod_R1 = {{1.129, 1.109, 1.091, 1.076, 1.062, 1.050, 1.040, 1.030}};
  static constexpr std::array<double, nRadii> widthmod_R2 = {{1.015, 1.007, 1.002, 1.000, 1.001, 1.006, 1.013, 1.023}};
  static constexpr std::array<double, nRadii> widthmod_R3 = {{1.044, 1.064, 1.087, 1.115, 1.147, 1.186, 1.232, 1.288}};

  std::array<double, nRadii> spacing_R1_e = {};
  std::array<double, nRadii> spacing_R1 = {};
  std::array<double, nRadii> spacing_R2 = {};
  std::array<double, nRadii> spacing_R3 = {};

  //bottom left - 1a
  double x1a_R1_e[nStripes_R1][nRadii], y1a_R1_e[nStripes_R1][nRadii];
  double x1a_R1[nStripes_R1][nRadii], y1a_R1[nStripes_R1][nRadii];
  double x1a_R2[nStripes_R2][nRadii], y1a_R2[nStripes_R2][nRadii];
  double x1a_R3[nStripes_R3][nRadii], y1a_R3[nStripes_R3][nRadii];

  //bottom right - 1b
  double x1b_R1_e[nStripes_R1][nRadii], y1b_R1_e[nStripes_R1][nRadii];
  double x1b_R1[nStripes_R1][nRadii], y1b_R1[nStripes_R1][nRadii];
  double x1b_R2[nStripes_R2][nRadii], y1b_R2[nStripes_R2][nRadii];
  double x1b_R3[nStripes_R3][nRadii], y1b_R3[nStripes_R3][nRadii];

  //top left - 2a
  double x2a_R1_e[nStripes_R1][nRadii], y2a_R1_e[nStripes_R1][nRadii];
  double x2a_R1[nStripes_R1][nRadii], y2a_R1[nStripes_R1][nRadii];
  double x2a_R2[nStripes_R2][nRadii], y2a_R2[nStripes_R2][nRadii];
  double x2a_R3[nStripes_R3][nRadii], y2a_R3[nStripes_R3][nRadii];

  //top right - 2b
  double x2b_R1_e[nStripes_R1][nRadii], y2b_R1_e[nStripes_R1][nRadii];
  double x2b_R1[nStripes_R1][nRadii], y2b_R1[nStripes_R1][nRadii];
  double x2b_R2[nStripes_R2][nRadii], y2b_R2[nStripes_R2][nRadii];
  double x2b_R3[nStripes_R3][nRadii], y2b_R3[nStripes_R3][nRadii];

  //left midpoint - 3a
  double x3a_R1_e[nStripes_R1][nRadii], y3a_R1_e[nStripes_R1][nRadii];
  double x3a_R1[nStripes_R1][nRadii], y3a_R1[nStripes_R1][nRadii];
  double x3a_R2[nStripes_R2][nRadii], y3a_R2[nStripes_R2][nRadii];
  double x3a_R3[nStripes_R3][nRadii], y3a_R3[nStripes_R3][nRadii];

  //right midpoint - 3b
  double x3b_R1_e[nStripes_R1][nRadii], y3b_R1_e[nStripes_R1][nRadii];
  double x3b_R1[nStripes_R1][nRadii], y3b_R1[nStripes_R1][nRadii];
  double x3b_R2[nStripes_R2][nRadii], y3b_R2[nStripes_R2][nRadii];
  double x3b_R3[nStripes_R3][nRadii], y3b_R3[nStripes_R3][nRadii];

  //Check which stripes get removed
  std::array<int, nRadii> nGoodStripes_R1_e = {};
  std::array<int, nRadii> nGoodStripes_R1 = {};
  std::array<int, nRadii> nGoodStripes_R2 = {};
  std::array<int, nRadii> nGoodStripes_R3 = {};

  ///min stripe index
  static constexpr std::array<int, nRadii> keepThisAndAfter = {{1, 0, 1, 0, 1, 0, 1, 0}};

  ///max stripe index
  static constexpr std::array<int, nRadii> keepUntil_R1_e = {{4, 4, 5, 4, 5, 5, 5, 5}};
  static constexpr std::array<int, nRadii> keepUntil_R1 = {{5, 5, 6, 5, 6, 5, 6, 5}};
  static constexpr std::array<int, nRadii> keepUntil_R2 = {{7, 7, 8, 7, 8, 8, 8, 8}};
  static constexpr std::array<int, nRadii> keepUntil_R3 = {{11, 10, 11, 11, 11, 11, 12, 11}};

  std::array<int, nRadii> nStripesIn_R1_e = {};
  std::array<int, nRadii> nStripesIn_R1 = {};
  std::array<int, nRadii> nStripesIn_R2 = {};
  std::array<int, nRadii> nStripesIn_R3 = {};
  std::array<int, nRadii> nStripesBefore_R1_e = {};
  std::array<int, nRadii> nStripesBefore_R1 = {};
  std::array<int, nRadii> nStripesBefore_R2 = {};
  std::array<int, nRadii> nStripesBefore_R3 = {};

  static constexpr int nStripesPerPetal = 213;
  static constexpr int nPetals = 18;
  static constexpr int nTotStripes = nStripesPerPetal * nPetals;

  /// mean number of electrons per stripe
  int electrons_per_stripe = 300;

  // number of electrons per deposited GeV in TPC gas
  /** 
   * it is used to convert a given number of electrons into an energy 
   * as expected by G4Hit. The energy is then converted back to a number of electrons
   * inside PHG4TpcElectronDrift
   */
  double electrons_per_gev = NAN;

  /// delay between central membrane hits and trigger time (ns)
  int m_centralMembraneDelay = 0;

  void CalculateVertices(
      int nStripes, int nPads,
      const std::array<double, nRadii>& R,
      std::array<double, nRadii>& spacing,
      double x1a[][nRadii], double y1a[][nRadii],
      double x1b[][nRadii], double y1b[][nRadii],
      double x2a[][nRadii], double y2a[][nRadii],
      double x2b[][nRadii], double y2b[][nRadii],
      double x3a[][nRadii], double y3a[][nRadii],
      double x3b[][nRadii], double y3b[][nRadii],
      double padfrac,
      double str_width[][nRadii],
      const std::array<double, nRadii>& widthmod,
      std::array<int, nRadii>& nGoodStripes,
      const std::array<int, nRadii>& keepUntil,
      std::array<int, nRadii>& nStripesIn,
      std::array<int, nRadii>& nStripesBefore);

  int SearchModule(int nStripes,
                   const double x1a[][nRadii], const double x1b[][nRadii],
                   const double x2a[][nRadii], const double x2b[][nRadii],
                   const double y1a[][nRadii], const double y1b[][nRadii],
                   const double y2a[][nRadii], const double y2b[][nRadii],
                   const double x3a[][nRadii], const double y3a[][nRadii],
                   const double x3b[][nRadii], const double y3b[][nRadii],
                   double x, double y, const std::array<int, nRadii>& nGoodStripes) const;

  PHG4Hit* GetPHG4HitFromStripe(int petalID, int moduleID, int radiusID, int stripeID, int nElectrons) const;
};

#endif
