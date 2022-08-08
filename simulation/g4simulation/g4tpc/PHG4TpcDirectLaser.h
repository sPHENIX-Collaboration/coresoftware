#ifndef G4TPC_PHG4TPCDIRECTLASER_H
#define G4TPC_PHG4TPCDIRECTLASER_H

#include <fun4all/SubsysReco.h>

#include <phparameter/PHParameterInterface.h>

#include <TVector3.h>

#include <cmath>
#include <string>  // for string, allocator
#include <vector>  // for vector

class PHG4HitContainer;
class SvtxTrackMap;
class PHG4TruthInfoContainer;
class PHCompositeNode;

class PHG4TpcDirectLaser : public SubsysReco, public PHParameterInterface
{
 public:
  /// constructor
  PHG4TpcDirectLaser(const std::string &name = "PHG4TpcDirectLaser");

  /// destructor
  ~PHG4TpcDirectLaser() override = default;

  /// run initialization
  int InitRun(PHCompositeNode *) override;

  /// per event processing
  int process_event(PHCompositeNode *) override;

  /// default parameters
  void SetDefaultParameters() override;

  /// detector name
  void Detector(const std::string &d)
  {
    detector = d;
  }

  /// define steps along phi
  void SetPhiStepping(int n, double min, double max);

  /// define steps along theta
  void SetThetaStepping(int n, double min, double max);

  /// get total number of steps
  int GetNpatternSteps() const
  {
    return nPhiSteps * nThetaSteps;
  };

  /// set current patter step
  void SetCurrentPatternStep(int value)
  {
    currentPatternStep = value;
  }

  /// advance automatically through patterns
  void SetDirectLaserAuto(bool value)
  {
    m_autoAdvanceDirectLaser = value;
  };

  void SetArbitraryThetaPhi(double theta, double phi)
  {
    arbitrary_theta = theta;
    arbitrary_phi = phi;
  }

 private:
  /// define lasers
  /* by default there are 4 lasers on each side of the TPC */
  void SetupLasers();

  /// aim lasers to a given theta and phi angle
  void AimToThetaPhi(double theta, double phi);

  /// aim lasers to a give step
  void AimToPatternStep(int n);

  /// aim to next step
  void AimToNextPatternStep();

  /// stores laser position and direction along z
  class Laser
  {
   public:
    /// laser position
    TVector3 m_position;

    /// laser phi position
    double m_phi = 0;

    /// laser direction along z
    int m_direction = 1;
  };

  /// append track in given angular direction and for a given laser
  void AppendLaserTrack(double theta, double phi, const Laser &);

  /// detector name
  std::string detector = "TPC";

  /// g4hitnode name
  std::string hitnodename;

  /// lasers
  std::vector<Laser> m_lasers;

  /// number of electrons deposited per cm laser track
  int electrons_per_cm = 300;

  // number of electrons per deposited GeV in TPC gas
  /** 
   * it is used to convert a given number of electrons into an energy 
   * as expected by G4Hit. The energy is then converted back to a number of electrons
   * inside PHG4TpcElectronDrift
   */
  double electrons_per_gev = NAN;

  double arbitrary_theta = -30.0;  // degrees
  double arbitrary_phi = -30.0;    // degrees

  ///@name default phi and theta steps
  //@{
  int nPhiSteps = 1;
  int nThetaSteps = 1;
  int nTotalSteps = 1;
  double minPhi = 0;
  double maxPhi = 0;
  double minTheta = 0;
  double maxTheta = 0;
  //@}

  // current patter step
  int currentPatternStep = 0;

  /// set to true to change direct laser tracks from one event to the other
  bool m_autoAdvanceDirectLaser = false;

  /// g4hit container
  PHG4HitContainer *m_g4hitcontainer = nullptr;

  //! truth information
  PHG4TruthInfoContainer *m_g4truthinfo = nullptr;

  /// track map, used to store track parameters
  std::string m_track_map_name = "SvtxTrackMap";
  SvtxTrackMap *m_track_map = nullptr;
};

#endif
