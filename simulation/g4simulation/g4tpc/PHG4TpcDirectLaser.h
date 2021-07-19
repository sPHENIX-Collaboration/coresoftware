#ifndef G4TPC_PHG4TPCDIRECTLASER_H
#define G4TPC_PHG4TPCDIRECTLASER_H

#include <fun4all/SubsysReco.h>

#include <TVector3.h>

class PHG4HitContainer;

class PHG4TpcDirectLaser: public SubsysReco 
{
  public:
  
  /// constructor
  PHG4TpcDirectLaser( const std::string& name = "PHG4TpcDirectLaser" ); 
  ~PHG4TpcDirectLaser() override = default;

  /// run initialization
  int InitRun(PHCompositeNode *) override;

  /// per event processing
  int process_event(PHCompositeNode *) override;

  /// detector name
  void Detector(const std::string &d)
  { detector = d; }
  
  /// define steps along phi
  void SetPhiStepping(int n, float min,float max);
  
  /// define steps along theta
  void SetThetaStepping(int n, float min,float max);
  
  /// get total number of steps
  int GetNpatternSteps() const 
  {return nPhiSteps*nThetaSteps;};

  /// set current patter step
  void SetCurrentPatternStep( int value )
  { currentPatternStep = value; }
  
  /// advance automatically through patterns
  void SetDirectLaserAuto(bool value)
  {m_autoAdvanceDirectLaser = value;};
  
  private:
  
  /// define lasers
  /* by default there are 4 lasers on each side of the TPC */
  void SetupLasers();
  
  /// aim lasers to a given theta and phi angle
  void AimToThetaPhi(float theta, float phi);
  
  /// aim lasers to a give step
  void AimToPatternStep(int n);

  /// aim to next step
  void AimToNextPatternStep();
  
  TVector3 GetCmStrike(TVector3 start, TVector3 direction) const;

  TVector3 GetFieldcageStrike(TVector3 start, TVector3 direction) const;
  
  TVector3 GetCylinderStrike(TVector3 s, TVector3 v, float radius) const;

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
  void AppendLaserTrack(float theta, float phi, const Laser&);
    
  ///@name units
  //@{
  static constexpr double mm = 0.10;
  static constexpr double cm = 1.0;
  //@}
    
  /// length of generated G4Hits along laser track
  static constexpr float maxHitLength=1.0;//1cm.

  // inner and outer radii of field cages/TPC
  static constexpr double begin_CM = 20.*cm ;
  static constexpr double end_CM = 78.*cm;
  
  //half the thickness of the CM;
  static constexpr double halfwidth_CM = 0.5*cm;
  
  /// detector name
  std::string detector = "TPC";
  
  /// g4hitnode name
  std::string hitnodename;

  /// lasers
  std::vector<Laser> m_lasers;
  
  /// g4hit container
  PHG4HitContainer* m_g4hitcontainer = nullptr;
  
  ///@name default phi and theta steps
  //@{
  int nPhiSteps=1;
  int nThetaSteps=1;
  int nTotalSteps=1;
  float minPhi=0;
  float maxPhi=0;
  float minTheta=0;
  float maxTheta=0;
  //@}

  // current patter step
  int currentPatternStep=0;

  /// set to true to change direct laser tracks from one event to the other
  bool m_autoAdvanceDirectLaser=false;
  
};


#endif
