#include "PHG4TpcDirectLaser.h"

#include <phparameter/PHParameterInterface.h>  // for PHParameterInterface

#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4HitDefs.h>  // for get_volume_id
#include <g4main/PHG4Hitv1.h>
#include <g4main/PHG4Particlev3.h>
#include <g4main/PHG4TruthInfoContainer.h>
#include <g4main/PHG4VtxPointv1.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>  // for SubsysReco

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>  // for PHObject
#include <phool/getClass.h>
#include <phool/phool.h>  // for PHWHERE

#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/SvtxTrackMap_v2.h>
#include <trackbase_historic/SvtxTrack_v2.h>

#include <TVector3.h>  // for TVector3, operator*

#include <gsl/gsl_const_mksa.h>  // for the speed of light

#include <cassert>
#include <iostream>  // for operator<<, basic_os...
#include <optional>

namespace
{
  using PHG4Particle_t = PHG4Particlev3;
  using PHG4VtxPoint_t = PHG4VtxPointv1;
  using PHG4Hit_t = PHG4Hitv1;

  // utility
  template <class T>
  inline constexpr T square(const T& x)
  {
    return x * x;
  }

  // unique detector id for all direct lasers
  static const int detId = PHG4HitDefs::get_volume_id("PHG4TpcDirectLaser");

  ///@name units
  //@{
  static constexpr double cm = 1.0;
  //@}

  /// speed of light, in cm per ns
  static constexpr double speed_of_light = GSL_CONST_MKSA_SPEED_OF_LIGHT * 1e-7;

  /// length of generated G4Hits along laser track
  static constexpr double maxHitLength = 1. * cm;

  /// TPC half length
  static constexpr double halflength_tpc = 105.5 * cm;

  // inner and outer radii of field cages/TPC
  static constexpr double begin_CM = 20. * cm;
  static constexpr double end_CM = 78. * cm;

  //half the thickness of the CM;
  static constexpr double halfwidth_CM = 0.5 * cm;

  //_____________________________________________________________
  std::optional<TVector3> central_membrane_intersection(TVector3 start, TVector3 direction)
  {
    const double end = start.z() > 0 ? halfwidth_CM : -halfwidth_CM;
    const double dist = end - start.z();

    // if line is vertical, it will never intercept the endcap
    if (direction.z() == 0) return std::nullopt;

    // check that distance and direction have the same sign
    if (dist * direction.z() < 0) return std::nullopt;

    const double direction_scale = dist / direction.z();
    return start + direction * direction_scale;
  }

  //_____________________________________________________________
  std::optional<TVector3> endcap_intersection(TVector3 start, TVector3 direction)
  {
    const double end = start.z() > 0 ? halflength_tpc : -halflength_tpc;
    const double dist = end - start.z();

    // if line is vertical, it will never intercept the endcap
    if (direction.z() == 0) return std::nullopt;

    // check that distance and direction have the same sign
    if (dist * direction.z() < 0) return std::nullopt;

    const double direction_scale = dist / direction.z();
    return start + direction * direction_scale;
  }

  //_____________________________________________________________
  std::optional<TVector3> cylinder_line_intersection(TVector3 s, TVector3 v, double radius)
  {
    const double R2 = square(radius);

    //Generalized Parameters for collision with cylinder of radius R:
    //from quadratic formula solutions of when a vector intersects a circle:
    const double a = square(v.x()) + square(v.y());
    const double b = 2 * (v.x() * s.x() + v.y() * s.y());
    const double c = square(s.x()) + square(s.y()) - R2;

    const double rootterm = square(b) - 4 * a * c;

    /*
     * if a==0 then we are parallel and will have no solutions.
     * if the rootterm is negative, we will have no real roots,
     * we are outside the cylinder and pointing skew to the cylinder such that we never cross.
     */
    if (rootterm < 0 || a == 0) return std::nullopt;

    //Find the (up to) two points where we collide with the cylinder:
    const double sqrtterm = std::sqrt(rootterm);
    const double t1 = (-b + sqrtterm) / (2 * a);
    const double t2 = (-b - sqrtterm) / (2 * a);

    /*
    * if either of the t's are nonzero, we have a collision
    * the collision closest to the start (hence with the smallest t that is greater than zero) is the one that happens.
    */
    const double& min_t = (t2 < t1 && t2 > 0) ? t2 : t1;
    return s + v * min_t;
  }

  //_____________________________________________________________
  std::optional<TVector3> field_cage_intersection(TVector3 start, TVector3 direction)
  {
    const auto ofc_strike = cylinder_line_intersection(start, direction, end_CM);
    const auto ifc_strike = cylinder_line_intersection(start, direction, begin_CM);

    // if either of the two intersection is invalid, return the other
    if (!ifc_strike) return ofc_strike;
    if (!ofc_strike) return ifc_strike;

    // both intersection are valid, calculate signed distance to start z
    const auto ifc_dist = (ifc_strike->Z() - start.Z()) / direction.Z();
    const auto ofc_dist = (ofc_strike->Z() - start.Z()) / direction.Z();

    if (ifc_dist < 0)
      return (ofc_dist > 0) ? ofc_strike : std::nullopt;
    else if (ofc_dist < 0)
      return ifc_strike;
    else
      return (ifc_dist < ofc_dist) ? ifc_strike : ofc_strike;
  }

  /// TVector3 stream
  inline std::ostream& operator<<(std::ostream& out, const TVector3& vector)
  {
    out << "( " << vector.x() << ", " << vector.y() << ", " << vector.z() << ")";
    return out;
  }

}  // namespace

//_____________________________________________________________
PHG4TpcDirectLaser::PHG4TpcDirectLaser(const std::string& name)
  : SubsysReco(name)
  , PHParameterInterface(name)
{
  InitializeParameters();
}

//_____________________________________________________________
int PHG4TpcDirectLaser::InitRun(PHCompositeNode* topNode)
{
  // g4 truth info
  m_g4truthinfo = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");
  if (!m_g4truthinfo)
  {
    std::cout << "Fun4AllDstPileupMerger::load_nodes - creating node G4TruthInfo" << std::endl;

    PHNodeIterator iter(topNode);
    auto dstNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST"));
    if (!dstNode)
    {
      std::cout << PHWHERE << "DST Node missing, aborting." << std::endl;
      return Fun4AllReturnCodes::ABORTRUN;
    }

    m_g4truthinfo = new PHG4TruthInfoContainer();
    dstNode->addNode(new PHIODataNode<PHObject>(m_g4truthinfo, "G4TruthInfo", "PHObject"));
  }

  // load and check G4Hit node
  hitnodename = "G4HIT_" + detector;
  auto* g4hit = findNode::getClass<PHG4HitContainer>(topNode, hitnodename.c_str());
  if (!g4hit)
  {
    std::cout << Name() << " Could not locate G4HIT node " << hitnodename << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  // find or create track map
  /* it is used to store laser parameters on a per event basis */
  m_track_map = findNode::getClass<SvtxTrackMap>(topNode, m_track_map_name);
  if (!m_track_map)
  {
    // find DST node and check
    PHNodeIterator iter(topNode);
    auto dstNode = static_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST"));
    if (!dstNode)
    {
      std::cout << PHWHERE << "DST Node missing, aborting." << std::endl;
      return Fun4AllReturnCodes::ABORTRUN;
    }

    // find or create SVTX node
    iter = PHNodeIterator(dstNode);
    auto node = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "SVTX"));
    if (!node) dstNode->addNode(node = new PHCompositeNode("SVTX"));

    // add track node
    m_track_map = new SvtxTrackMap_v2;
    node->addNode(new PHIODataNode<PHObject>(m_track_map, m_track_map_name, "PHObject"));
  }

  // setup parameters
  UpdateParametersWithMacro();
  electrons_per_cm = get_int_param("electrons_per_cm");
  electrons_per_gev = get_double_param("electrons_per_gev");

  // setup lasers
  SetupLasers();

  // print configuration
  std::cout << "PHG4TpcDirectLaser::InitRun - m_autoAdvanceDirectLaser: " << m_autoAdvanceDirectLaser << std::endl;
  std::cout << "PHG4TpcDirectLaser::InitRun - phi steps: " << nPhiSteps << " min: " << minPhi << " max: " << maxPhi << std::endl;
  std::cout << "PHG4TpcDirectLaser::InitRun - theta steps: " << nThetaSteps << " min: " << minTheta << " max: " << maxTheta << std::endl;
  std::cout << "PHG4TpcDirectLaser::InitRun - nTotalSteps: " << nTotalSteps << std::endl;

  std::cout << "PHG4TpcDirectLaser::InitRun - electrons_per_cm: " << electrons_per_cm << std::endl;
  std::cout << "PHG4TpcDirectLaser::InitRun - electrons_per_gev " << electrons_per_gev << std::endl;

  return Fun4AllReturnCodes::EVENT_OK;
}

//_____________________________________________________________
int PHG4TpcDirectLaser::process_event(PHCompositeNode* topNode)
{
  // g4 input event
  m_g4truthinfo = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");
  assert(m_g4truthinfo);

  // load g4hit container
  m_g4hitcontainer = findNode::getClass<PHG4HitContainer>(topNode, hitnodename.c_str());
  assert(m_g4hitcontainer);

  // load track map
  m_track_map = findNode::getClass<SvtxTrackMap>(topNode, m_track_map_name);
  assert(m_track_map);

  if (m_autoAdvanceDirectLaser)
  {
    AimToNextPatternStep();
  }
  else
  {
    // use arbitrary direction
    AimToThetaPhi(M_PI / 180. * arbitrary_theta, M_PI / 180. * arbitrary_phi);
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

//_____________________________________________________________
void PHG4TpcDirectLaser::SetDefaultParameters()
{
  // same gas parameters as in PHG4TpcElectronDrift::SetDefaultParameters

  // Data on gasses @20 C and 760 Torr from the following source:
  // http://www.slac.stanford.edu/pubs/icfa/summer98/paper3/paper3.pdf
  // diffusion and drift velocity for 400kV for NeCF4 50/50 from calculations:
  // http://skipper.physics.sunysb.edu/~prakhar/tpc/HTML_Gases/split.html
  static constexpr double Ne_dEdx = 1.56;    // keV/cm
  static constexpr double CF4_dEdx = 7.00;   // keV/cm
  static constexpr double Ne_NTotal = 43;    // Number/cm
  static constexpr double CF4_NTotal = 100;  // Number/cm
  static constexpr double Tpc_NTot = 0.5 * Ne_NTotal + 0.5 * CF4_NTotal;
  static constexpr double Tpc_dEdx = 0.5 * Ne_dEdx + 0.5 * CF4_dEdx;
  static constexpr double Tpc_ElectronsPerKeV = Tpc_NTot / Tpc_dEdx;

  // number of electrons per deposited GeV in TPC gas
  set_default_double_param("electrons_per_gev", Tpc_ElectronsPerKeV * 1e6);

  // number of electrons deposited by laser per cm
  set_default_int_param("electrons_per_cm", 72);
}

//_____________________________________________________________
void PHG4TpcDirectLaser::SetPhiStepping(int n, double min, double max)
{
  if (n < 0 || max < min)
  {
    std::cout << PHWHERE << " - invalid" << std::endl;
    return;
  }
  nPhiSteps = n;
  minPhi = min;
  maxPhi = max;
  nTotalSteps = nThetaSteps * nPhiSteps;
  return;
}
//_____________________________________________________________
void PHG4TpcDirectLaser::SetThetaStepping(int n, double min, double max)
{
  if (n < 0 || max < min)
  {
    std::cout << PHWHERE << " - invalid" << std::endl;
    return;
  }
  nThetaSteps = n;
  minTheta = min;
  maxTheta = max;
  nTotalSteps = nThetaSteps * nPhiSteps;

  return;
}

//_____________________________________________________________
void PHG4TpcDirectLaser::SetupLasers()
{
  // clear previous lasers
  m_lasers.clear();

  // position of first laser at positive z
  const TVector3 position_base(60 * cm, 0., halflength_tpc);

  // add lasers
  for (int i = 0; i < 8; ++i)
  {
    Laser laser;

    // set laser direction
    /*
     * first four lasers are on positive z readout plane, and shoot towards negative z
     * next four lasers are on negative z readout plane and shoot towards positive z
     */
    laser.m_position = position_base;
    if (i < 4)
    {
      laser.m_position.SetZ(position_base.z());
      laser.m_direction = -1;
    }
    else
    {
      laser.m_position.SetZ(-position_base.z());
      laser.m_direction = 1;
    }

    // rotate around z
    laser.m_phi = M_PI / 2 * i;
    laser.m_position.RotateZ(laser.m_phi);

    // append
    m_lasers.push_back(laser);
  }
}

//_____________________________________________________________
void PHG4TpcDirectLaser::AimToNextPatternStep()
{
  if (nTotalSteps >= 1)
  {
    AimToPatternStep(currentPatternStep);
    ++currentPatternStep;
  }
}

//_____________________________________________________________
void PHG4TpcDirectLaser::AimToThetaPhi(double theta, double phi)
{
  if (Verbosity())
  {
    std::cout << "PHG4TpcDirectLaser::AimToThetaPhi - theta: " << theta << " phi: " << phi << std::endl;
  }

  // all lasers
  for (const auto& laser : m_lasers)
  {
    AppendLaserTrack(theta, phi, laser);
  }
}

//_____________________________________________________________
void PHG4TpcDirectLaser::AimToPatternStep(int n)
{
  //trim against overflows
  n = n % nTotalSteps;

  if (Verbosity())
  {
    std::cout << "PHG4TpcDirectLaser::AimToPatternStep - step: " << n << "/" << nTotalSteps << std::endl;
  }

  // store as current pattern
  currentPatternStep = n;

  // calculate theta
  const int thetaStep = n / nPhiSteps;
  const double theta = minTheta + thetaStep * (maxTheta - minTheta) / nThetaSteps;

  // calculate phi
  const int phiStep = n % nPhiSteps;
  const double phi = minPhi + phiStep * (maxPhi - minPhi) / nPhiSteps;

  // generate laser tracks
  AimToThetaPhi(theta, phi);

  return;
}

//_____________________________________________________________
void PHG4TpcDirectLaser::AppendLaserTrack(double theta, double phi, const PHG4TpcDirectLaser::Laser& laser)
{
  if (!m_g4hitcontainer)
  {
    std::cout << PHWHERE << "invalid g4hit container. aborting" << std::endl;
    return;
  }

  // store laser position
  const auto& pos = laser.m_position;

  // define track direction
  const auto& direction = laser.m_direction;
  TVector3 dir(0, 0, direction);

  //adjust direction
  dir.RotateY(theta * direction);
  dir.RotateZ(phi);

  // also rotate by laser azimuth
  dir.RotateZ(laser.m_phi);

  // print
  if (Verbosity())
  {
    std::cout << "PHG4TpcDirectLaser::AppendLaserTrack - position: " << pos << " direction: " << dir << std::endl;
  }

  // dummy momentum
  static constexpr double total_momentum = 1;

  // mc track id
  int trackid = -1;

  // create truth vertex and particle
  if (m_g4truthinfo)
  {
    // add vertex
    const auto vtxid = m_g4truthinfo->maxvtxindex() + 1;
    const auto vertex = new PHG4VtxPoint_t(pos.x(), pos.y(), pos.z(), 0, vtxid);
    m_g4truthinfo->AddVertex(vtxid, vertex);

    // increment track id
    trackid = m_g4truthinfo->maxtrkindex() + 1;

    // create new g4particle
    auto particle = new PHG4Particle_t();
    particle->set_track_id(trackid);
    particle->set_vtx_id(vtxid);
    particle->set_parent_id(0);
    particle->set_primary_id(trackid);
    particle->set_px(total_momentum * dir.x());
    particle->set_py(total_momentum * dir.y());
    particle->set_pz(total_momentum * dir.z());

    m_g4truthinfo->AddParticle(trackid, particle);
  }

  // store in SvtxTrack map
  if (m_track_map)
  {
    SvtxTrack_v2 track;
    track.set_x(pos.x());
    track.set_y(pos.y());
    track.set_z(pos.z());

    // total momentum is irrelevant. What matters is the direction
    track.set_px(total_momentum * dir.x());
    track.set_py(total_momentum * dir.y());
    track.set_pz(total_momentum * dir.z());

    // insert in map
    m_track_map->insert(&track);

    if (Verbosity())
    {
      std::cout << "PHG4TpcDirectLaser::AppendLaserTrack - position: " << pos << " direction: " << dir << std::endl;
    }
  }

  // find collision point
  /*
   * intersection to either central membrane or endcaps
   * if the position along beam and laser direction have the same sign, it will intercept the endcap
   * otherwise will intercept the central membrane
   */
  const auto plane_strike = (pos.z() * dir.z() > 0) ? endcap_intersection(pos, dir) : central_membrane_intersection(pos, dir);

  // field cage intersection
  const auto fc_strike = field_cage_intersection(pos, dir);

  // if none of the strikes is valid, there is no valid information found.
  if (!(plane_strike || fc_strike)) return;

  // decide relevant end of laser
  /* chose field cage intersection if valid, and if either plane intersection is invalid or happens on a larger z along the laser direction) */
  const TVector3& strike = (fc_strike && (!plane_strike || fc_strike->z() / dir.z() < plane_strike->z() / dir.z())) ? *fc_strike : *plane_strike;

  //find length
  TVector3 delta = (strike - pos);
  double fullLength = delta.Mag();
  int nHitSteps = fullLength / maxHitLength + 1;

  TVector3 start = pos;
  TVector3 end = start;
  TVector3 step = dir * (maxHitLength / (dir.Mag()));
  double stepLength = 0;

  if (Verbosity())
  {
    std::cout << "PHG4TpcDirectLaser::AppendLaserTrack -"
              << " fullLength: " << fullLength
              << " nHitSteps: " << nHitSteps
              << std::endl;
  }

  for (int i = 0; i < nHitSteps; i++)
  {
    start = end;  //new starting point is the previous ending point.
    if (i + 1 == nHitSteps)
    {
      //last step is the remainder size
      end = strike;
      delta = start - end;
      stepLength = delta.Mag();
    }
    else
    {
      //all other steps are uniform length
      end = start + step;
      stepLength = step.Mag();
    }

    //from phg4tpcsteppingaction.cc
    auto hit = new PHG4Hit_t;
    hit->set_trkid(trackid);
    hit->set_layer(99);

    //here we set the entrance values in cm
    hit->set_x(0, start.X() / cm);
    hit->set_y(0, start.Y() / cm);
    hit->set_z(0, start.Z() / cm);
    hit->set_t(0, (start - pos).Mag() / speed_of_light);

    hit->set_x(1, end.X() / cm);
    hit->set_y(1, end.Y() / cm);
    hit->set_z(1, end.Z() / cm);
    hit->set_t(1, (end - pos).Mag() / speed_of_light);

    // momentum
    hit->set_px(0, dir.X());  // GeV
    hit->set_py(0, dir.Y());
    hit->set_pz(0, dir.Z());

    hit->set_px(1, dir.X());
    hit->set_py(1, dir.Y());
    hit->set_pz(1, dir.Z());

    const double totalE = electrons_per_cm * stepLength / electrons_per_gev;

    hit->set_eion(totalE);
    hit->set_edep(totalE);
    m_g4hitcontainer->AddHit(detId, hit);
  }

  return;
}
