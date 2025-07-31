// $Id: $

/*!
 * \file PHHepMCGenHelper.cc
 * \brief
 * \author Jin Huang <jhuang@bnl.gov>
 * \version $Revision:   $
 * \date $Date: $
 */

#include "PHHepMCGenHelper.h"

#include "PHHepMCGenEvent.h"
#include "PHHepMCGenEventMap.h"
#include "PHHepMCGenEventv1.h"

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>    // for PHIODataNode
#include <phool/PHNode.h>          // for PHNode
#include <phool/PHNodeIterator.h>  // for PHNodeIterator
#include <phool/PHObject.h>        // for PHObject
#include <phool/PHRandomSeed.h>
#include <phool/getClass.h>
#include <phool/phool.h>  // for PHWHERE

#include <HepMC/SimpleVector.h>  // for FourVector

#include <CLHEP/Units/PhysicalConstants.h>
#include <CLHEP/Units/SystemOfUnits.h>
#include <CLHEP/Vector/Boost.h>
#include <CLHEP/Vector/LorentzRotation.h>
#include <CLHEP/Vector/LorentzVector.h>
#include <CLHEP/Vector/Rotation.h>

#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>

#include <cassert>
#include <cmath>
#include <cstdlib>  // for exit
#include <iostream>
#include <limits>

PHHepMCGenHelper::PHHepMCGenHelper()
  : RandomGenerator(gsl_rng_alloc(gsl_rng_mt19937))
{
  unsigned int seed = PHRandomSeed();  // fixed seed is handled in this function
  gsl_rng_set(RandomGenerator, seed);
}

PHHepMCGenHelper::~PHHepMCGenHelper()
{
  gsl_rng_free(RandomGenerator);
}

//! init interface nodes
int PHHepMCGenHelper::create_node_tree(PHCompositeNode *topNode)
{
  PHNodeIterator iter(topNode);
  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    std::cout << PHWHERE << "DST Node missing doing nothing" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  _geneventmap = findNode::getClass<PHHepMCGenEventMap>(dstNode, "PHHepMCGenEventMap");
  if (!_geneventmap)
  {
    _geneventmap = new PHHepMCGenEventMap();
    PHIODataNode<PHObject> *newmapnode = new PHIODataNode<PHObject>(_geneventmap, "PHHepMCGenEventMap", "PHObject");
    dstNode->addNode(newmapnode);
  }

  assert(_geneventmap);

  return Fun4AllReturnCodes::EVENT_OK;
}

//! choice of reference version of the PHHepMCGenEvent
const PHHepMCGenEvent *PHHepMCGenHelper::get_PHHepMCGenEvent_template()
{
  // choice of version of PHHepMCGenEvent
  const static PHHepMCGenEventv1 mc_evnet_template;

  return &mc_evnet_template;
}

//! send HepMC::GenEvent to DST tree. This function takes ownership of evt
PHHepMCGenEvent *PHHepMCGenHelper::insert_event(HepMC::GenEvent *evt)
{
  assert(evt);
  assert(_geneventmap);

  PHHepMCGenEvent *genevent = _geneventmap->insert_event(_embedding_id, get_PHHepMCGenEvent_template());

  genevent->addEvent(evt);
  HepMC2Lab_boost_rotation_translation(genevent);

  return genevent;
}

void PHHepMCGenHelper::move_vertex(PHHepMCGenEvent *genevent)
{
  assert(genevent);

  assert(_vertex_width_x >= 0);
  assert(_vertex_width_y >= 0);
  assert(_vertex_width_z >= 0);
  assert(_vertex_width_t >= 0);

  assert(not _reuse_vertex);  // logic check

  // not reusing vertex so smear with the vertex parameters
  genevent->moveVertex(
      (smear(_vertex_x, _vertex_width_x, _vertex_func_x)),
      (smear(_vertex_y, _vertex_width_y, _vertex_func_y)),
      (smear(_vertex_z, _vertex_width_z, _vertex_func_z)),
      (smear(_vertex_t, _vertex_width_t, _vertex_func_t)));
}

//! use m_beam_bunch_width to calculate horizontal and vertical collision width
//! \param[in] hv_index 0: horizontal. 1: vertical
//! https://github.com/eic/documents/blob/d06b5597a0a89dcad215bab50fe3eefa17a097a5/reports/general/Note-Simulations-BeamEffects.pdf
double PHHepMCGenHelper::get_collision_width(unsigned int hv_index)
{
  assert((hv_index == 0) or (hv_index == 1));

  const double widthA = m_beam_bunch_width.first[hv_index];
  const double widthB = m_beam_bunch_width.second[hv_index];

  return widthA * widthB / sqrt((widthA * widthA) + (widthB * widthB));
}

//! generate vertx with bunch interaction according to
//! https://github.com/eic/documents/blob/d06b5597a0a89dcad215bab50fe3eefa17a097a5/reports/general/Note-Simulations-BeamEffects.pdf
//! \return pair of bunch local z position for beam A and beam B
std::pair<double, double> PHHepMCGenHelper::generate_vertx_with_bunch_interaction(PHHepMCGenEvent *genevent)
{
  const std::pair<double, double> bunch_zs(
      smear(
          0,                            //  central vertical angle shift
          m_beam_bunch_width.first[2],  // vertical angle smear
          Gaus),
      smear(
          0,                             //  central vertical angle shift
          m_beam_bunch_width.second[2],  // vertical angle smear
          Gaus));

  CLHEP::Hep3Vector beamA_center = pair2Hep3Vector(m_beam_direction_theta_phi.first);
  CLHEP::Hep3Vector beamB_center = pair2Hep3Vector(m_beam_direction_theta_phi.second);
  //  const static CLHEP::Hep3Vector z_axis(0, 0, 1);
  const static CLHEP::Hep3Vector y_axis(0, 1, 0);

  // the final longitudinal vertex smear axis
  CLHEP::Hep3Vector beamCenterDiffAxis = (beamA_center - beamB_center);
  assert(beamCenterDiffAxis.mag() > CLHEP::Hep3Vector::getTolerance());
  beamCenterDiffAxis = beamCenterDiffAxis / beamCenterDiffAxis.mag();

  CLHEP::Hep3Vector vec_crossing = beamA_center - 0.5 * (beamA_center - beamB_center);

  CLHEP::Hep3Vector vec_longitudinal_collision = beamCenterDiffAxis * (bunch_zs.first + bunch_zs.second) / 2.;
  double ct_collision = 0.5 * (-bunch_zs.first + bunch_zs.second) / beamCenterDiffAxis.dot(beamA_center);
  double t_collision = ct_collision * CLHEP::cm / CLHEP::c_light / CLHEP::ns;
  CLHEP::Hep3Vector vec_crossing_collision = ct_collision * vec_crossing;  // shift of collision to crossing dierction

  CLHEP::Hep3Vector horizontal_axis = y_axis.cross(beamCenterDiffAxis);
  assert(horizontal_axis.mag() > CLHEP::Hep3Vector::getTolerance());
  horizontal_axis = horizontal_axis / horizontal_axis.mag();

  CLHEP::Hep3Vector vertical_axis = beamCenterDiffAxis.cross(horizontal_axis);
  assert(vertical_axis.mag() > CLHEP::Hep3Vector::getTolerance());
  vertical_axis = vertical_axis / vertical_axis.mag();

  CLHEP::Hep3Vector vec_horizontal_collision_vertex_smear = horizontal_axis *
                                                            smear(
                                                                0,
                                                                get_collision_width(0),
                                                                Gaus);
  CLHEP::Hep3Vector vec_vertical_collision_vertex_smear = vertical_axis *
                                                          smear(
                                                              0,
                                                              get_collision_width(1),
                                                              Gaus);

  CLHEP::Hep3Vector vec_collision_vertex =
      vec_horizontal_collision_vertex_smear +
      vec_vertical_collision_vertex_smear +  //
      vec_crossing_collision + vec_longitudinal_collision;

  genevent->set_collision_vertex(HepMC::FourVector(
      vec_collision_vertex.x(),
      vec_collision_vertex.y(),
      vec_collision_vertex.z(),
      t_collision));

  if (m_verbosity)
  {
    std::cout << __PRETTY_FUNCTION__
              << ":"
              << "bunch_zs.first  = " << bunch_zs.first << ", "
              << "bunch_zs.second = " << bunch_zs.second << ", "
              << "cos(theta/2) = " << beamCenterDiffAxis.dot(beamA_center) << ", " << std::endl

              << "beamCenterDiffAxis = " << beamCenterDiffAxis << ", "
              << "vec_crossing = " << vec_crossing << ", "
              << "horizontal_axis = " << horizontal_axis << ", "
              << "vertical_axis = " << vertical_axis << ", " << std::endl

              << "vec_longitudinal_collision = " << vec_longitudinal_collision << ", "
              << "vec_crossing_collision = " << vec_crossing_collision << ", "
              << "vec_vertical_collision_vertex_smear = " << vec_vertical_collision_vertex_smear << ", "
              << "vec_horizontal_collision_vertex_smear = " << vec_horizontal_collision_vertex_smear << ", " << std::endl
              << "vec_collision_vertex = " << vec_collision_vertex << ", " << std::endl

              << "ct_collision = " << ct_collision << ", "
              << "t_collision = " << t_collision << ", "
              << std::endl;
  }

  return bunch_zs;
}

CLHEP::Hep3Vector PHHepMCGenHelper::pair2Hep3Vector(const std::pair<double, double> &theta_phi)
{
  const double &theta = theta_phi.first;
  const double &phi = theta_phi.second;

  return CLHEP::Hep3Vector(
      sin(theta) * cos(phi),
      sin(theta) * sin(phi),
      cos(theta));
}

//! move vertex in translation,boost,rotation according to vertex settings
void PHHepMCGenHelper::HepMC2Lab_boost_rotation_translation(PHHepMCGenEvent *genevent)
{
  if (m_verbosity)
  {
    Print();
  }

  assert(genevent);

  if (_reuse_vertex)
  {
    // just copy over the vertex boost_rotation_translation

    assert(_geneventmap);

    PHHepMCGenEvent *vtx_evt =
        _geneventmap->get(_reuse_vertex_embedding_id);

    if (!vtx_evt)
    {
      std::cout << "PHHepMCGenHelper::HepMC2Lab_boost_rotation_translation - Fatal Error - the requested source subevent with embedding ID "
                << _reuse_vertex_embedding_id << " does not exist. Current HepMCEventMap:";
      _geneventmap->identify();
      exit(1);
    }

    // copy boost_rotation_translation

    genevent->moveVertex(
        vtx_evt->get_collision_vertex().x(),
        vtx_evt->get_collision_vertex().y(),
        vtx_evt->get_collision_vertex().z(),
        vtx_evt->get_collision_vertex().t());

    genevent->set_boost_beta_vector(vtx_evt->get_boost_beta_vector());
    genevent->set_rotation_vector(vtx_evt->get_rotation_vector());
    genevent->set_rotation_angle(vtx_evt->get_rotation_angle());

    if (m_verbosity)
    {
      std::cout << __PRETTY_FUNCTION__ << ": copied boost rotation shift of the collision" << std::endl;
      genevent->identify();
    }
    return;
  }  //!   if (_reuse_vertex)

  // now handle the collision vertex first, in the head-on collision frame
  // this is used as input to the Crab angle correction
  std::pair<double, double> beam_bunch_zs;
  if (m_use_beam_bunch_sim)
  {
    // bunch interaction simulation
    beam_bunch_zs = generate_vertx_with_bunch_interaction(genevent);
  }
  else
  {
    // vertex distribution simulation
    move_vertex(genevent);
    const double init_vertex_longitudinal = genevent->get_collision_vertex().z();
    beam_bunch_zs.first = beam_bunch_zs.second = init_vertex_longitudinal;
  }

  // boost-rotation from beam angles

  const static CLHEP::Hep3Vector z_axis(0, 0, 1);

  CLHEP::Hep3Vector beamA_center = pair2Hep3Vector(m_beam_direction_theta_phi.first);
  CLHEP::Hep3Vector beamB_center = pair2Hep3Vector(m_beam_direction_theta_phi.second);

  if (m_verbosity)
  {
    std::cout << __PRETTY_FUNCTION__ << ": " << std::endl;
    std::cout << "beamA_center = " << beamA_center << std::endl;
    std::cout << "beamB_center = " << beamB_center << std::endl;
  }

  assert(fabs(beamB_center.mag2() - 1) < CLHEP::Hep3Vector::getTolerance());
  assert(fabs(beamB_center.mag2() - 1) < CLHEP::Hep3Vector::getTolerance());

  if (beamA_center.dot(beamB_center) > -0.5)
  {
    std::cout << "PHHepMCGenHelper::HepMC2Lab_boost_rotation_translation - WARNING -"
              << "Beam A and Beam B are not near back to back. "
              << "Please double check beam direction setting at set_beam_direction_theta_phi()."
              << "beamA_center = " << beamA_center << ","
              << "beamB_center = " << beamB_center << ","
              << " Current setting:";

    Print();
  }

  // function to do angular shifts relative to central beam angle
  auto smear_beam_divergence = [&, this](
                                   const CLHEP::Hep3Vector &beam_center,
                                   const std::pair<double, double> &divergence_hv,
                                   const std::pair<double, double> &beam_angular_z_coefficient_hv)
  {
    const double &x_divergence = divergence_hv.first;
    const double &y_divergence = divergence_hv.second;

    // y direction in accelerator
    static const CLHEP::Hep3Vector accelerator_plane(0, 1, 0);

    //    CLHEP::Hep3Vector beam_direction(beam_center);
    CLHEP::HepRotation x_smear_in_accelerator_plane(
        accelerator_plane,
        smear(
            beam_bunch_zs.first * beam_angular_z_coefficient_hv.first,  //  central horizontal angle shift
            x_divergence,                                               // horizontal angle smear
            Gaus));
    CLHEP::HepRotation y_smear_out_accelerator_plane(
        accelerator_plane.cross(beam_center),
        smear(
            beam_bunch_zs.second * beam_angular_z_coefficient_hv.second,  //  central vertical angle shift
            y_divergence,                                                 // vertical angle smear
            Gaus));

    return y_smear_out_accelerator_plane * x_smear_in_accelerator_plane * beam_center;
  };

  CLHEP::Hep3Vector beamA_vec = smear_beam_divergence(beamA_center,
                                                      m_beam_angular_divergence_hv.first,
                                                      m_beam_angular_z_coefficient_hv.first);
  CLHEP::Hep3Vector beamB_vec = smear_beam_divergence(beamB_center,
                                                      m_beam_angular_divergence_hv.second,
                                                      m_beam_angular_z_coefficient_hv.second);

  if (m_verbosity)
  {
    std::cout << __PRETTY_FUNCTION__ << ": " << std::endl;
    std::cout << "beamA_vec = " << beamA_vec << std::endl;
    std::cout << "beamB_vec = " << beamB_vec << std::endl;
  }

  assert(fabs(beamA_vec.mag2() - 1) < CLHEP::Hep3Vector::getTolerance());
  assert(fabs(beamB_vec.mag2() - 1) < CLHEP::Hep3Vector::getTolerance());

  // apply minimal beam energy shift rotation and boost
  CLHEP::Hep3Vector boost_axis = beamA_vec + beamB_vec;
  if (boost_axis.mag2() > CLHEP::Hep3Vector::getTolerance())
  {
    // non-zero boost

    // split the boost to half for each beam for minimal beam  energy shift
    genevent->set_boost_beta_vector(0.5 * boost_axis);

    if (m_verbosity)
    {
      std::cout << __PRETTY_FUNCTION__ << ": non-zero boost " << std::endl;
    }
  }  //    if (cos_rotation_angle> CLHEP::Hep3Vector::getTolerance())
  else
  {
    genevent->set_boost_beta_vector(CLHEP::Hep3Vector(0, 0, 0));
    if (m_verbosity)
    {
      std::cout << __PRETTY_FUNCTION__ << ": zero boost " << std::endl;
    }
  }

  // rotation to collision to along z-axis with beamA pointing to +z
  CLHEP::Hep3Vector beamDiffAxis = (beamA_vec - beamB_vec);
  if (beamDiffAxis.mag2() < CLHEP::Hep3Vector::getTolerance())
  {
    std::cout << "PHHepMCGenHelper::HepMC2Lab_boost_rotation_translation - Fatal error -"
              << "Beam A and Beam B are too close to each other in direction "
              << "Please double check beam direction and divergence setting. "
              << "beamA_vec = " << beamA_vec << ","
              << "beamB_vec = " << beamB_vec << ","
              << " Current setting:";

    Print();

    exit(1);
  }

  beamDiffAxis = beamDiffAxis / beamDiffAxis.mag();
  double cos_rotation_angle_to_z = beamDiffAxis.dot(z_axis);
  if (m_verbosity)
  {
    std::cout << __PRETTY_FUNCTION__ << ": check rotation ";
    std::cout << "cos_rotation_angle_to_z= " << cos_rotation_angle_to_z << std::endl;
  }

  if (1 - cos_rotation_angle_to_z < CLHEP::Hep3Vector::getTolerance())
  {
    // no rotation
    genevent->set_rotation_vector(z_axis);
    genevent->set_rotation_angle(0);

    if (m_verbosity)
    {
      std::cout << __PRETTY_FUNCTION__ << ": no rotation " << std::endl;
    }
  }
  else if (cos_rotation_angle_to_z + 1 < CLHEP::Hep3Vector::getTolerance())
  {
    // you got beam flipped
    genevent->set_rotation_vector(CLHEP::Hep3Vector(0, 1, 0));
    genevent->set_rotation_angle(M_PI);
    if (m_verbosity)
    {
      std::cout << __PRETTY_FUNCTION__ << ": reverse beam direction " << std::endl;
    }
  }
  else
  {
    // need a rotation
    CLHEP::Hep3Vector rotation_axis = (beamA_vec - beamB_vec).cross(z_axis);
    const double rotation_angle_to_z = -acos(cos_rotation_angle_to_z);

    genevent->set_rotation_vector(rotation_axis);
    genevent->set_rotation_angle(rotation_angle_to_z);

    if (m_verbosity)
    {
      std::cout << __PRETTY_FUNCTION__ << ": has rotation " << std::endl;
    }
  }  //  if (boost_axis.mag2() > CLHEP::Hep3Vector::getTolerance())

  if (m_verbosity)
  {
    std::cout << __PRETTY_FUNCTION__ << ": final boost rotation shift of the collision" << std::endl;
    genevent->identify();
  }
}

void PHHepMCGenHelper::set_vertex_distribution_function(VTXFUNC x, VTXFUNC y, VTXFUNC z, VTXFUNC t)
{
  if (m_use_beam_bunch_sim)
  {
    std::cout << __PRETTY_FUNCTION__ << " Fatal Error: "
              << "m_use_beam_bunch_sim = " << m_use_beam_bunch_sim << ". Expect to simulate bunch interaction instead of applying vertex distributions"
              << std::endl;
    exit(1);
  }
  _vertex_func_x = x;
  _vertex_func_y = y;
  _vertex_func_z = z;
  _vertex_func_t = t;
  return;
}

void PHHepMCGenHelper::set_vertex_distribution_mean(const double x, const double y, const double z, const double t)
{
  if (m_use_beam_bunch_sim)
  {
    std::cout << __PRETTY_FUNCTION__ << " Fatal Error: "
              << "m_use_beam_bunch_sim = " << m_use_beam_bunch_sim << ". Expect to simulate bunch interaction instead of applying vertex distributions"
              << std::endl;
    exit(1);
  }

  _vertex_x = x;
  _vertex_y = y;
  _vertex_z = z;
  _vertex_t = t;
  return;
}

void PHHepMCGenHelper::set_vertex_distribution_width(const double x, const double y, const double z, const double t)
{
  if (m_use_beam_bunch_sim)
  {
    std::cout << __PRETTY_FUNCTION__ << " Fatal Error: "
              << "m_use_beam_bunch_sim = " << m_use_beam_bunch_sim << ". Expect to simulate bunch interaction instead of applying vertex distributions"
              << std::endl;
    exit(1);
  }

  _vertex_width_x = x;
  _vertex_width_y = y;
  _vertex_width_z = z;
  _vertex_width_t = t;
  return;
}

void PHHepMCGenHelper::set_beam_bunch_width(const std::vector<double> &beamA, const std::vector<double> &beamB)
{
  if (not m_use_beam_bunch_sim)
  {
    std::cout << __PRETTY_FUNCTION__ << " Fatal Error: "
              << "m_use_beam_bunch_sim = " << m_use_beam_bunch_sim << ". Expect not to simulate bunch interaction but applying vertex distributions"
              << std::endl;
    exit(1);
  }

  m_beam_bunch_width.first = beamA;
  m_beam_bunch_width.second = beamB;
}

double PHHepMCGenHelper::smear(const double position,
                               const double width,
                               VTXFUNC dist) const
{
  assert(width >= 0);

  double res = position;

  if (width == 0)
  {
    return res;
  }

  if (dist == Uniform)
  {
    res = (position - width) + 2 * gsl_rng_uniform_pos(RandomGenerator) * width;
  }
  else if (dist == Gaus)
  {
    res = position + gsl_ran_gaussian(RandomGenerator, width);
  }
  else
  {
    std::cout << "PHHepMCGenHelper::smear - FATAL Error - unknown vertex function " << dist << std::endl;
    exit(10);
  }
  return res;
}

void PHHepMCGenHelper::CopySettings(PHHepMCGenHelper &helper_dest)
{
  // allow copy of vertex distributions
  helper_dest.use_beam_bunch_sim(false);
  helper_dest.set_vertex_distribution_width(_vertex_width_x, _vertex_width_y, _vertex_width_z, _vertex_width_t);
  helper_dest.set_vertex_distribution_function(_vertex_func_x, _vertex_func_y, _vertex_func_z, _vertex_func_t);
  helper_dest.set_vertex_distribution_mean(_vertex_x, _vertex_y, _vertex_z, _vertex_t);

  // allow copy of bunch distributions
  helper_dest.use_beam_bunch_sim(true);
  helper_dest.set_beam_bunch_width(m_beam_bunch_width.first, m_beam_bunch_width.second);

  // final bunch settings
  helper_dest.use_beam_bunch_sim(m_use_beam_bunch_sim);

  helper_dest.set_beam_direction_theta_phi(
      m_beam_direction_theta_phi.first.first,
      m_beam_direction_theta_phi.first.second,
      m_beam_direction_theta_phi.second.first,
      m_beam_direction_theta_phi.second.second);
  helper_dest.set_beam_angular_divergence_hv(
      m_beam_angular_divergence_hv.first.first,
      m_beam_angular_divergence_hv.first.second,
      m_beam_angular_divergence_hv.second.first,
      m_beam_angular_divergence_hv.second.second);
  helper_dest.set_beam_angular_z_coefficient_hv(
      m_beam_angular_z_coefficient_hv.first.first,
      m_beam_angular_z_coefficient_hv.first.second,
      m_beam_angular_z_coefficient_hv.second.first,
      m_beam_angular_z_coefficient_hv.second.second);

  helper_dest.set_beam_angular_z_coefficient_hv(
      m_beam_angular_z_coefficient_hv.first.first,
      m_beam_angular_z_coefficient_hv.first.second,
      m_beam_angular_z_coefficient_hv.second.first,
      m_beam_angular_z_coefficient_hv.second.second);

  return;
}

void PHHepMCGenHelper::CopySettings(PHHepMCGenHelper *helper_dest)
{
  if (helper_dest)
  {
    CopySettings(*helper_dest);
  }
  else
  {
    std::cout << "PHHepMCGenHelper::CopySettings - fatal error - invalid input class helper_dest which is nullptr!" << std::endl;
    exit(1);
  }
}

void PHHepMCGenHelper::CopyHelperSettings(PHHepMCGenHelper *helper_src)
{
  if (helper_src)
  {
    helper_src->CopySettings(this);
  }
  else
  {
    std::cout << "PHHepMCGenHelper::CopyHelperSettings - fatal error - invalid input class helper_src which is nullptr!" << std::endl;
    exit(1);
  }
}

void PHHepMCGenHelper::Print(const std::string & /*what*/) const
{
  static std::map<VTXFUNC, std::string> vtxfunc = {{VTXFUNC::Uniform, "Uniform"}, {VTXFUNC::Gaus, "Gaus"}};

  std::cout << "Vertex distribution width x: " << _vertex_width_x
            << ", y: " << _vertex_width_y
            << ", z: " << _vertex_width_z
            << ", t: " << _vertex_width_t
            << std::endl;

  std::cout << "Vertex distribution function x: " << vtxfunc[_vertex_func_x]
            << ", y: " << vtxfunc[_vertex_func_y]
            << ", z: " << vtxfunc[_vertex_func_z]
            << ", t: " << vtxfunc[_vertex_func_t]
            << std::endl;

  std::cout << "Beam direction: A  theta-phi = " << m_beam_direction_theta_phi.first.first
            << ", " << m_beam_direction_theta_phi.first.second << std::endl;
  std::cout << "Beam direction: B  theta-phi = " << m_beam_direction_theta_phi.second.first
            << ", " << m_beam_direction_theta_phi.second.second << std::endl;

  std::cout << "Beam divergence: A X-Y = " << m_beam_angular_divergence_hv.first.first
            << ", " << m_beam_angular_divergence_hv.first.second << std::endl;
  std::cout << "Beam divergence: B X-Y = " << m_beam_angular_divergence_hv.second.first
            << ", " << m_beam_angular_divergence_hv.second.second << std::endl;

  std::cout << "Beam angle shift as linear function of longitudinal vertex position : A X-Y = " << m_beam_angular_z_coefficient_hv.first.first
            << ", " << m_beam_angular_z_coefficient_hv.first.second << std::endl;
  std::cout << "Beam angle shift as linear function of longitudinal vertex position: B X-Y = " << m_beam_angular_z_coefficient_hv.second.first
            << ", " << m_beam_angular_z_coefficient_hv.second.second << std::endl;

  std::cout << "m_use_beam_bunch_sim = " << m_use_beam_bunch_sim << std::endl;

  std::cout << "Beam bunch A width = ["
            << m_beam_bunch_width.first[0] << ", " << m_beam_bunch_width.first[1] << ", " << m_beam_bunch_width.first[2] << "] cm" << std::endl;
  std::cout << "Beam bunch B width = ["
            << m_beam_bunch_width.second[0] << ", " << m_beam_bunch_width.second[1] << ", " << m_beam_bunch_width.second[2] << "] cm" << std::endl;

  return;
}
