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

#include <CLHEP/Vector/Boost.h>
#include <CLHEP/Vector/LorentzRotation.h>
#include <CLHEP/Vector/LorentzVector.h>
#include <CLHEP/Vector/Rotation.h>
#include <CLHEP/Vector/ThreeVector.h>

#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>

#include <cassert>
#include <cmath>
#include <cstdlib>  // for exit
#include <iostream>
#include <limits>

using namespace std;

PHHepMCGenHelper::PHHepMCGenHelper()
  : _vertex_func_x(Gaus)
  , _vertex_func_y(Gaus)
  , _vertex_func_z(Gaus)
  , _vertex_func_t(Gaus)
  , _vertex_x(0)
  , _vertex_y(0)
  , _vertex_z(0)
  , _vertex_t(0)
  , _vertex_width_x(0)
  , _vertex_width_y(0)
  , _vertex_width_z(0)
  , _vertex_width_t(0)
  , _embedding_id(0)
  , _reuse_vertex(false)
  , _reuse_vertex_embedding_id(numeric_limits<int>::min())
  , _geneventmap(nullptr)
{
  RandomGenerator = gsl_rng_alloc(gsl_rng_mt19937);
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
    cout << PHWHERE << "DST Node missing doing nothing" << endl;
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
const PHHepMCGenEvent *PHHepMCGenHelper::get_PHHepMCGenEvent_template() const
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

  if (_reuse_vertex)
  {
    assert(_geneventmap);

    PHHepMCGenEvent *vtx_evt =
        _geneventmap->get(_reuse_vertex_embedding_id);

    if (!vtx_evt)
    {
      cout << "PHHepMCGenHelper::move_vertex - Fatal Error - the requested source subevent with embedding ID "
           << _reuse_vertex_embedding_id << " does not exist. Current HepMCEventMap:";
      _geneventmap->identify();
      exit(11);
    }

    genevent->moveVertex(
        vtx_evt->get_collision_vertex().x(),
        vtx_evt->get_collision_vertex().y(),
        vtx_evt->get_collision_vertex().z(),
        vtx_evt->get_collision_vertex().t());
  }
  else
  {
    // not reusing vertex so smear with the vertex parameters
    genevent->moveVertex(
        (smear(_vertex_x, _vertex_width_x, _vertex_func_x)),
        (smear(_vertex_y, _vertex_width_y, _vertex_func_y)),
        (smear(_vertex_z, _vertex_width_z, _vertex_func_z)),
        (smear(_vertex_t, _vertex_width_t, _vertex_func_t)));
  }
}

//! move vertex in translation,boost,rotation according to vertex settings
void PHHepMCGenHelper::HepMC2Lab_boost_rotation_translation(PHHepMCGenEvent *genevent)
{
  if (m_verbosity)
  {
    Print();
  }

  assert(genevent);

  // now handle the collision vertex first, in the head-on collision frame
  // this is used as input to the Crab angle correction
  move_vertex(genevent);
  const double init_vertex_longitudinal = genevent->get_collision_vertex().z();

  // boost-rotation from beam angles

  const static CLHEP::Hep3Vector z_axis(0, 0, 1);

  // function to convert spherical coordinate to Hep3Vector in x-y-z
  auto pair2Hep3Vector = [](const std::pair<double, double> &theta_phi) {
    const double &theta = theta_phi.first;
    const double &phi = theta_phi.second;

    return CLHEP::Hep3Vector(
        sin(theta) * cos(phi),
        sin(theta) * sin(phi),
        cos(theta));
  };

  CLHEP::Hep3Vector beamA_center = pair2Hep3Vector(m_beam_direction_theta_phi.first);
  CLHEP::Hep3Vector beamB_center = pair2Hep3Vector(m_beam_direction_theta_phi.second);

  if (m_verbosity)
  {
    cout << __PRETTY_FUNCTION__ << ": " << endl;
    cout << "beamA_center = " << beamA_center << endl;
    cout << "beamB_center = " << beamB_center << endl;
  }

  assert(fabs(beamB_center.mag2() - 1) < CLHEP::Hep3Vector::getTolerance());
  assert(fabs(beamB_center.mag2() - 1) < CLHEP::Hep3Vector::getTolerance());

  if (beamA_center.dot(beamB_center) > -0.5)
  {
    cout << "PHHepMCGenHelper::HepMC2Lab_boost_rotation_translation - WARNING -"
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
                                   const std::pair<double, double> &beam_angular_z_coefficient_hv) {
    const double &x_divergence = divergence_hv.first;
    const double &y_divergence = divergence_hv.second;

    // y direction in accelerator
    static const CLHEP::Hep3Vector accelerator_plane(0, 1, 0);

    CLHEP::Hep3Vector beam_direction(beam_center);
    CLHEP::HepRotation x_smear_in_accelerator_plane(
        accelerator_plane,
        smear(
            init_vertex_longitudinal * beam_angular_z_coefficient_hv.first,  //  central horizontal angle shift
            x_divergence,                                                    // horizontal angle smear
            Gaus));
    CLHEP::HepRotation y_smear_out_accelerator_plane(
        accelerator_plane.cross(beam_center),
        smear(
            init_vertex_longitudinal * beam_angular_z_coefficient_hv.second,  //  central vertical angle shift
            y_divergence,                                                     // vertical angle smear
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
    cout << __PRETTY_FUNCTION__ << ": " << endl;
    cout << "beamA_vec = " << beamA_vec << endl;
    cout << "beamB_vec = " << beamB_vec << endl;
  }

  assert(fabs(beamA_vec.mag2() - 1) < CLHEP::Hep3Vector::getTolerance());
  assert(fabs(beamB_vec.mag2() - 1) < CLHEP::Hep3Vector::getTolerance());

  // apply minimal beam energy shift rotation and boost
  CLHEP::Hep3Vector boost_axis = beamA_vec + beamB_vec;
  if (boost_axis.mag2() > CLHEP::Hep3Vector::getTolerance())
  {
    //non-zero boost

    // split the boost to half for each beam for minimal beam  energy shift
    genevent->set_boost_beta_vector(0.5 * boost_axis);

    if (m_verbosity)
    {
      cout << __PRETTY_FUNCTION__ << ": non-zero boost " << endl;
    }
  }  //    if (cos_rotation_angle> CLHEP::Hep3Vector::getTolerance())
  else
  {
    genevent->set_boost_beta_vector(CLHEP::Hep3Vector(0, 0, 0));
    if (m_verbosity)
    {
      cout << __PRETTY_FUNCTION__ << ": zero boost " << endl;
    }
  }

  //rotation to collision to along z-axis with beamA pointing to +z
  CLHEP::Hep3Vector beamDiffAxis = (beamA_vec - beamB_vec);
  if (beamDiffAxis.mag2() < CLHEP::Hep3Vector::getTolerance())
  {
    cout << "PHHepMCGenHelper::HepMC2Lab_boost_rotation_translation - Fatal error -"
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
    cout << __PRETTY_FUNCTION__ << ": check rotation ";
    cout << "cos_rotation_angle_to_z= " << cos_rotation_angle_to_z << endl;
  }

  if (1 - cos_rotation_angle_to_z < CLHEP::Hep3Vector::getTolerance())
  {
    //no rotation
    genevent->set_rotation_vector(z_axis);
    genevent->set_rotation_angle(0);

    if (m_verbosity)
    {
      cout << __PRETTY_FUNCTION__ << ": no rotation " << endl;
    }
  }
  else if (cos_rotation_angle_to_z + 1 < CLHEP::Hep3Vector::getTolerance())
  {
    // you got beam flipped
    genevent->set_rotation_vector(CLHEP::Hep3Vector(0, 1, 0));
    genevent->set_rotation_angle(M_PI);
    if (m_verbosity)
    {
      cout << __PRETTY_FUNCTION__ << ": reverse beam direction " << endl;
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
      cout << __PRETTY_FUNCTION__ << ": has rotation " << endl;
    }
  }  //  if (boost_axis.mag2() > CLHEP::Hep3Vector::getTolerance())

  // rotate the collision vertex z direction to middle of the beam angles
  if (not _reuse_vertex)
  {
    // the final longitudinal vertex smear axis
    CLHEP::Hep3Vector beamCenterDiffAxis = (beamA_center - beamB_center);
    beamCenterDiffAxis = beamCenterDiffAxis / beamCenterDiffAxis.mag();

    double cos_rotation_center_angle_to_z = beamCenterDiffAxis.dot(z_axis);

    if (1 - fabs(cos_rotation_center_angle_to_z) < CLHEP::Hep3Vector::getTolerance())
    {
      // new axis is basically beam axis

      if (m_verbosity)
      {
        cout << __PRETTY_FUNCTION__
             << ": collision longitudinal axis is very close to z-axis. No additional rotation of vertexes: "
             << "cos_rotation_center_angle_to_z = " << cos_rotation_center_angle_to_z
             << endl;
      }
      else
      {
        // need a rotation
        CLHEP::Hep3Vector rotation_axis = beamCenterDiffAxis.cross(z_axis);
        const double rotation_angle_to_z = -acos(cos_rotation_center_angle_to_z);
        const CLHEP::HepRotation rotation(rotation_axis, rotation_angle_to_z);

        const HepMC::FourVector init_4vertex = genevent->get_collision_vertex();
        CLHEP::Hep3Vector init_3vertex(
            init_4vertex.x(),
            init_4vertex.y(),
            init_4vertex.z());

        CLHEP::Hep3Vector final_3vertex = rotation * init_3vertex;

        genevent->set_collision_vertex(HepMC::FourVector(
            final_3vertex.x(),
            final_3vertex.y(),
            final_3vertex.z(),
            init_4vertex.t()));

        if (m_verbosity)
        {
          cout << __PRETTY_FUNCTION__
               << ": collision longitudinal axis is rotated: "
               << "cos_rotation_center_angle_to_z = " << cos_rotation_center_angle_to_z << ", "
               << "rotation_axis = " << rotation_axis << ", "
               << "init_3vertex = " << init_3vertex << ", "
               << "final_3vertex = " << final_3vertex << ", "
               << endl;
        }
      }
    }
  }

  if (m_verbosity)
  {
    cout << __PRETTY_FUNCTION__ << ": final boost rotation shift of the collision" << endl;
    genevent->identify();
  }
}

void PHHepMCGenHelper::set_vertex_distribution_function(VTXFUNC x, VTXFUNC y, VTXFUNC z, VTXFUNC t)
{
  _vertex_func_x = x;
  _vertex_func_y = y;
  _vertex_func_z = z;
  _vertex_func_t = t;
  return;
}

void PHHepMCGenHelper::set_vertex_distribution_mean(const double x, const double y, const double z, const double t)
{
  _vertex_x = x;
  _vertex_y = y;
  _vertex_z = z;
  _vertex_t = t;
  return;
}

void PHHepMCGenHelper::set_vertex_distribution_width(const double x, const double y, const double z, const double t)
{
  _vertex_width_x = x;
  _vertex_width_y = y;
  _vertex_width_z = z;
  _vertex_width_t = t;
  return;
}

double PHHepMCGenHelper::smear(const double position,
                               const double width,
                               VTXFUNC dist) const
{
  assert(width >= 0);

  double res = position;

  if (width == 0)
    return res;

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
    cout << "PHHepMCGenHelper::smear - FATAL Error - unknown vertex function " << dist << endl;
    exit(10);
  }
  return res;
}

void PHHepMCGenHelper::CopySettings(PHHepMCGenHelper &helper_dest)
{
  helper_dest.set_vertex_distribution_width(_vertex_width_x, _vertex_width_y, _vertex_width_z, _vertex_width_t);
  helper_dest.set_vertex_distribution_function(_vertex_func_x, _vertex_func_y, _vertex_func_z, _vertex_func_t);
  helper_dest.set_vertex_distribution_mean(_vertex_x, _vertex_y, _vertex_z, _vertex_t);

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

  return;
}

void PHHepMCGenHelper::CopySettings(PHHepMCGenHelper *helper_dest)
{
  if (helper_dest)
    CopySettings(*helper_dest);
  else
  {
    cout << "PHHepMCGenHelper::CopySettings - fatal error - invalid input class helper_dest which is nullptr!" << endl;
    exit(1);
  }
}

void PHHepMCGenHelper::CopyHelperSettings(PHHepMCGenHelper *helper_src)
{
  if (helper_src)
    helper_src->CopySettings(this);
  else
  {
    cout << "PHHepMCGenHelper::CopyHelperSettings - fatal error - invalid input class helper_src which is nullptr!" << endl;
    exit(1);
  }
}

void PHHepMCGenHelper::Print(const std::string &what) const
{
  static map<VTXFUNC, string> vtxfunc = {{VTXFUNC::Uniform, "Uniform"}, {VTXFUNC::Gaus, "Gaus"}};

  cout << "Vertex distribution width x: " << _vertex_width_x
       << ", y: " << _vertex_width_y
       << ", z: " << _vertex_width_z
       << ", t: " << _vertex_width_t
       << endl;

  cout << "Vertex distribution function x: " << vtxfunc[_vertex_func_x]
       << ", y: " << vtxfunc[_vertex_func_y]
       << ", z: " << vtxfunc[_vertex_func_z]
       << ", t: " << vtxfunc[_vertex_func_t]
       << endl;

  cout << "Beam direction: A  theta-phi = " << m_beam_direction_theta_phi.first.first
       << ", " << m_beam_direction_theta_phi.first.second << endl;
  cout << "Beam direction: B  theta-phi = " << m_beam_direction_theta_phi.second.first
       << ", " << m_beam_direction_theta_phi.second.second << endl;

  cout << "Beam divergence: A X-Y = " << m_beam_angular_divergence_hv.first.first
       << ", " << m_beam_angular_divergence_hv.first.second << endl;
  cout << "Beam divergence: B X-Y = " << m_beam_angular_divergence_hv.second.first
       << ", " << m_beam_angular_divergence_hv.second.second << endl;

  cout << "Beam angle shift as linear function of longitudinal vertex position : A X-Y = " << m_beam_angular_z_coefficient_hv.first.first
       << ", " << m_beam_angular_z_coefficient_hv.first.second << endl;
  cout << "Beam angle shift as linear function of longitudinal vertex position: B X-Y = " << m_beam_angular_z_coefficient_hv.second.first
       << ", " << m_beam_angular_z_coefficient_hv.second.second << endl;

  return;
}
