#include "PHG4SimpleEventGenerator.h"

#include "PHG4InEvent.h"
#include "PHG4Particle.h"  // for PHG4Particle
#include "PHG4Particlev2.h"
#include "PHG4Utils.h"

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHDataNode.h>      // for PHDataNode
#include <phool/PHNode.h>          // for PHNode
#include <phool/PHNodeIterator.h>  // for PHNodeIterator
#include <phool/PHObject.h>        // for PHObject
#include <phool/getClass.h>
#include <phool/phool.h>  // for PHWHERE

#include <TSystem.h>

#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>  // for gsl_rng_uniform_pos

#include <cassert>
#include <cmath>
#include <cstdlib>
#include <iostream>  // for operator<<, endl, basic_ostream
#include <memory>    // for allocator_traits<>::value_type

PHG4SimpleEventGenerator::PHG4SimpleEventGenerator(const std::string &name)
  : PHG4ParticleGeneratorBase(name)
{
  return;
}

void PHG4SimpleEventGenerator::add_particles(const std::string &name, const unsigned int num)
{
  _particle_names.push_back(std::make_pair(name, num));
  return;
}

void PHG4SimpleEventGenerator::add_particles(const int pid, const unsigned int num)
{
  _particle_codes.push_back(std::make_pair(pid, num));
  return;
}

void PHG4SimpleEventGenerator::set_eta_range(const double min, const double max)
{
  if (min > max)
  {
    std::cout << "not setting eta bc etamin " << min << " > etamax: " << max << std::endl;
    gSystem->Exit(1);
  }
  m_EtaMin = min;
  m_EtaMax = max;
  m_ThetaMin = NAN;
  m_ThetaMax = NAN;
  return;
}

void PHG4SimpleEventGenerator::set_theta_range(const double min, const double max)
{
  if (min > max)
  {
    std::cout << __PRETTY_FUNCTION__ << " thetamin " << min << " > thetamax: " << max << std::endl;
    gSystem->Exit(1);
  }
  if (min < 0 || max > M_PI)
  {
    std::cout << __PRETTY_FUNCTION__ << " min or max outside range (range is 0 to pi) min: " << min << ", max: " << max << std::endl;
    gSystem->Exit(1);
  }
  m_ThetaMin = min;
  m_ThetaMax = max;
  m_EtaMin = NAN;
  m_EtaMax = NAN;
  return;
}

void PHG4SimpleEventGenerator::set_phi_range(const double min, const double max)
{
  if (min > max)
  {
    std::cout << __PRETTY_FUNCTION__ << " phimin " << min << " > phimax: " << max << std::endl;
    gSystem->Exit(1);
    return;
  }
  if (min < -M_PI || max > M_PI)
  {
    std::cout << __PRETTY_FUNCTION__ << "min or max outside range (range is -pi to pi), min: " << min << ", max: " << max << std::endl;
    gSystem->Exit(1);
  }

  m_PhiMin = min;
  m_PhiMax = max;
  return;
}

void PHG4SimpleEventGenerator::set_power_law_n(const double n)
{
  m_powerLawN = n;
}

void PHG4SimpleEventGenerator::set_pt_range(const double min, const double max, const double pt_gaus_width)
{
  if (min > max)
  {
    std::cout << __PRETTY_FUNCTION__ << " ptmin " << min << " > ptmax: " << max << std::endl;
    gSystem->Exit(1);
  }
  if (min < 0 || max < 0 || pt_gaus_width < 0)
  {
    std::cout << __PRETTY_FUNCTION__ << " values need to be >= 0, min: " << min
              << ", max: " << max << ", pt_gaus_width: " << pt_gaus_width << std::endl;
    gSystem->Exit(1);
  }

  m_Pt_Min = min;
  m_Pt_Max = max;
  m_Pt_GausWidth = pt_gaus_width;
  m_P_Min = NAN;
  m_P_Max = NAN;
  m_P_GausWidth = NAN;
  return;
}

void PHG4SimpleEventGenerator::set_p_range(const double min, const double max, const double p_gaus_width)
{
  if (min > max)
  {
    std::cout << __PRETTY_FUNCTION__ << " pmin " << min << " > pmax: " << max << std::endl;
    gSystem->Exit(1);
  }
  if (min < 0 || max < 0 || p_gaus_width < 0)
  {
    std::cout << __PRETTY_FUNCTION__ << " values need to be >= 0, min: " << min
              << ", max: " << max << ", p_gaus_width: " << p_gaus_width << std::endl;
    gSystem->Exit(1);
  }
  m_Pt_Min = NAN;
  m_Pt_Max = NAN;
  m_Pt_GausWidth = NAN;
  m_P_Min = min;
  m_P_Max = max;
  m_P_GausWidth = p_gaus_width;
  return;
}

void PHG4SimpleEventGenerator::set_vertex_distribution_function(FUNCTION x, FUNCTION y, FUNCTION z)
{
  m_VertexFunc_x = x;
  m_VertexFunc_y = y;
  m_VertexFunc_z = z;
  return;
}

void PHG4SimpleEventGenerator::set_vertex_distribution_mean(const double x, const double y, const double z)
{
  m_Vertex_x = x;
  m_Vertex_y = y;
  m_Vertex_z = z;
  return;
}

void PHG4SimpleEventGenerator::set_vertex_distribution_width(const double x, const double y, const double z)
{
  m_VertexWidth_x = x;
  m_VertexWidth_y = y;
  m_VertexWidth_z = z;
  return;
}

void PHG4SimpleEventGenerator::set_existing_vertex_offset_vector(const double x, const double y, const double z)
{
  m_VertexOffset_x = x;
  m_VertexOffset_y = y;
  m_VertexOffset_z = z;
  return;
}

void PHG4SimpleEventGenerator::set_vertex_size_parameters(const double mean, const double width)
{
  m_VertexSizeMean = mean;
  m_VertexSizeWidth = width;
  return;
}

int PHG4SimpleEventGenerator::InitRun(PHCompositeNode *topNode)
{
  if (m_FunctionNames.find(m_VertexFunc_x) == m_FunctionNames.end())
  {
    std::cout << PHWHERE << "::Error - unknown x vertex distribution function requested" << std::endl;
    gSystem->Exit(1);
  }
  if (m_FunctionNames.find(m_VertexFunc_y) == m_FunctionNames.end())
  {
    std::cout << PHWHERE << "::Error - unknown y vertex distribution function requested" << std::endl;
    gSystem->Exit(1);
  }
  if (m_FunctionNames.find(m_VertexFunc_z) == m_FunctionNames.end())
  {
    std::cout << PHWHERE << "::Error - unknown z vertex distribution function requested" << std::endl;
    gSystem->Exit(1);
  }

  m_InEvent = findNode::getClass<PHG4InEvent>(topNode, "PHG4INEVENT");
  if (!m_InEvent)
  {
    PHNodeIterator iter(topNode);
    PHCompositeNode *dstNode;
    dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));

    m_InEvent = new PHG4InEvent();
    PHDataNode<PHObject> *newNode = new PHDataNode<PHObject>(m_InEvent, "PHG4INEVENT", "PHObject");
    dstNode->addNode(newNode);
  }

  if (Verbosity() > 0)
  {
    std::cout << "================ PHG4SimpleEventGenerator::InitRun() ======================" << std::endl;
    std::cout << " Random seed = " << get_seed() << std::endl;
    std::cout << " Particles:" << std::endl;
    for (unsigned int i = 0; i < _particle_codes.size(); ++i)
    {
      std::cout << "    " << _particle_codes[i].first << ", count = " << _particle_codes[i].second << std::endl;
    }
    for (unsigned int i = 0; i < _particle_names.size(); ++i)
    {
      std::cout << "    " << _particle_names[i].first << ", count = " << _particle_names[i].second << std::endl;
    }
    if (get_reuse_existing_vertex())
    {
      std::cout << " Vertex Distribution: Set to reuse a previously generated sim vertex" << std::endl;
      std::cout << " Vertex offset vector (x,y,z) = (" << m_VertexOffset_x << "," << m_VertexOffset_y << "," << m_VertexOffset_z << ")" << std::endl;
    }
    else
    {
      std::cout << " Vertex Distribution Function (x,y,z) = ("
                << m_FunctionNames.find(m_VertexFunc_x)->second << ","
                << m_FunctionNames.find(m_VertexFunc_y)->second << ","
                << m_FunctionNames.find(m_VertexFunc_z)->second << ")"
                << std::endl;
      std::cout << " Vertex mean (x,y,z) = (" << m_Vertex_x << "," << m_Vertex_y << "," << m_Vertex_z << ")" << std::endl;
      std::cout << " Vertex width (x,y,z) = (" << m_VertexWidth_x << "," << m_VertexWidth_y << "," << m_VertexWidth_z << ")" << std::endl;
    }
    std::cout << " Vertex size function (r) = ("
              << m_FunctionNames.find(m_VertexSizeFunc_r)->second << ")"
              << std::endl;
    std::cout << " Vertex size (mean) = (" << m_VertexSizeMean << ")" << std::endl;
    std::cout << " Vertex size (width) = (" << m_VertexSizeWidth << ")" << std::endl;
    if (std::isfinite(m_EtaMin) && std::isfinite(m_EtaMax))
    {
      std::cout << " Eta range = " << m_EtaMin << " - " << m_EtaMax << std::endl;
    }
    if (std::isfinite(m_ThetaMin) && std::isfinite(m_ThetaMax))
    {
      std::cout << " Theta range = " << m_ThetaMin << " - " << m_ThetaMax
                << ", deg: " << m_ThetaMin / M_PI * 180. << " - " << m_ThetaMax / M_PI * 180. << std::endl;
    }
    std::cout << " Phi range = " << m_PhiMin << " - " << m_PhiMax
              << ", deg: " << m_PhiMin / M_PI * 180. << " - " << m_PhiMax / M_PI * 180. << std::endl;
    if (std::isfinite(m_Pt_Min) && std::isfinite(m_Pt_Max))
    {
      std::cout << " pT range = " << m_Pt_Min << " - " << m_Pt_Max << std::endl;
    }
    if (std::isfinite(m_P_Min) && std::isfinite(m_P_Max))
    {
      std::cout << " p range = " << m_P_Min << " - " << m_P_Max << std::endl;
    }
    std::cout << " t0 = " << get_t0() << std::endl;
    std::cout << "===========================================================================" << std::endl;
  }

  // the definition table should be filled now, so convert codes into names
  for (unsigned int i = 0; i < _particle_codes.size(); ++i)
  {
    int pdgcode = _particle_codes[i].first;
    unsigned int count = _particle_codes[i].second;
    std::string pdgname = get_pdgname(pdgcode);
    _particle_names.push_back(std::make_pair(pdgname, count));
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int PHG4SimpleEventGenerator::process_event(PHCompositeNode *topNode)
{
  if (Verbosity() > 0)
  {
    std::cout << "====================== PHG4SimpleEventGenerator::process_event() =====================" << std::endl;
    std::cout << "PHG4SimpleEventGenerator::process_event - reuse_existing_vertex = " << get_reuse_existing_vertex() << std::endl;
  }

  if (!ReuseExistingVertex(topNode))
  {
    // generate a new vertex point
    set_vtx(smearvtx(m_Vertex_x, m_VertexWidth_x, m_VertexFunc_x),
            smearvtx(m_Vertex_y, m_VertexWidth_y, m_VertexFunc_y),
            smearvtx(m_Vertex_z, m_VertexWidth_z, m_VertexFunc_z));
  }
  set_vtx(get_vtx_x() + m_VertexOffset_x,
          get_vtx_y() + m_VertexOffset_y,
          get_vtx_z() + m_VertexOffset_z);

  if (Verbosity() > 0)
  {
    std::cout << "PHG4SimpleEventGenerator::process_event - vertex center" << get_reuse_existing_vertex()
              << get_vtx_x() << ", " << get_vtx_y() << ", " << get_vtx_z() << " cm"
              << std::endl;
  }

  int vtxindex = -1;
  int trackid = -1;
  for (unsigned int i = 0; i < _particle_names.size(); ++i)
  {
    std::string pdgname = _particle_names[i].first;
    int pdgcode = get_pdgcode(pdgname);
    unsigned int nparticles = _particle_names[i].second;

    for (unsigned int j = 0; j < nparticles; ++j)
    {
      if ((m_VertexSizeWidth > 0.0) || (m_VertexSizeMean != 0.0))
      {
        double r = smearvtx(m_VertexSizeMean, m_VertexSizeWidth, m_VertexSizeFunc_r);

        double x = 0.0;
        double y = 0.0;
        double z = 0.0;
        gsl_ran_dir_3d(RandomGenerator(), &x, &y, &z);
        x *= r;
        y *= r;
        z *= r;

        vtxindex = m_InEvent->AddVtx(get_vtx_x() + x, get_vtx_y() + y, get_vtx_z() + z, get_t0());
      }
      else if ((i == 0) && (j == 0))
      {
        vtxindex = m_InEvent->AddVtx(get_vtx_x(), get_vtx_y(), get_vtx_z(), get_t0());
      }

      ++trackid;

      double eta;
      if (!std::isnan(m_EtaMin) && !std::isnan(m_EtaMax))
      {
        eta = (m_EtaMax - m_EtaMin) * gsl_rng_uniform_pos(RandomGenerator()) + m_EtaMin;
      }
      else if (!std::isnan(m_ThetaMin) && !std::isnan(m_ThetaMax))
      {
        double theta = (m_ThetaMax - m_ThetaMin) * gsl_rng_uniform_pos(RandomGenerator()) + m_ThetaMin;
        eta = PHG4Utils::get_eta(theta);
      }
      else
      {
        std::cout << PHWHERE << "Error: neither eta range or theta range was specified" << std::endl;
        std::cout << "That should not happen, please inform the software group howthis happened" << std::endl;
        exit(-1);
      }

      double phi = (m_PhiMax - m_PhiMin) * gsl_rng_uniform_pos(RandomGenerator()) + m_PhiMin;      

      double pt;
      
      if (!std::isnan(m_P_Min) && !std::isnan(m_P_Max) && !std::isnan(m_P_GausWidth))
      {
	pt = ((m_P_Max - m_P_Min) * gsl_rng_uniform_pos(RandomGenerator()) + m_P_Min + gsl_ran_gaussian(RandomGenerator(), m_P_GausWidth)) / cosh(eta);
	if(!std::isnan(m_powerLawN))
	{
	  double y = gsl_rng_uniform_pos(RandomGenerator());
	  double x1 = pow(m_Pt_Max, m_powerLawN+1);
	  double x0 = pow(m_Pt_Min, m_powerLawN+1);
	  pt = pow((x1-x0)*y + x0,1./(m_powerLawN+1.));
	}
      }
      else if (!std::isnan(m_Pt_Min) && !std::isnan(m_Pt_Max) && !std::isnan(m_Pt_GausWidth))
      {
        pt = (m_Pt_Max - m_Pt_Min) * gsl_rng_uniform_pos(RandomGenerator()) + m_Pt_Min + gsl_ran_gaussian(RandomGenerator(), m_Pt_GausWidth);
	if(!std::isnan(m_powerLawN))
	{
	  double y = gsl_rng_uniform_pos(RandomGenerator());
	  double x1 = pow(m_Pt_Max, m_powerLawN+1);
	  double x0 = pow(m_Pt_Min, m_powerLawN+1);
	  pt = pow((x1-x0)*y + x0,1./(m_powerLawN+1.));
	}

      }
      else
      {
        std::cout << PHWHERE << "Error: neither a p range or pt range was specified" << std::endl;
        exit(-1);
      }

      double px = pt * cos(phi);
      double py = pt * sin(phi);
      double pz = pt * sinh(eta);
      double m = get_mass(pdgcode);
      double e = sqrt(px * px + py * py + pz * pz + m * m);

      PHG4Particle *particle = new PHG4Particlev2();
      particle->set_track_id(trackid);
      particle->set_vtx_id(vtxindex);
      particle->set_parent_id(0);
      particle->set_name(pdgname);
      particle->set_pid(pdgcode);
      particle->set_px(px);
      particle->set_py(py);
      particle->set_pz(pz);
      particle->set_e(e);

      m_InEvent->AddParticle(vtxindex, particle);
      if (EmbedFlag() != 0) m_InEvent->AddEmbeddedParticle(particle, EmbedFlag());
    }
  }

  if (Verbosity() > 0)
  {
    m_InEvent->identify();
    std::cout << "======================================================================================" << std::endl;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

double
PHG4SimpleEventGenerator::smearvtx(const double position, const double width, FUNCTION dist) const
{
  double res = position;
  if (dist == Uniform)
  {
    res = (position - width) + 2 * gsl_rng_uniform_pos(RandomGenerator()) * width;
  }
  else if (dist == Gaus)
  {
    res = position + gsl_ran_gaussian(RandomGenerator(), width);
  }
  else
  {
    std::cout << __PRETTY_FUNCTION__ << " invalid distribution function " << dist
              << " (" << m_FunctionNames.find(dist)->second << ")" << std::endl;
    gSystem->Exit(1);
  }
  return res;
}
