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

#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHDataNode.h>
#include <phool/PHRandomSeed.h>
#include <phool/getClass.h>

#include <HepMC/GenEvent.h>
#include <HepMC/IO_GenEvent.h>

#include <gsl/gsl_const.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>

#include <cassert>
#include <climits>
#include <iostream>

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
  gsl_rng_set(RandomGenerator, seed );
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

//! send HepMC::GenEvent to DST tree. This function takes ownership of evt
PHHepMCGenEvent *PHHepMCGenHelper::insert_event(HepMC::GenEvent *evt)
{
  assert(evt);
  assert(_geneventmap);

  PHHepMCGenEvent *genevent = _geneventmap->insert_event(_embedding_id);

  genevent->addEvent(evt);
  move_vertex(genevent);

  return genevent;
}

void PHHepMCGenHelper::move_vertex(PHHepMCGenEvent *genevent)
{
  assert(genevent);

  assert(_vertex_width_x >= 0);

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

  genevent->moveVertex(
      (smear(_vertex_x, _vertex_width_x, _vertex_func_x)),
      (smear(_vertex_y, _vertex_width_y, _vertex_func_y)),
      (smear(_vertex_z, _vertex_width_z, _vertex_func_z)),
      (smear(_vertex_t, _vertex_width_t, _vertex_func_t)));
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
