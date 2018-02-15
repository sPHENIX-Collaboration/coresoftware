#include "PHG4InnerHcalSubsystem.h"
#include "PHG4HcalDefs.h"
#include "PHG4InnerHcalDetector.h"
#include "PHG4InnerHcalSteppingAction.h"

#include <phparameter/PHParameters.h>

#include <g4main/PHG4HitContainer.h>

#include <pdbcalbase/PdbParameterMap.h>

#include <phool/getClass.h>

#include <Geant4/globals.hh>

#include <boost/foreach.hpp>

#include <set>
#include <sstream>

using namespace std;

//_______________________________________________________________________
PHG4InnerHcalSubsystem::PHG4InnerHcalSubsystem(const std::string &name, const int lyr)
  : PHG4DetectorSubsystem(name, lyr)
  , detector_(nullptr)
  , steppingAction_(nullptr)
{
  InitializeParameters();
}

//_______________________________________________________________________
int PHG4InnerHcalSubsystem::InitRunSubsystem(PHCompositeNode *topNode)
{
  PHNodeIterator iter(topNode);
  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));

  // create detector
  detector_ = new PHG4InnerHcalDetector(topNode, GetParams(), Name());
  detector_->SuperDetector(SuperDetector());
  detector_->OverlapCheck(CheckOverlap());
  set<string> nodes;
  if (GetParams()->get_int_param("active"))
  {
    PHNodeIterator dstIter(dstNode);
    PHCompositeNode *DetNode = dynamic_cast<PHCompositeNode *>(dstIter.findFirst("PHCompositeNode", SuperDetector()));
    if (!DetNode)
    {
      DetNode = new PHCompositeNode(SuperDetector());
      dstNode->addNode(DetNode);
    }

    ostringstream nodename;
    if (SuperDetector() != "NONE")
    {
      nodename << "G4HIT_" << SuperDetector();
    }
    else
    {
      nodename << "G4HIT_" << Name();
    }
    nodes.insert(nodename.str());
    if (GetParams()->get_int_param("absorberactive"))
    {
      nodename.str("");
      if (SuperDetector() != "NONE")
      {
        nodename << "G4HIT_ABSORBER_" << SuperDetector();
      }
      else
      {
        nodename << "G4HIT_ABSORBER_" << Name();
      }
      nodes.insert(nodename.str());
    }
    BOOST_FOREACH (string node, nodes)
    {
      PHG4HitContainer *g4_hits = findNode::getClass<PHG4HitContainer>(topNode, node.c_str());
      if (!g4_hits)
      {
        g4_hits = new PHG4HitContainer(node);
        DetNode->addNode(new PHIODataNode<PHObject>(g4_hits, node.c_str(), "PHObject"));
      }
    }

    // create stepping action
    steppingAction_ = new PHG4InnerHcalSteppingAction(detector_, GetParams());
  }
  else
  {
    // if this is a black hole it does not have to be active
    if (GetParams()->get_int_param("blackhole"))
    {
      steppingAction_ = new PHG4InnerHcalSteppingAction(detector_, GetParams());
    }
  }
  return 0;
}

//_______________________________________________________________________
int PHG4InnerHcalSubsystem::process_event(PHCompositeNode *topNode)
{
  // pass top node to stepping action so that it gets
  // relevant nodes needed internally
  if (steppingAction_)
  {
    steppingAction_->SetInterfacePointers(topNode);
  }
  return 0;
}

void PHG4InnerHcalSubsystem::Print(const string &what) const
{
  cout << Name() << " Parameters: " << endl;
  GetParams()->Print();
  if (detector_)
  {
    detector_->Print(what);
  }
  if (steppingAction_)
  {
    steppingAction_->Print(what);
  }

  return;
}

//_______________________________________________________________________
PHG4Detector *PHG4InnerHcalSubsystem::GetDetector(void) const
{
  return detector_;
}

void PHG4InnerHcalSubsystem::SetDefaultParameters()
{
  set_default_double_param(PHG4HcalDefs::innerrad, 117.27);
  set_default_double_param("light_balance_inner_corr", NAN);
  set_default_double_param("light_balance_inner_radius", NAN);
  set_default_double_param("light_balance_outer_corr", NAN);
  set_default_double_param("light_balance_outer_radius", NAN);
  set_default_double_param(PHG4HcalDefs::outerrad, 134.42);
  set_default_double_param("place_x", 0.);
  set_default_double_param("place_y", 0.);
  set_default_double_param("place_z", 0.);
  set_default_double_param("rot_x", 0.);
  set_default_double_param("rot_y", 0.);
  set_default_double_param("rot_z", 0.);
  set_default_double_param("scinti_eta_coverage", 1.1);
  set_default_double_param("scinti_gap_neighbor", 0.1);
  set_default_double_param("scinti_inner_gap", 0.85);
  set_default_double_param("scinti_outer_gap", 1.22 * (5.0 / 4.0));
// some math issue in the code subtracts 0.4mm so the scintillator
// does not end at 133.09 as per drawing but at 133.05
// adding 0.4mm compensates for this (so 133.13 gives the desired 133.09
  set_default_double_param("scinti_outer_radius", 133.13);
  set_default_double_param("scinti_tile_thickness", 0.7);
  set_default_double_param("size_z", 175.94 * 2);
  set_default_double_param("steplimits", NAN);
  set_default_double_param("tilt_angle", 36.15);  // engineering drawing
// corresponds very closely to 4 crossinge (35.5497 deg)

  set_default_int_param("light_scint_model", 1);
// if ncross is set (and tilt_angle is NAN) tilt_angle is calculated 
// from number of crossings
  set_default_int_param("ncross", 0); 
  set_default_int_param(PHG4HcalDefs::n_towers, 64);
  set_default_int_param(PHG4HcalDefs::scipertwr, 4);
  set_default_int_param(PHG4HcalDefs::n_scinti_tiles, 12);

  set_default_string_param("material", "SS310");
}

void PHG4InnerHcalSubsystem::SetLightCorrection(const double inner_radius, const double inner_corr, const double outer_radius, const double outer_corr)
{
  set_double_param("light_balance_inner_corr", inner_corr);
  set_double_param("light_balance_inner_radius", inner_radius);
  set_double_param("light_balance_outer_corr", outer_corr);
  set_double_param("light_balance_outer_radius", outer_radius);
  return;
}
