#include "PHG4TpcSubsystem.h"
#include "PHG4TpcDetector.h"
#include "PHG4TpcDisplayAction.h"
#include "PHG4TpcSteppingAction.h"

#include <phparameter/PHParameters.h>

#include <g4main/PHG4HitContainer.h>

#include <pdbcalbase/PdbParameterMap.h>

#include <phool/getClass.h>
#include <phool/PHCompositeNode.h>

#include <Geant4/globals.hh>

#include <boost/foreach.hpp>

#include <set>
#include <sstream>

using namespace std;

//_______________________________________________________________________
PHG4TpcSubsystem::PHG4TpcSubsystem(const std::string &name, const int lyr)
  : PHG4DetectorSubsystem(name, lyr)
  , detector_(nullptr)
  , steppingAction_(nullptr)
  , m_DisplayAction(nullptr)
{
  InitializeParameters();
}

//_______________________________________________________________________
PHG4TpcSubsystem::~PHG4TpcSubsystem()
{
  delete m_DisplayAction;
}

//_______________________________________________________________________
int PHG4TpcSubsystem::InitRunSubsystem(PHCompositeNode *topNode)
{
  PHNodeIterator iter(topNode);
  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));

  // create display settings before detector (detector adds its volumes to it)
  m_DisplayAction = new PHG4TpcDisplayAction(Name());
  // create detector
  detector_ = new PHG4TpcDetector(this, topNode, GetParams(), Name());
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
    steppingAction_ = new PHG4TpcSteppingAction(detector_, GetParams());
  }
  else
  {
    // if this is a black hole it does not have to be active
    if (GetParams()->get_int_param("blackhole"))
    {
      steppingAction_ = new PHG4TpcSteppingAction(detector_, GetParams());
    }
  }
  return 0;
}

//_______________________________________________________________________
int PHG4TpcSubsystem::process_event(PHCompositeNode *topNode)
{
  // pass top node to stepping action so that it gets
  // relevant nodes needed internally
  if (steppingAction_)
  {
    steppingAction_->SetInterfacePointers(topNode);
  }
  return 0;
}

void PHG4TpcSubsystem::Print(const string &what) const
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
PHG4Detector *PHG4TpcSubsystem::GetDetector(void) const
{
  return detector_;
}

void PHG4TpcSubsystem::SetDefaultParameters()
{
  set_default_double_param("gas_inner_radius", 21.);
  set_default_double_param("gas_outer_radius", 77.);
  set_default_double_param("place_x", 0.);
  set_default_double_param("place_y", 0.);
  set_default_double_param("place_z", 0.);
  set_default_double_param("rot_x", 0.);
  set_default_double_param("rot_y", 0.);
  set_default_double_param("rot_z", 0.);
  set_default_double_param("tpc_length", 211.);

  set_default_double_param("steplimits", NAN);

  set_default_string_param("tpc_gas", "sPHENIX_TPC_Gas");

  // material budget:
  // Cu (all layers): 0.5 oz cu per square foot, 1oz == 0.0347mm --> 0.5 oz ==  0.00347cm/2.
  // Kapton insulation 18 layers of * 5mil = 18*0.0127=0.2286
  // 250 um FR4 (Substrate for Cu layers)
  // HoneyComb (nomex) 1/2 inch=0.5*2.54 cm
  set_default_string_param("cage_layer_1_material", "G4_Cu");
  set_default_double_param("cage_layer_1_thickness", 0.00347 / 2.);

  set_default_string_param("cage_layer_2_material", "FR4");
  set_default_double_param("cage_layer_2_thickness", 0.025);

  set_default_string_param("cage_layer_3_material", "NOMEX");
  set_default_double_param("cage_layer_3_thickness", 0.5 * 2.54);

  set_default_string_param("cage_layer_4_material", "G4_Cu");
  set_default_double_param("cage_layer_4_thickness", 0.00347 / 2.);

  set_default_string_param("cage_layer_5_material", "FR4");
  set_default_double_param("cage_layer_5_thickness", 0.025);

  set_default_string_param("cage_layer_6_material", "G4_KAPTON");
  set_default_double_param("cage_layer_6_thickness", 0.2286);

  set_default_string_param("cage_layer_7_material", "G4_Cu");
  set_default_double_param("cage_layer_7_thickness", 0.00347 / 2.);

  set_default_string_param("cage_layer_8_material", "G4_KAPTON");
  set_default_double_param("cage_layer_8_thickness", 0.05);  // 50 um

  set_default_string_param("cage_layer_9_material", "G4_Cu");
  set_default_double_param("cage_layer_9_thickness", 0.00347 / 2.);

  // Thomas K Hemmick <Thomas.Hemmick@stonybrook.edu>
  //  The total thickness along Zed would be 5.6 millimeters (+/- 2.8 mm around Zed=0).
  //  The outer surfaces would have 0.005 inches (125 um) FR4 coated with a negligible thickness of Al. (revised to Au as below)
  //  The interior would be some stiffener of either honeycomb or rohacell.  The range of radiation lengths for this material are:
  //  Large cell honeycomb:  1450 cm  (0.028 g/cm^3 density)
  //  rohacell:  760 cm (0.052 g/cm^3 density)
  //  Close cell honeycomb:  635 cm (0.064 g/cm^3 density)
  //  I think a calculation just for the rohacell would be more than sufficient.
  set_default_string_param("window_core_material", "ROHACELL_FOAM_51");
  set_default_double_param("window_thickness", 0.56);  // overall thickness
  //I just checked with PC manufacturers and we can get 8.9 micron thick copper in reasonably large sheets.
  // At normal incidence, 8.9 microns is 0.06% of a radiation length.
  set_default_string_param("window_surface1_material", "G4_Cu");
  set_default_double_param("window_surface1_thickness", 8.9e-4);  // 8.9  um outter shell thickness be default
  // The FR4 should be either 5 or 10 mils thick.  10 mils is 254 microns and 5 mils is 0.127 microns.  I think either of these is mechanically fine...
  set_default_string_param("window_surface2_material", "FR4");
  set_default_double_param("window_surface2_thickness", 0.0127);  // 127  um 2nd shell thickness be default
}
