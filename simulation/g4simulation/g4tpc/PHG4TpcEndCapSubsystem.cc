//____________________________________________________________________________..
//
// This is the interface to the framework. You only need to define the parameters
// you use for your detector in the SetDefaultParameters() method here
// The place to do this is marked by //implement your own here//
// The parameters have no units, they need to be converted in the
// PHG4TpcEndCapDetector::ConstructMe() method
// but the convention is as mentioned cm and deg
//____________________________________________________________________________..
//
#include "PHG4TpcEndCapSubsystem.h"

#include "PHG4TpcEndCapDetector.h"
#include "PHG4TpcEndCapDisplayAction.h"
#include "PHG4TpcEndCapSteppingAction.h"

#include <phparameter/PHParameters.h>

#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4SteppingAction.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>
#include <phool/getClass.h>

using namespace std;

//_______________________________________________________________________
PHG4TpcEndCapSubsystem::PHG4TpcEndCapSubsystem(const std::string &name)
  : PHG4DetectorSubsystem(name)
  , m_Detector(nullptr)
  , m_SteppingAction(nullptr)
  , m_DisplayAction(nullptr)
{
  // call base class method which will set up parameter infrastructure
  // and call our SetDefaultParameters() method
  InitializeParameters();
}

PHG4TpcEndCapSubsystem::~PHG4TpcEndCapSubsystem()
{
  if (m_DisplayAction)
    delete m_DisplayAction;
}

//_______________________________________________________________________
int PHG4TpcEndCapSubsystem::InitRunSubsystem(PHCompositeNode *topNode)
{
  PHNodeIterator iter(topNode);
  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
  PHNodeIterator dstIter(dstNode);

  // create display settings before detector (detector adds its volumes to it)
  m_DisplayAction = new PHG4TpcEndCapDisplayAction(Name());

  if (GetParams()->get_int_param("active"))
  {
    PHCompositeNode *DetNode = dynamic_cast<PHCompositeNode *>(dstIter.findFirst("PHCompositeNode", Name()));
    if (!DetNode)
    {
      DetNode = new PHCompositeNode(Name());
      dstNode->addNode(DetNode);
    }
    string g4hitnodename = "G4HIT_" + Name();
    PHG4HitContainer *g4_hits = findNode::getClass<PHG4HitContainer>(DetNode, g4hitnodename);
    if (!g4_hits)
    {
      g4_hits = new PHG4HitContainer(g4hitnodename);
      DetNode->addNode(new PHIODataNode<PHObject>(g4_hits, g4hitnodename, "PHObject"));
    }
  }
  // create detector
  m_Detector = new PHG4TpcEndCapDetector(this, topNode, GetParams(), Name());
  m_Detector->OverlapCheck(CheckOverlap());
  // create stepping action if detector is active
  if (GetParams()->get_int_param("active"))
  {
    m_SteppingAction = new PHG4TpcEndCapSteppingAction(m_Detector, GetParams());
  }
  return 0;
}
//_______________________________________________________________________
int PHG4TpcEndCapSubsystem::process_event(PHCompositeNode *topNode)
{
  // pass top node to stepping action so that it gets
  // relevant nodes needed internally
  if (m_SteppingAction)
  {
    m_SteppingAction->SetInterfacePointers(topNode);
  }
  return 0;
}
//_______________________________________________________________________
void PHG4TpcEndCapSubsystem::Print(const string &what) const
{
  if (m_Detector)
  {
    m_Detector->Print(what);
  }
  return;
}

//_______________________________________________________________________
PHG4Detector *PHG4TpcEndCapSubsystem::GetDetector(void) const
{
  return m_Detector;
}

//_______________________________________________________________________
void PHG4TpcEndCapSubsystem::SetDefaultParameters()
{
  set_default_int_param("construction_verbosity", 0);
  // sizes are in cm
  // angles are in deg
  // units should be converted to G4 units when used
  //implement your own here//
  set_default_double_param("place_x", 0.);
  set_default_double_param("place_y", 0.);
  set_default_double_param("place_z", 0.);
  set_default_double_param("rot_x", 0.);
  set_default_double_param("rot_y", 0.);
  set_default_double_param("rot_z", 0.);

  set_default_double_param("envelop_r_min", 20.);
  set_default_double_param("envelop_r_max", 79.);
  set_default_double_param("envelop_front_surface_z", 211. / 2.);

  set_default_int_param("n_GEM_layers", 4);

  const double inch_to_cm = 2.54;

  set_default_int_param("n_sectors", 12);
  set_default_int_param("n_radial_modules", 3);

  set_default_string_param("wagon_wheel_material", "G4_Al");

  set_default_double_param("wagon_wheel_sector_phi_offset_degree", 0);

  set_default_double_param("wagon_wheel_front_frame_thickness", inch_to_cm * .38);
  set_default_double_param("wagon_wheel_front_frame_spoke_width", inch_to_cm * 1.12);

  set_default_double_param("wagon_wheel_front_frame_R_inner", inch_to_cm * 15.74 / 2.);
  set_default_double_param("wagon_wheel_front_frame_R_outer", inch_to_cm * 30.7);

  set_default_double_param("wagon_wheel_front_frame_R_R1_inner", inch_to_cm * 9.81);
  set_default_double_param("wagon_wheel_front_frame_R_R1_outer", inch_to_cm * 15.47);

  set_default_double_param("wagon_wheel_front_frame_R_R2_inner", inch_to_cm * 16.59);
  set_default_double_param("wagon_wheel_front_frame_R_R2_outer", inch_to_cm * 22.24);

  set_default_double_param("wagon_wheel_front_frame_R_R3_inner", inch_to_cm * 23.36);
  set_default_double_param("wagon_wheel_front_frame_R_R3_outer", inch_to_cm * 29.02);

  set_default_double_param("wagon_wheel_rim_outer_Rin", inch_to_cm * 29.58);
  set_default_double_param("wagon_wheel_rim_outer_Rout", inch_to_cm * 60.49 / 2.);
  set_default_double_param("wagon_wheel_rim_outer_thickness", inch_to_cm * (4.5 - .38));

  set_default_double_param("wagon_wheel_spoke_width", inch_to_cm * .36);
  set_default_double_param("wagon_wheel_spoke_height_inner", inch_to_cm * (1.715 - .38));
  set_default_double_param("wagon_wheel_spoke_height_outer", inch_to_cm * (4.5 - .38));
  set_default_double_param("wagon_wheel_spoke_R_inner", inch_to_cm * 9);
  set_default_double_param("wagon_wheel_spoke_R_outer", inch_to_cm * 29.5);

  set_default_int_param("electronics_enable", 1);
  set_default_int_param("electronics_nFEE_R1", 6);
  set_default_int_param("electronics_nFEE_R2", 8);
  set_default_int_param("electronics_nFEE_R3", 12);

  set_default_double_param("electronics_FEE_depth", inch_to_cm * 15.74 / 2.);
  set_default_double_param("electronics_FEE_Cu_thickness", 35e-4);
  set_default_double_param("electronics_FEE_PCB_thickness", 0.2);
  set_default_double_param("electronics_FEE_Al_thickness", 0.2);

  set_default_string_param("electronics_cooling_block_material", "G4_Al");
  set_default_double_param("electronics_cooling_block_thickness", inch_to_cm * 3.5);

  set_default_double_param("electronics_cooling_block_R_inner", inch_to_cm * 9.26);
  set_default_double_param("electronics_cooling_block_R_outer", inch_to_cm * 29.57);

  set_default_double_param("electronics_cooling_block_R_R1_inner", inch_to_cm * 9.81);
  set_default_double_param("electronics_cooling_block_R_R1_outer", inch_to_cm * 15.47);

  set_default_double_param("electronics_cooling_block_R_R2_inner", inch_to_cm * 16.59);
  set_default_double_param("electronics_cooling_block_R_R2_outer", inch_to_cm * 22.24);

  set_default_double_param("electronics_cooling_block_R_R3_inner", inch_to_cm * 23.36);
  set_default_double_param("electronics_cooling_block_R_R3_outer", inch_to_cm * 29.02);
}
