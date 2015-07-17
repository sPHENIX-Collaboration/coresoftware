// $$Id: PHG4SpacalDetector.cc,v 1.7 2015/02/10 15:39:26 pinkenbu Exp $$

/*!
 * \file ${file_name}
 * \brief
 * \author Jin Huang <jhuang@bnl.gov>
 * \version $$Revision: 1.7 $$
 * \date $$Date: 2015/02/10 15:39:26 $$
 */
#include "PHG4SpacalDetector.h"
#include "PHG4CylinderGeomContainer.h"
#include "PHG4CylinderGeom_Spacalv1.h"

#include <g4main/PHG4PhenixDetector.h>
#include <g4main/PHG4Utils.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <fun4all/getClass.h>

#include <Geant4/G4Box.hh>
#include <Geant4/G4Colour.hh>
#include <Geant4/G4Cons.hh>
#include <Geant4/G4Trap.hh>
#include <Geant4/G4ExtrudedSolid.hh>
#include <Geant4/G4LogicalVolume.hh>
#include <Geant4/G4UserLimits.hh>
#include <Geant4/G4Material.hh>
#include <Geant4/G4PhysicalConstants.hh>
#include <Geant4/G4PVPlacement.hh>
#include <Geant4/G4SubtractionSolid.hh>
#include <Geant4/G4SystemOfUnits.hh>
#include <Geant4/G4Tubs.hh>
#include <Geant4/G4TwoVector.hh>
#include <Geant4/G4UnionSolid.hh>
#include <Geant4/G4VisAttributes.hh>

#include <cassert>
#include <cmath>
#include <sstream>

using namespace std;

//_______________________________________________________________
//note this inactive thickness is ~1.5% of a radiation length
PHG4SpacalDetector::PHG4SpacalDetector(PHCompositeNode *Node,
    const std::string &dnam, SpacalGeom_t * geom, const int lyr) :
    PHG4Detector(Node, dnam), _region(NULL), cylinder_solid(NULL), cylinder_logic(
        NULL), cylinder_physi(NULL), active(0), absorberactive(0), layer(
        lyr), _geom(geom)
{

  if (_geom == NULL)
    {
      cout <<"PHG4SpacalDetector::Constructor - Fatal Error - invalid geometry object!"<<endl;
      exit(1);
    }

  step_limits = new G4UserLimits(_geom->get_calo_step_size() * cm);

  clading_step_limits = new G4UserLimits(
      _geom->get_fiber_clading_step_size() * cm);

  fiber_core_step_limits = new G4UserLimits(
      _geom->get_fiber_core_step_size() * cm);
}

PHG4SpacalDetector::~PHG4SpacalDetector(void)
{
  // deleting NULL pointers is allowed (results in NOOP) 
  // so checking for not null before deleting is not needed
    delete step_limits;
    delete clading_step_limits;
    delete fiber_core_step_limits;
}

//_______________________________________________________________
int
PHG4SpacalDetector::IsInCylinderActive(const G4VPhysicalVolume * volume)
{
  //  cout << "checking detector" << endl;
  if (active && fiber_core_vol.find(volume) != fiber_core_vol.end())
    {
//      return fiber_core_vol.find(volume)->second;
      return FIBER_CORE;
    }
  if (absorberactive)
    {
      if (fiber_vol.find(volume) != fiber_vol.end())
        return FIBER_CLADING;

      if (calo_vol.find(volume) != calo_vol.end())
        return ABSORBER;

      if (block_vol.find(volume) != block_vol.end())
        return ABSORBER;

    }
  return INACTIVE;
}

//_______________________________________________________________
void
PHG4SpacalDetector::Construct(G4LogicalVolume* logicWorld)
{
  assert(_geom);

  if ((_geom->get_zmin() * cm + _geom->get_zmax() * cm) / 2
      != _geom->get_zpos() * cm)
    {
      cout
          << "PHG4SpacalDetector::Construct - ERROR - not yet support unsymmetric system. Let me know if you need it. - Jin"
          << endl;
      _geom->Print();
      exit(-1);
    }
  if (_geom->get_zmin() * cm >= _geom->get_zmax() * cm)
    {
      cout << "PHG4SpacalDetector::Construct - ERROR - zmin >= zmax!" << endl;
      _geom->Print();
      exit(-1);
    }

//  G4Tubs* _cylinder_solid = new G4Tubs(G4String(GetName().c_str()),
//      _geom->get_radius() * cm, _geom->get_max_radius() * cm,
//      _geom->get_length() * cm / 2.0, 0, twopi);
  G4Trap* _cylinder_solid;
  G4Transform3D block_trans ;

  {
    //Processed PostionSeeds 1 from 1 1
    _cylinder_solid = new G4Trap( /*const G4String& pName*/G4String(GetName().c_str()),
     6.751951 * cm,// G4double pDz,
      0,0,// G4double pTheta, G4double pPhi,
      2.247390*cm, 9.653108*cm , 9.642601*cm, // G4double pDy1, G4double pDx1, G4double pDx2,
       0,// G4double pAlp1,
      2.567902*cm, 10.983502*cm , 10.971494*cm,// G4double pDy2, G4double pDx3, G4double pDx4,
      0// G4double pAlp2 //
        );
     block_trans = G4Translate3D(0.000000*cm, 105.088944*cm , 2.483168*cm) * G4RotateX3D(-1.547066*rad);

//  cylinder_solid = _cylinder_solid;

      G4Material * cylinder_mat = G4Material::GetMaterial(_geom->get_absorber_mat());
      assert(cylinder_mat);

      cylinder_logic = new G4LogicalVolume(_cylinder_solid, cylinder_mat,
          G4String(GetName().c_str()), 0, 0, 0);
      G4VisAttributes* VisAtt = new G4VisAttributes();
      PHG4Utils::SetColour(VisAtt, "W_Epoxy");
      VisAtt->SetVisibility(true);
      VisAtt->SetForceSolid(
          (not _geom->is_virualize_fiber())
              and (not _geom->is_azimuthal_seg_visible()));
      cylinder_logic->SetVisAttributes(VisAtt);

      cylinder_physi = new G4PVPlacement(block_trans, cylinder_logic,
          G4String(GetName().c_str()), logicWorld, false, 0, overlapcheck);

    }


  {
    //Processed PostionSeeds 3 from 3 1
    _cylinder_solid = new G4Trap( /*const G4String& pName*/G4String(GetName().c_str()),
     6.751976 * cm,// G4double pDz,
      0,0,// G4double pTheta, G4double pPhi,
      2.248520*cm, 9.607852*cm , 9.658661*cm, // G4double pDy1, G4double pDx1, G4double pDx2,
       0,// G4double pAlp1,
      2.566732*cm, 10.925488*cm , 10.983502*cm,// G4double pDy2, G4double pDx3, G4double pDx4,
      0// G4double pAlp2 //
        );
     block_trans = G4Translate3D(0.000000*cm, 104.898063*cm , 12.114086*cm) * G4RotateX3D(-1.455797*rad);

      G4Material * cylinder_mat = G4Material::GetMaterial(_geom->get_absorber_mat());
      assert(cylinder_mat);

      cylinder_logic = new G4LogicalVolume(_cylinder_solid, cylinder_mat,
          G4String(GetName().c_str()), 0, 0, 0);
      G4VisAttributes* VisAtt = new G4VisAttributes();
      PHG4Utils::SetColour(VisAtt, "W_Epoxy");
      VisAtt->SetVisibility(true);
      VisAtt->SetForceSolid(
          (not _geom->is_virualize_fiber())
              and (not _geom->is_azimuthal_seg_visible()));
      cylinder_logic->SetVisAttributes(VisAtt);

      cylinder_physi = new G4PVPlacement(block_trans, cylinder_logic,
          G4String(GetName().c_str()), logicWorld, false, 0, overlapcheck);

    }

  {

    //Processed PostionSeeds 5 from 5 1
    _cylinder_solid = new G4Trap( /*const G4String& pName*/G4String(GetName().c_str()),
     6.751971 * cm,// G4double pDz,
      0,0,// G4double pTheta, G4double pPhi,
      2.249652*cm, 9.584687*cm , 9.675123*cm, // G4double pDy1, G4double pDx1, G4double pDx2,
       0,// G4double pAlp1,
      2.565431*cm, 10.880358*cm , 10.983502*cm,// G4double pDy2, G4double pDx3, G4double pDx4,
      0// G4double pAlp2 //
        );
     block_trans = G4Translate3D(0.000000*cm, 104.766491*cm , 21.838412*cm) * G4RotateX3D(-1.365252*rad);

    G4Material * cylinder_mat = G4Material::GetMaterial(_geom->get_absorber_mat());
    assert(cylinder_mat);

    cylinder_logic = new G4LogicalVolume(_cylinder_solid, cylinder_mat,
        G4String(GetName().c_str()), 0, 0, 0);
    G4VisAttributes* VisAtt = new G4VisAttributes();
    PHG4Utils::SetColour(VisAtt, "W_Epoxy");
    VisAtt->SetVisibility(true);
    VisAtt->SetForceSolid(
        (not _geom->is_virualize_fiber())
            and (not _geom->is_azimuthal_seg_visible()));
    cylinder_logic->SetVisAttributes(VisAtt);

    cylinder_physi = new G4PVPlacement(block_trans, cylinder_logic,
        G4String(GetName().c_str()), logicWorld, false, 0, overlapcheck);

  }
  {

    //Processed PostionSeeds 7 from 7 1
    _cylinder_solid = new G4Trap( /*const G4String& pName*/G4String(GetName().c_str()),
     6.751918 * cm,// G4double pDz,
      0,0,// G4double pTheta, G4double pPhi,
      2.251626*cm, 9.573204*cm , 9.701866*cm, // G4double pDy1, G4double pDx1, G4double pDx2,
       0,// G4double pAlp1,
      2.560214*cm, 10.837179*cm , 10.983502*cm,// G4double pDy2, G4double pDx3, G4double pDx4,
      0// G4double pAlp2 //
        );
     block_trans = G4Translate3D(0.000000*cm, 104.695625*cm , 31.745365*cm) * G4RotateX3D(-1.276437*rad);

    G4Material * cylinder_mat = G4Material::GetMaterial(_geom->get_absorber_mat());
    assert(cylinder_mat);

    cylinder_logic = new G4LogicalVolume(_cylinder_solid, cylinder_mat,
        G4String(GetName().c_str()), 0, 0, 0);
    G4VisAttributes* VisAtt = new G4VisAttributes();
    PHG4Utils::SetColour(VisAtt, "W_Epoxy");
    VisAtt->SetVisibility(true);
    VisAtt->SetForceSolid(
        (not _geom->is_virualize_fiber())
            and (not _geom->is_azimuthal_seg_visible()));
    cylinder_logic->SetVisAttributes(VisAtt);

    cylinder_physi = new G4PVPlacement(block_trans, cylinder_logic,
        G4String(GetName().c_str()), logicWorld, false, 0, overlapcheck);

  }
  {

    //Processed PostionSeeds 9 from 9 1
    _cylinder_solid = new G4Trap( /*const G4String& pName*/G4String(GetName().c_str()),
     6.751878 * cm,// G4double pDz,
      0,0,// G4double pTheta, G4double pPhi,
      2.253432*cm, 9.572854*cm , 9.737865*cm, // G4double pDy1, G4double pDx1, G4double pDx2,
       0,// G4double pAlp1,
      2.552634*cm, 10.796576*cm , 10.983502*cm,// G4double pDy2, G4double pDx3, G4double pDx4,
      0// G4double pAlp2 //
        );
     block_trans = G4Translate3D(0.000000*cm, 104.683052*cm , 41.925557*cm) * G4RotateX3D(-1.189890*rad);

    G4Material * cylinder_mat = G4Material::GetMaterial(_geom->get_absorber_mat());
    assert(cylinder_mat);

    cylinder_logic = new G4LogicalVolume(_cylinder_solid, cylinder_mat,
        G4String(GetName().c_str()), 0, 0, 0);
    G4VisAttributes* VisAtt = new G4VisAttributes();
    PHG4Utils::SetColour(VisAtt, "W_Epoxy");
    VisAtt->SetVisibility(true);
    VisAtt->SetForceSolid(
        (not _geom->is_virualize_fiber())
            and (not _geom->is_azimuthal_seg_visible()));
    cylinder_logic->SetVisAttributes(VisAtt);

    cylinder_physi = new G4PVPlacement(block_trans, cylinder_logic,
        G4String(GetName().c_str()), logicWorld, false, 0, overlapcheck);

  }
  {

    //Processed PostionSeeds 11 from 11 1
    _cylinder_solid = new G4Trap( /*const G4String& pName*/G4String(GetName().c_str()),
     6.751840 * cm,// G4double pDz,
      0,0,// G4double pTheta, G4double pPhi,
      2.256511*cm, 9.582510*cm , 9.781544*cm, // G4double pDy1, G4double pDx1, G4double pDx2,
       0,// G4double pAlp1,
      2.546341*cm, 10.758926*cm , 10.983502*cm,// G4double pDy2, G4double pDx3, G4double pDx4,
      0// G4double pAlp2 //
        );
     block_trans = G4Translate3D(0.000000*cm, 104.722867*cm , 52.466304*cm) * G4RotateX3D(-1.106480*rad);

    G4Material * cylinder_mat = G4Material::GetMaterial(_geom->get_absorber_mat());
    assert(cylinder_mat);

    cylinder_logic = new G4LogicalVolume(_cylinder_solid, cylinder_mat,
        G4String(GetName().c_str()), 0, 0, 0);
    G4VisAttributes* VisAtt = new G4VisAttributes();
    PHG4Utils::SetColour(VisAtt, "W_Epoxy");
    VisAtt->SetVisibility(true);
    VisAtt->SetForceSolid(
        (not _geom->is_virualize_fiber())
            and (not _geom->is_azimuthal_seg_visible()));
    cylinder_logic->SetVisAttributes(VisAtt);

    cylinder_physi = new G4PVPlacement(block_trans, cylinder_logic,
        G4String(GetName().c_str()), logicWorld, false, 0, overlapcheck);

  }
  {


    //Processed PostionSeeds 13 from 13 1
    _cylinder_solid = new G4Trap( /*const G4String& pName*/G4String(GetName().c_str()),
     6.751779 * cm,// G4double pDz,
      0,0,// G4double pTheta, G4double pPhi,
      2.259203*cm, 9.601148*cm , 9.831603*cm, // G4double pDy1, G4double pDx1, G4double pDx2,
       0,// G4double pAlp1,
      2.537383*cm, 10.724678*cm , 10.983502*cm,// G4double pDy2, G4double pDx3, G4double pDx4,
      0// G4double pAlp2 //
        );
     block_trans = G4Translate3D(0.000000*cm, 104.810306*cm , 63.468504*cm) * G4RotateX3D(-1.026390*rad);
    G4Material * cylinder_mat = G4Material::GetMaterial(_geom->get_absorber_mat());
    assert(cylinder_mat);

    cylinder_logic = new G4LogicalVolume(_cylinder_solid, cylinder_mat,
        G4String(GetName().c_str()), 0, 0, 0);
    G4VisAttributes* VisAtt = new G4VisAttributes();
    PHG4Utils::SetColour(VisAtt, "W_Epoxy");
    VisAtt->SetVisibility(true);
    VisAtt->SetForceSolid(
        (not _geom->is_virualize_fiber())
            and (not _geom->is_azimuthal_seg_visible()));
    cylinder_logic->SetVisAttributes(VisAtt);

    cylinder_physi = new G4PVPlacement(block_trans, cylinder_logic,
        G4String(GetName().c_str()), logicWorld, false, 0, overlapcheck);

  }
  {

    //Processed PostionSeeds 15 from 15 1
    _cylinder_solid = new G4Trap( /*const G4String& pName*/G4String(GetName().c_str()),
     6.827956 * cm,// G4double pDz,
      0,0,// G4double pTheta, G4double pPhi,
      2.260314*cm, 9.615107*cm , 9.873956*cm, // G4double pDy1, G4double pDx1, G4double pDx2,
       0,// G4double pAlp1,
      2.527163*cm, 10.694108*cm , 10.983502*cm,// G4double pDy2, G4double pDx3, G4double pDx4,
      0// G4double pAlp2 //
        );
     block_trans = G4Translate3D(0.000000*cm, 104.875648*cm , 74.978133*cm) * G4RotateX3D(-0.950304*rad);

    G4Material * cylinder_mat = G4Material::GetMaterial(_geom->get_absorber_mat());
    assert(cylinder_mat);

    cylinder_logic = new G4LogicalVolume(_cylinder_solid, cylinder_mat,
        G4String(GetName().c_str()), 0, 0, 0);
    G4VisAttributes* VisAtt = new G4VisAttributes();
    PHG4Utils::SetColour(VisAtt, "W_Epoxy");
    VisAtt->SetVisibility(true);
    VisAtt->SetForceSolid(
        (not _geom->is_virualize_fiber())
            and (not _geom->is_azimuthal_seg_visible()));
    cylinder_logic->SetVisAttributes(VisAtt);

    cylinder_physi = new G4PVPlacement(block_trans, cylinder_logic,
        G4String(GetName().c_str()), logicWorld, false, 0, overlapcheck);

  }
  {

    //Processed PostionSeeds 17 from 17 1
    _cylinder_solid = new G4Trap( /*const G4String& pName*/G4String(GetName().c_str()),
     7.050129 * cm,// G4double pDz,
      0,0,// G4double pTheta, G4double pPhi,
      2.256143*cm, 9.614732*cm , 9.898473*cm, // G4double pDy1, G4double pDx1, G4double pDx2,
       0,// G4double pAlp1,
      2.514461*cm, 10.667289*cm , 10.983502*cm,// G4double pDy2, G4double pDx3, G4double pDx4,
      0// G4double pAlp2 //
        );
     block_trans = G4Translate3D(0.000000*cm, 104.868853*cm , 87.029925*cm) * G4RotateX3D(-0.878224*rad);

    G4Material * cylinder_mat = G4Material::GetMaterial(_geom->get_absorber_mat());
    assert(cylinder_mat);

    cylinder_logic = new G4LogicalVolume(_cylinder_solid, cylinder_mat,
        G4String(GetName().c_str()), 0, 0, 0);
    G4VisAttributes* VisAtt = new G4VisAttributes();
    PHG4Utils::SetColour(VisAtt, "W_Epoxy");
    VisAtt->SetVisibility(true);
    VisAtt->SetForceSolid(
        (not _geom->is_virualize_fiber())
            and (not _geom->is_azimuthal_seg_visible()));
    cylinder_logic->SetVisAttributes(VisAtt);

    cylinder_physi = new G4PVPlacement(block_trans, cylinder_logic,
        G4String(GetName().c_str()), logicWorld, false, 0, overlapcheck);

  }
  {


    //Processed PostionSeeds 19 from 19 1
    _cylinder_solid = new G4Trap( /*const G4String& pName*/G4String(GetName().c_str()),
     7.329495 * cm,// G4double pDz,
      0,0,// G4double pTheta, G4double pPhi,
      2.251071*cm, 9.614857*cm , 9.920437*cm, // G4double pDy1, G4double pDx1, G4double pDx2,
       0,// G4double pAlp1,
      2.501607*cm, 10.643899*cm , 10.983502*cm,// G4double pDy2, G4double pDx3, G4double pDx4,
      0// G4double pAlp2 //
        );
     block_trans = G4Translate3D(0.000000*cm, 104.865551*cm , 99.766946*cm) * G4RotateX3D(-0.810340*rad);
    G4Material * cylinder_mat = G4Material::GetMaterial(_geom->get_absorber_mat());
    assert(cylinder_mat);

    cylinder_logic = new G4LogicalVolume(_cylinder_solid, cylinder_mat,
        G4String(GetName().c_str()), 0, 0, 0);
    G4VisAttributes* VisAtt = new G4VisAttributes();
    PHG4Utils::SetColour(VisAtt, "W_Epoxy");
    VisAtt->SetVisibility(true);
    VisAtt->SetForceSolid(
        (not _geom->is_virualize_fiber())
            and (not _geom->is_azimuthal_seg_visible()));
    cylinder_logic->SetVisAttributes(VisAtt);

    cylinder_physi = new G4PVPlacement(block_trans, cylinder_logic,
        G4String(GetName().c_str()), logicWorld, false, 0, overlapcheck);

  }
  {

    //Processed PostionSeeds 21 from 21 1
    _cylinder_solid = new G4Trap( /*const G4String& pName*/G4String(GetName().c_str()),
     7.685134 * cm,// G4double pDz,
      0,0,// G4double pTheta, G4double pPhi,
      2.247516*cm, 9.612530*cm , 9.937324*cm, // G4double pDy1, G4double pDx1, G4double pDx2,
       0,// G4double pAlp1,
      2.491419*cm, 10.623460*cm , 10.983502*cm,// G4double pDy2, G4double pDx3, G4double pDx4,
      0// G4double pAlp2 //
        );
     block_trans = G4Translate3D(0.000000*cm, 104.850629*cm , 113.268570*cm) * G4RotateX3D(-0.746991*rad);

    G4Material * cylinder_mat = G4Material::GetMaterial(_geom->get_absorber_mat());
    assert(cylinder_mat);

    cylinder_logic = new G4LogicalVolume(_cylinder_solid, cylinder_mat,
        G4String(GetName().c_str()), 0, 0, 0);
    G4VisAttributes* VisAtt = new G4VisAttributes();
    PHG4Utils::SetColour(VisAtt, "W_Epoxy");
    VisAtt->SetVisibility(true);
    VisAtt->SetForceSolid(
        (not _geom->is_virualize_fiber())
            and (not _geom->is_azimuthal_seg_visible()));
    cylinder_logic->SetVisAttributes(VisAtt);

    cylinder_physi = new G4PVPlacement(block_trans, cylinder_logic,
        G4String(GetName().c_str()), logicWorld, false, 0, overlapcheck);

  }
  {

    //Processed PostionSeeds 23 from 23 1
    _cylinder_solid = new G4Trap( /*const G4String& pName*/G4String(GetName().c_str()),
     8.104392 * cm,// G4double pDz,
      0,0,// G4double pTheta, G4double pPhi,
      2.238046*cm, 9.611255*cm , 9.951909*cm, // G4double pDy1, G4double pDx1, G4double pDx2,
       0,// G4double pAlp1,
      2.481231*cm, 10.605823*cm , 10.983502*cm,// G4double pDy2, G4double pDx3, G4double pDx4,
      0// G4double pAlp2 //
        );
    block_trans = G4Translate3D(0.000000*cm, 104.839643*cm , 127.643953*cm) * G4RotateX3D(-0.687482*rad);

    G4Material * cylinder_mat = G4Material::GetMaterial(_geom->get_absorber_mat());
    assert(cylinder_mat);

    cylinder_logic = new G4LogicalVolume(_cylinder_solid, cylinder_mat,
        G4String(GetName().c_str()), 0, 0, 0);
    G4VisAttributes* VisAtt = new G4VisAttributes();
    PHG4Utils::SetColour(VisAtt, "W_Epoxy");
    VisAtt->SetVisibility(true);
    VisAtt->SetForceSolid(
        (not _geom->is_virualize_fiber())
            and (not _geom->is_azimuthal_seg_visible()));
    cylinder_logic->SetVisAttributes(VisAtt);

    cylinder_physi = new G4PVPlacement(block_trans, cylinder_logic,
        G4String(GetName().c_str()), logicWorld, false, 0, overlapcheck);

  }


  {

    //Processed PostionSeeds 24 from 24 1
    _cylinder_solid = new G4Trap( /*const G4String& pName*/G4String(GetName().c_str()),
     8.326711 * cm,// G4double pDz,
      0,0,// G4double pTheta, G4double pPhi,
      2.238994*cm, 9.612956*cm , 9.961690*cm, // G4double pDy1, G4double pDx1, G4double pDx2,
       0,// G4double pAlp1,
      2.468562*cm, 10.599019*cm , 10.983502*cm,// G4double pDy2, G4double pDx3, G4double pDx4,
      0// G4double pAlp2 //
        );
     block_trans = G4Translate3D(0.000000*cm, 104.851518*cm , 135.224711*cm) * G4RotateX3D(-0.658682*rad);

    G4Material * cylinder_mat = G4Material::GetMaterial(_geom->get_absorber_mat());
    assert(cylinder_mat);

    cylinder_logic = new G4LogicalVolume(_cylinder_solid, cylinder_mat,
        G4String(GetName().c_str()), 0, 0, 0);
    G4VisAttributes* VisAtt = new G4VisAttributes();
    PHG4Utils::SetColour(VisAtt, "W_Epoxy");
    VisAtt->SetVisibility(true);
    VisAtt->SetForceSolid(
        (not _geom->is_virualize_fiber())
            and (not _geom->is_azimuthal_seg_visible()));
    cylinder_logic->SetVisAttributes(VisAtt);

    cylinder_physi = new G4PVPlacement(block_trans, cylinder_logic,
        G4String(GetName().c_str()), logicWorld, false, 0, overlapcheck);

  }

  {

    //Processed PostionSeeds 22 from 22 1
    _cylinder_solid = new G4Trap( /*const G4String& pName*/G4String(GetName().c_str()),
     7.875754 * cm,// G4double pDz,
      0,0,// G4double pTheta, G4double pPhi,
      2.248351*cm, 9.613706*cm , 9.947556*cm, // G4double pDy1, G4double pDx1, G4double pDx2,
       0,// G4double pAlp1,
      2.473651*cm, 10.616205*cm , 10.983502*cm,// G4double pDy2, G4double pDx3, G4double pDx4,
      0// G4double pAlp2 //
        );
     block_trans = G4Translate3D(0.000000*cm, 104.861169*cm , 120.364250*cm) * G4RotateX3D(-0.716799*rad);

    G4Material * cylinder_mat = G4Material::GetMaterial(_geom->get_absorber_mat());
    assert(cylinder_mat);

    cylinder_logic = new G4LogicalVolume(_cylinder_solid, cylinder_mat,
        G4String(GetName().c_str()), 0, 0, 0);
    G4VisAttributes* VisAtt = new G4VisAttributes();
    PHG4Utils::SetColour(VisAtt, "W_Epoxy");
    VisAtt->SetVisibility(true);
    VisAtt->SetForceSolid(
        (not _geom->is_virualize_fiber())
            and (not _geom->is_azimuthal_seg_visible()));
    cylinder_logic->SetVisAttributes(VisAtt);

    cylinder_physi = new G4PVPlacement(block_trans, cylinder_logic,
        G4String(GetName().c_str()), logicWorld, false, 0, overlapcheck);

  }

  {

    //Processed PostionSeeds 20 from 20 1
    _cylinder_solid = new G4Trap( /*const G4String& pName*/G4String(GetName().c_str()),
     7.507337 * cm,// G4double pDz,
      0,0,// G4double pTheta, G4double pPhi,
      2.255989*cm, 9.613056*cm , 9.929544*cm, // G4double pDy1, G4double pDx1, G4double pDx2,
       0,// G4double pAlp1,
      2.483913*cm, 10.635043*cm , 10.983502*cm,// G4double pDy2, G4double pDx3, G4double pDx4,
      0// G4double pAlp2 //
        );
     block_trans = G4Translate3D(0.000000*cm, 104.861614*cm , 106.422317*cm) * G4RotateX3D(-0.778053*rad);

    G4Material * cylinder_mat = G4Material::GetMaterial(_geom->get_absorber_mat());
    assert(cylinder_mat);

    cylinder_logic = new G4LogicalVolume(_cylinder_solid, cylinder_mat,
        G4String(GetName().c_str()), 0, 0, 0);
    G4VisAttributes* VisAtt = new G4VisAttributes();
    PHG4Utils::SetColour(VisAtt, "W_Epoxy");
    VisAtt->SetVisibility(true);
    VisAtt->SetForceSolid(
        (not _geom->is_virualize_fiber())
            and (not _geom->is_azimuthal_seg_visible()));
    cylinder_logic->SetVisAttributes(VisAtt);

    cylinder_physi = new G4PVPlacement(block_trans, cylinder_logic,
        G4String(GetName().c_str()), logicWorld, false, 0, overlapcheck);

  }

//  std::pair<G4LogicalVolume *,G4Transform3D> psec = Construct_AzimuthalSeg();
//  G4LogicalVolume *sec_logic = psec.first;
//  const G4Transform3D & sec_trans = psec.second;
//
//  int n_sec_construct =
//      (_geom->is_virualize_fiber()) ? 1 : (_geom->is_azimuthal_seg_visible()?10:_geom->get_azimuthal_n_sec());
//
//  if (_geom->is_virualize_fiber() or _geom->is_azimuthal_seg_visible())
//    {
//
//      cout
//          << "PHG4SpacalDetector::Construct - WARNING - "
//          <<"only construct "<<n_sec_construct<<" sectors for visualization purpose!!!"
//          << endl;
//
//    }
//  else
//    {
//      cout
//          << "PHG4SpacalDetector::Construct - INFO - "
//          <<"start construction of "<<n_sec_construct<<" sectors."
//          << endl;
//
//    }
//
//  for (int sec = 0; sec < n_sec_construct; ++sec)
//    {
//
//      const double rot = twopi / (double) (_geom->get_azimuthal_n_sec())
//          * ((double) (sec) - n_sec_construct/2);
//
//      G4Transform3D sec_place = G4RotateZ3D(rot) * sec_trans;
//
//      stringstream name;
//      name << GetName() << "_sec" << sec;
//
//      G4PVPlacement * calo_phys = new G4PVPlacement(sec_place, sec_logic,
//          G4String(name.str().c_str()), cylinder_logic, false, sec,
//          overlapcheck);
//      calo_vol[calo_phys] = sec;
//
//    }
//  _geom->set_nscint(_geom->get_nscint() * n_sec_construct);
//
//  if (active)
//    {
//      ostringstream geonode;
//      if (superdetector != "NONE")
//        {
//          geonode << "CYLINDERGEOM_" << superdetector;
//        }
//      else
//        {
//          geonode << "CYLINDERGEOM_" << detector_type << "_" << layer;
//        }
//      PHG4CylinderGeomContainer *geo = findNode::getClass<
//          PHG4CylinderGeomContainer>(topNode, geonode.str().c_str());
//      if (!geo)
//        {
//          geo = new PHG4CylinderGeomContainer();
//          PHNodeIterator iter(topNode);
//          PHCompositeNode *runNode =
//              dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode",
//                  "RUN"));
//          PHIODataNode<PHObject> *newNode = new PHIODataNode<PHObject>(geo,
//              geonode.str().c_str(), "PHObject");
//          runNode->addNode(newNode);
//        }
//      // here in the detector class we have internal units, convert to cm
//      // before putting into the geom object
//      PHG4CylinderGeom *mygeom = clone_geom();
//      geo->AddLayerGeom(layer, mygeom);
//      //    geo->identify();
//    }
//
//
//  if (absorberactive)
//    {
//      ostringstream geonode;
//      if (superdetector != "NONE")
//        {
//          geonode << "CYLINDERGEOM_ABSORBER_" << superdetector;
//        }
//      else
//        {
//          geonode << "CYLINDERGEOM_ABSORBER_" << detector_type << "_" << layer;
//        }
//      PHG4CylinderGeomContainer *geo = findNode::getClass<
//          PHG4CylinderGeomContainer>(topNode, geonode.str().c_str());
//      if (!geo)
//        {
//          geo = new PHG4CylinderGeomContainer();
//          PHNodeIterator iter(topNode);
//          PHCompositeNode *runNode =
//              dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode",
//                  "RUN"));
//          PHIODataNode<PHObject> *newNode = new PHIODataNode<PHObject>(geo,
//              geonode.str().c_str(), "PHObject");
//          runNode->addNode(newNode);
//        }
//      // here in the detector class we have internal units, convert to cm
//      // before putting into the geom object
//      PHG4CylinderGeom *mygeom =  clone_geom();
//      geo->AddLayerGeom(layer, mygeom);
//      //    geo->identify();
//    }

  if (_geom->get_construction_verbose() >= 1)
    {
      cout << "PHG4SpacalDetector::Construct::" << GetName()
          << " - Completed. Print Geometry:" << endl;
      Print();
    }
}

std::pair<G4LogicalVolume *,G4Transform3D>
PHG4SpacalDetector::Construct_AzimuthalSeg()
{

  G4Tubs* sec_solid = new G4Tubs(G4String(GetName() + string("_sec")),
      _geom->get_radius() * cm, _geom->get_max_radius() * cm,
      _geom->get_length() * cm / 2.0, 0, twopi / _geom->get_azimuthal_n_sec());

  G4Material * cylinder_mat = G4Material::GetMaterial(_geom->get_absorber_mat());
  assert(cylinder_mat);

  G4LogicalVolume * sec_logic = new G4LogicalVolume(sec_solid, cylinder_mat,
      G4String(G4String(GetName() + string("_sec"))), 0, 0, step_limits);

  G4VisAttributes* VisAtt = new G4VisAttributes();
  VisAtt->SetColor(.1, .1, .1, .5);
  VisAtt->SetVisibility(_geom->is_virualize_fiber());
  VisAtt->SetForceSolid(false);
  sec_logic->SetVisAttributes(VisAtt);

  const double fiber_length = _geom->get_thickness() * cm
      - 2 * _geom->get_fiber_outer_r() * cm;
  G4LogicalVolume *fiber_logic = Construct_Fiber(fiber_length, string(""));

  int fiber_count = 0;
//  double z_step = _geom->get_fiber_distance() * cm * sqrt(3) / 2.;
  double z_step = _geom->get_z_distance() * cm;
  G4double z = _geom->get_zmin() * cm - _geom->get_zpos() * cm + z_step;
  while (z < _geom->get_zmax() * cm - _geom->get_zpos() * cm - z_step)
    {

      const double rot = twopi / _geom->get_azimuthal_n_sec()
          * ((fiber_count % 2 == 0) ? 1. / 4. : 3. / 4.);

      G4Transform3D fiber_place(
          G4RotateZ3D(rot) * G4TranslateZ3D(z)
              * G4TranslateX3D(_geom->get_half_radius() * cm)
              * G4RotateY3D(halfpi));

      stringstream name;
      name << GetName() << "_fiber_" << fiber_count;

      G4PVPlacement * fiber_physi = new G4PVPlacement(fiber_place, fiber_logic,
          G4String(name.str().c_str()), sec_logic, false, fiber_count,
          overlapcheck);
      fiber_vol[fiber_physi] = fiber_count;

      z += z_step;
      fiber_count++;
    }
  _geom->set_nscint(fiber_count);

  cout << "PHG4SpacalDetector::Construct_AzimuthalSeg::" << GetName()
      << " - constructed " << fiber_count << " fibers" << endl;
  cout << "\t" << "_geom->get_fiber_distance() = " << _geom->get_fiber_distance()
      << endl;
  cout << "\t" << "fiber_length = " << fiber_length / cm << endl;
  cout << "\t" << "z_step = " << z_step << endl;
  cout << "\t" << "_geom->get_azimuthal_bin() = " << _geom->get_azimuthal_n_sec()
      << endl;
  cout << "\t" << "_geom->get_azimuthal_distance() = "
      << _geom->get_azimuthal_distance() << endl;
  cout << "\t" << "_geom->is_virualize_fiber() = " << _geom->is_virualize_fiber()
      << endl;

  return make_pair(sec_logic,G4Transform3D::Identity);
}

G4LogicalVolume *
PHG4SpacalDetector::Construct_Fiber(const G4double length, const string & id)
{

  G4Tubs* fiber_solid = new G4Tubs(G4String(GetName() + string("_fiber") + id),
      0, _geom->get_fiber_outer_r() * cm, length / 2.0, 0, twopi);

  G4Material * clading_mat = G4Material::GetMaterial(
      _geom->get_fiber_clading_mat());
  assert(clading_mat);

  G4LogicalVolume * fiber_logic = new G4LogicalVolume(fiber_solid, clading_mat,
      G4String(G4String(GetName() + string("_fiber") + id)), 0, 0,
      clading_step_limits);

    {
      G4VisAttributes* VisAtt = new G4VisAttributes();
      PHG4Utils::SetColour(VisAtt, "G4_POLYSTYRENE");
      VisAtt->SetVisibility(_geom->is_virualize_fiber());
      VisAtt->SetForceSolid(_geom->is_virualize_fiber());
      fiber_logic->SetVisAttributes(VisAtt);
    }

  G4Tubs* core_solid = new G4Tubs(
      G4String(GetName() + string("_fiber_core") + id), 0,
      _geom->get_fiber_core_diameter() * cm / 2, length / 2.0, 0, twopi);

  G4Material * core_mat = G4Material::GetMaterial(_geom->get_fiber_core_mat());
  assert(core_mat);

  G4LogicalVolume * core_logic = new G4LogicalVolume(core_solid, core_mat,
      G4String(G4String(GetName() + string("_fiber_core") + id)), 0, 0,
      fiber_core_step_limits);

    {
      G4VisAttributes* VisAtt = new G4VisAttributes();
      PHG4Utils::SetColour(VisAtt, "G4_POLYSTYRENE");
      VisAtt->SetVisibility(false);
      VisAtt->SetForceSolid(false);
      core_logic->SetVisAttributes(VisAtt);
    }

  G4PVPlacement * core_physi = new G4PVPlacement(0, G4ThreeVector(), core_logic,
      G4String(G4String(GetName() + string("_fiber_core") + id)), fiber_logic,
      false, 0, overlapcheck);
  fiber_core_vol[core_physi] = 0;

  return fiber_logic;
}

void
PHG4SpacalDetector::Print(const std::string &what) const
{
  cout << "PHG4SpacalDetector::Print::" << GetName() << " - Print Geometry:"
      << endl;
  _geom->Print();

  return;
}
