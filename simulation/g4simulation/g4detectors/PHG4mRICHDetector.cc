/*===============================================================*
 *                        March 2nd 2017                         *
 *         mRICH Detector created by Cheuk-Ping Wong @GSU        *
 *===============================================================*
 * March 2nd 2017. Modified from PHG4BlockDetector.cc. This code *
 *                 still need some small modifications for mRICH *
 *===============================================================*/
#include "PHG4mRICHDetector.h"
#include "PHG4Parameters.h"

#include <g4main/PHG4Utils.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/getClass.h>

#include <Geant4/globals.hh>
#include <Geant4/G4Material.hh>
#include <Geant4/G4NistManager.hh>
#include <Geant4/G4Box.hh>
#include <Geant4/G4Polyhedra.hh>
#include <Geant4/G4Polycone.hh>
#include <Geant4/G4Cons.hh>
#include <Geant4/G4Sphere.hh>
#include <Geant4/G4Tubs.hh>
#include <Geant4/G4VPhysicalVolume.hh>
#include <Geant4/G4LogicalVolume.hh>
#include <Geant4/G4SubtractionSolid.hh>
#include <Geant4/G4IntersectionSolid.hh>
#include <Geant4/G4PVPlacement.hh>
#include <Geant4/G4SystemOfUnits.hh>
#include <Geant4/G4VisAttributes.hh>
#include <Geant4/G4Colour.hh>
#include <Geant4/G4LogicalBorderSurface.hh>
#include <Geant4/G4OpticalSurface.hh>
#include <Geant4/G4RotationMatrix.hh>

#include <sstream>
#include <string>


using namespace std;
using namespace CLHEP;

//_______________________________________________________________
PHG4mRICHDetector::PHG4mRICHDetector( PHCompositeNode *Node, PHG4Parameters *parameters, const std::string &dnam, const int lyr):
  PHG4Detector(Node, dnam),
  params(parameters),
  //block_physi(NULL),
  layer(lyr)
{}

//_______________________________________________________________
bool PHG4mRICHDetector::IsInmRICH(G4VPhysicalVolume * volume) const
{
  if ( strcmp(volume->GetName(),"sensor")==0 ) 
  {
    return true;
  }
  return false;
}
//______________________________________________________________
void PHG4mRICHDetector::Construct( G4LogicalVolume* logicWorld )
{
  int single_mRICH=0;
  G4double bowlPar[4];
  
  if (single_mRICH) Construct_a_mRICH(logicWorld);
  else {
    //build_Space(logicWorld,bowlPar);
    G4LogicalVolume* space = build_Space(logicWorld,bowlPar);
    G4LogicalVolume* a_mRICH=Construct_a_mRICH(0);
    build_mRICH_wall(space,a_mRICH,bowlPar);
  }
}
//_______________________________________________________________
G4LogicalVolume* PHG4mRICHDetector::Construct_a_mRICH( G4LogicalVolume* logicWorld )
{
  mRichMaterialList* materialList=new mRichMaterialList();
  mRichParameter* parameters=new mRichParameter(materialList);

  /*holder box and hollow volume*/ G4VPhysicalVolume* hollowVol=build_holderBox(parameters,logicWorld);
  /*foam holder for aerogel     */ build_foamHolder(parameters,hollowVol->GetLogicalVolume());
  /*aerogel                     */ build_aerogel(parameters,hollowVol);
  /*lens                        */ //build_lens(parameters->GetLensPar("fresnelLens"), hollowVol->GetLogicalVolume());
  /*mirror                      */ //build_mirror(parameters,hollowVol);
  /*sensor plane                */ //build_sensor(parameters,hollowVol->GetLogicalVolume());
  /*readout electronics         */ //build_polyhedra(parameters->GetPolyPar("readout"),hollowVol->GetLogicalVolume());

  return hollowVol->GetMotherLogical();  //return detector holder box.
                                         //you have more than 1 daugthers,
                                         //but you can only have one mother.

  printf("============== detector built ================\n");
  
  
}

//________________________________________________________________________//

mRichMaterialList::mRichMaterialList()
{
  G4NistManager * nist = G4NistManager::Instance();

  //Gean4 Predefined
  Water= nist->FindOrBuildMaterial("G4_Water");
  Aluminum= nist->FindOrBuildMaterial("G4_Al");

  SetAir();
  SetAir_Opt();
  SetAcrylic();
  SetAerogel1();
  SetAerogel2();
  SetBorosilicate();
}
//________________________________________________________________________//
mRichMaterialList::~mRichMaterialList(){;}
//________________________________________________________________________//
void mRichMaterialList::SetAir()
{
  const int nEntries1=32;
  //using same photon energy array as aerogel 1
  G4double PhotonEnergy[nEntries1] =
    { 2.034*eV, 2.068*eV, 2.103*eV, 2.139*eV,     // 610, 600, 590, 580, (nm)  
      2.177*eV, 2.216*eV, 2.256*eV, 2.298*eV,     // 570, 560, 550, 540,
      2.341*eV, 2.386*eV, 2.433*eV, 2.481*eV,     // 530, 520, 510, 500,
      2.532*eV, 2.585*eV, 2.640*eV, 2.697*eV,     // 490, 480, 470, 460,
      2.757*eV, 2.820*eV, 2.885*eV, 2.954*eV,     // 450, 440, 430, 420,
      3.026*eV, 3.102*eV, 3.181*eV, 3.265*eV,     // 410, 400, 390, 380,  
      3.353*eV, 3.446*eV, 3.545*eV, 3.649*eV,     // 370, 360, 350, 340,
      3.760*eV, 3.877*eV, 4.002*eV, 4.136*eV };   // 330, 320, 310, 300.

  G4double AirRefractiveIndex[nEntries1] =
    { 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
      1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
      1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
      1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
      1.00, 1.00, 1.00, 1.00 };

  G4MaterialPropertiesTable* myMPT = new G4MaterialPropertiesTable();
  myMPT->AddProperty("RINDEX", PhotonEnergy, AirRefractiveIndex , nEntries1);

  G4NistManager * nist = G4NistManager::Instance();
  Air= nist->FindOrBuildMaterial("G4_AIR");
  Air->SetMaterialPropertiesTable(myMPT);
}
//________________________________________________________________________//
void mRichMaterialList::SetAir_Opt()
{
  const int nEntries_Air_Opt=2;

  G4Element* N = new G4Element("Nitrogen", "N",  7 , 14.01*g/mole);
  G4Element* O = new G4Element("Oxygen"  , "O",  8 , 16.00*g/mole);

  G4double PhotonEnergy_Air_Opt[nEntries_Air_Opt] = {2.034*eV, 4.136*eV};
  G4double RefractiveIndex_Air_Opt[nEntries_Air_Opt] = {1.00, 1.00};

  G4MaterialPropertiesTable* myMPT = new G4MaterialPropertiesTable();
  myMPT->AddProperty("RINDEX", PhotonEnergy_Air_Opt, RefractiveIndex_Air_Opt, nEntries_Air_Opt);

  Air_Opt = new G4Material("Air_Opt", 1.29*mg/cm3, 2);
  Air_Opt->AddElement(N, 70.*perCent);
  Air_Opt->AddElement(O, 30.*perCent);
  Air_Opt->SetMaterialPropertiesTable(myMPT);
}
//________________________________________________________________________//
void mRichMaterialList::SetAcrylic()
{
  const int nEntries1=32;

  G4Element* H = new G4Element("Hydrogen", "H",  1 , 1.01*g/mole);
  G4Element* O = new G4Element("Oxygen"  , "O",  8 , 16.00*g/mole);
  G4Element* C = new G4Element("Carbon"  ,"C" , 6., 12.01*g/mole);

  //same photon energy array as aerogel 1
  G4double PhotonEnergy[nEntries1] =
    { 2.034*eV, 2.068*eV, 2.103*eV, 2.139*eV,     // 610, 600, 590, 580, (nm)
      2.177*eV, 2.216*eV, 2.256*eV, 2.298*eV,     // 570, 560, 550, 540,
      2.341*eV, 2.386*eV, 2.433*eV, 2.481*eV,     // 530, 520, 510, 500,  
      2.532*eV, 2.585*eV, 2.640*eV, 2.697*eV,     // 490, 480, 470, 460,
      2.757*eV, 2.820*eV, 2.885*eV, 2.954*eV,     // 450, 440, 430, 420,  
      3.026*eV, 3.102*eV, 3.181*eV, 3.265*eV,     // 410, 400, 390, 380,
      3.353*eV, 3.446*eV, 3.545*eV, 3.649*eV,     // 370, 360, 350, 340,  
      3.760*eV, 3.877*eV, 4.002*eV, 4.136*eV };   // 330, 320, 310, 300.

  G4double AcRefractiveIndex[nEntries1] =
    { 1.4902, 1.4907, 1.4913, 1.4918, 1.4924,   // 610, 600, 590, 580, 570,
      1.4930, 1.4936, 1.4942, 1.4948, 1.4954,   // 560, 550, 540, 530, 520,  (this line is interpolated)
      1.4960, 1.4965, 1.4971, 1.4977, 1.4983,   // 510, 500, 490, 480, 470,
      1.4991, 1.5002, 1.5017, 1.5017, 1.5017,   // 460, 450, 440, 430, 420,
      1.5017, 1.5017, 1.5017, 1.5017, 1.5017,   // 410,  
      1.5017, 1.5017, 1.5017, 1.5017, 1.5017,   // 360,     look up values below 435 
      1.5017, 1.5017};

  G4double AcAbsorption[nEntries1] =
    {25.25*cm, 25.25*cm, 25.25*cm, 25.25*cm,
     25.25*cm, 25.25*cm, 25.25*cm, 25.25*cm,
     25.25*cm, 25.25*cm, 25.25*cm, 25.25*cm,
     25.25*cm, 25.25*cm, 25.25*cm, 25.25*cm,
     25.25*cm, 25.25*cm, 25.25*cm, 25.25*cm,
     25.25*cm, 00.667*cm, 00.037*cm, 00.333*cm,
     00.001*cm, 00.001*cm, 00.001*cm, 00.001*cm,
     00.001*cm, 00.001*cm, 00.001*cm, 00.001*cm};

  G4MaterialPropertiesTable* myMPT = new G4MaterialPropertiesTable();
  myMPT->AddProperty("RINDEX" ,PhotonEnergy, AcRefractiveIndex,nEntries1);
  myMPT->AddProperty("ABSLENGTH", PhotonEnergy, AcAbsorption,nEntries1);

  Acrylic=new G4Material("Acrylic", 1.19*g/cm3, 3);
  Acrylic->AddElement(C, 5);
  Acrylic->AddElement(H, 8);     // molecular ratios
  Acrylic->AddElement(O, 2);
  Acrylic->SetMaterialPropertiesTable(myMPT);
}
//________________________________________________________________________//
void mRichMaterialList::SetAerogel1()
{
  const int nEntries1=32;

  G4Element* O = new G4Element("Oxygen"  , "O",  8 , 16.00*g/mole);
  G4Element* Si= new G4Element("Silicon", "Si", 14, 28.00*g/mole);

  G4double Agel1PhotonEnergy[nEntries1] =
    { 2.034*eV, 2.068*eV, 2.103*eV, 2.139*eV,     // 610, 600, 590, 580, (nm) 
      2.177*eV, 2.216*eV, 2.256*eV, 2.298*eV,     // 570, 560, 550, 540,
      2.341*eV, 2.386*eV, 2.433*eV, 2.481*eV,     // 530, 520, 510, 500,     
      2.532*eV, 2.585*eV, 2.640*eV, 2.697*eV,     // 490, 480, 470, 460,            
      2.757*eV, 2.820*eV, 2.885*eV, 2.954*eV,     // 450, 440, 430, 420,
      3.026*eV, 3.102*eV, 3.181*eV, 3.265*eV,     // 410, 400, 390, 380,                               
      3.353*eV, 3.446*eV, 3.545*eV, 3.649*eV,     // 370, 360, 350, 340,
      3.760*eV, 3.877*eV, 4.002*eV, 4.136*eV };   // 330, 320, 310, 300.

  G4double Agel1RefractiveIndex[nEntries1] =
    { 1.02435, 1.0244,  1.02445, 1.0245,  1.02455,
      1.0246,  1.02465, 1.0247,  1.02475, 1.0248,
      1.02485, 1.02492, 1.025,   1.02505, 1.0251,
      1.02518, 1.02522, 1.02530, 1.02535, 1.0254,
      1.02545, 1.0255,  1.02555, 1.0256,  1.02568,
      1.02572, 1.0258,  1.02585, 1.0259,  1.02595,
      1.026,   1.02608};

  G4double Agel1Absorption[nEntries1] =    //from Hubert                                               
    { 3.448*m,  4.082*m,  6.329*m,  9.174*m, 12.346*m, 13.889*m,
      15.152*m, 17.241*m, 18.868*m, 20.000*m, 26.316*m, 35.714*m,
      45.455*m, 47.619*m, 52.632*m, 52.632*m, 55.556*m, 52.632*m,
      52.632*m, 47.619*m, 45.455*m, 41.667*m, 37.037*m, 33.333*m,
      30.000*m, 28.500*m, 27.000*m, 24.500*m, 22.000*m, 19.500*m,
      17.500*m, 14.500*m };

  G4double Agel1Rayleigh[nEntries1];
  //const G4double AerogelTypeAClarity = 0.00719*micrometer*micrometer*micrometer*micrometer/cm;
  const G4double AerogelTypeAClarity = 0.0020*micrometer*micrometer*micrometer*micrometer/cm;
  G4double Cparam    =  AerogelTypeAClarity*cm/(micrometer*micrometer*micrometer*micrometer);
  G4double PhotMomWaveConv = 1239*eV*nm;

  if(Cparam != 0.0 ) {
    for(int i=0; i<nEntries1; i++ ){
      G4double ephoton = Agel1PhotonEnergy[i];
      //In the following the 1000 is to convert form nm to micrometer
      G4double wphoton=(PhotMomWaveConv/ephoton)/(1000.0*nm);
      Agel1Rayleigh[i]=(std::pow(wphoton,4))/Cparam;
    }
  }

  G4MaterialPropertiesTable* myMPT = new G4MaterialPropertiesTable();
  myMPT->AddProperty("RINDEX", Agel1PhotonEnergy, Agel1RefractiveIndex, nEntries1);
  myMPT->AddProperty("ABSLENGTH", Agel1PhotonEnergy, Agel1Absorption,nEntries1);
  myMPT->AddProperty("RAYLEIGH", Agel1PhotonEnergy, Agel1Rayleigh, nEntries1);  //Need table of rayleigh Scattering!!!
  myMPT->AddConstProperty("SCINTILLATIONYIELD",0./MeV);
  myMPT->AddConstProperty("RESOLUTIONSCALE",1.0);

  Aerogel1 = new G4Material("Aerogel1", 0.02*g/cm3, 2);
  Aerogel1->AddElement(Si, 1);
  Aerogel1->AddElement(O, 2);

  Aerogel1->SetMaterialPropertiesTable(myMPT);
}
//________________________________________________________________________//
void mRichMaterialList::SetAerogel2()
{
  const int nEntries2=50;

  G4Element* O = new G4Element("Oxygen"  , "O",  8 , 16.00*g/mole);
  G4Element* Si= new G4Element("Silicon", "Si", 14, 28.00*g/mole);

  G4double Agel2PhotonEnergy[nEntries2]=
    {1.87855*eV,1.96673*eV,2.05490*eV,2.14308*eV,2.23126*eV,
     2.31943*eV,2.40761*eV,2.49579*eV,2.58396*eV,2.67214*eV,
     2.76032*eV,2.84849*eV,2.93667*eV,3.02485*eV,3.11302*eV,
     3.20120*eV,3.28938*eV,3.37755*eV,3.46573*eV,3.55391*eV,
     3.64208*eV,3.73026*eV,3.81844*eV,3.90661*eV,3.99479*eV,
     4.08297*eV,4.17114*eV,4.25932*eV,4.34750*eV,4.43567*eV,
     4.52385*eV,4.61203*eV,4.70020*eV,4.78838*eV,4.87656*eV,
     4.96473*eV,5.05291*eV,5.14109*eV,5.22927*eV,5.31744*eV,
     5.40562*eV,5.49380*eV,5.58197*eV,5.67015*eV,5.75833*eV,
     5.84650*eV,5.93468*eV,6.02286*eV,6.11103*eV,6.19921*eV };

  G4double Agel2RefractiveIndex[nEntries2] =
    {1.02825,1.02829,1.02834,1.02839,1.02844,
     1.02849,1.02854,1.02860,1.02866,1.02872,
     1.02878,1.02885,1.02892,1.02899,1.02906,
     1.02914,1.02921,1.02929,1.02938,1.02946,
     1.02955,1.02964,1.02974,1.02983,1.02993,
     1.03003,1.03014,1.03025,1.03036,1.03047,
     1.03059,1.03071,1.03084,1.03096,1.03109,
     1.03123,1.03137,1.03151,1.03166,1.03181,
     1.03196,1.03212,1.03228,1.03244,1.03261,
     1.03279,1.03297,1.03315,1.03334,1.03354};

  G4double Agel2Absorption[nEntries2] =      //from Marco                                               
    {17.5000*cm,17.7466*cm,17.9720*cm,18.1789*cm,18.3694*cm,
     18.5455*cm,18.7086*cm,18.8602*cm,19.0015*cm,19.1334*cm,
     19.2569*cm,19.3728*cm,19.4817*cm,19.5843*cm,19.6810*cm,
     19.7725*cm,19.8590*cm,19.9410*cm,20.0188*cm,20.0928*cm,
     18.4895*cm,16.0174*cm,13.9223*cm,12.1401*cm,10.6185*cm,
     9.3147*cm,8.1940*cm,7.2274*cm,6.3913*cm,5.6659*cm,
     5.0347*cm,4.4841*cm,4.0024*cm,3.5801*cm,3.2088*cm,
     2.8817*cm,2.5928*cm,2.3372*cm,2.1105*cm,1.9090*cm,
     1.7296*cm,1.5696*cm,1.4266*cm,1.2986*cm,1.1837*cm,
     1.0806*cm,0.9877*cm,0.9041*cm,0.8286*cm,0.7603*cm };

  G4double Agel2Rayleigh[nEntries2] =         //from Marco                                              
    {35.1384*cm, 29.24805*cm, 24.5418*cm, 20.7453*cm, 17.6553*cm,
     15.1197*cm, 13.02345*cm, 11.2782*cm, 9.81585*cm, 8.58285*cm,
     7.53765*cm, 6.6468*cm, 5.88375*cm, 5.22705*cm, 4.6596*cm,
     4.167*cm, 3.73785*cm, 3.36255*cm, 3.03315*cm, 2.7432*cm,
     2.487*cm, 2.26005*cm, 2.05845*cm, 1.87875*cm, 1.71825*cm,
     1.57455*cm, 1.44555*cm, 1.3296*cm, 1.2249*cm, 1.1304*cm,
     1.04475*cm, 0.9672*cm, 0.89655*cm, 0.83235*cm, 0.77385*cm,
     0.7203*cm, 0.67125*cm, 0.6264*cm, 0.58515*cm, 0.54735*cm,
     0.51255*cm, 0.48045*cm, 0.45075*cm, 0.4233*cm, 0.39795*cm,
     0.37455*cm, 0.3528*cm, 0.33255*cm, 0.3138*cm, 0.29625*cm};


  G4MaterialPropertiesTable* myMPT = new G4MaterialPropertiesTable();
  myMPT->AddProperty("RINDEX", Agel2PhotonEnergy, Agel2RefractiveIndex, nEntries2);
  myMPT->AddProperty("ABSLENGTH", Agel2PhotonEnergy, Agel2Absorption,nEntries2);
  myMPT->AddProperty("RAYLEIGH", Agel2PhotonEnergy, Agel2Rayleigh, nEntries2);
  myMPT->AddConstProperty("SCINTILLATIONYIELD",0./MeV);
  myMPT->AddConstProperty("RESOLUTIONSCALE",1.0);

  Aerogel2 = new G4Material("Aerogel2", 0.02*g/cm3, 2);
  Aerogel2->AddElement(Si, 1);
  Aerogel2->AddElement(O, 2);
  Aerogel2->SetMaterialPropertiesTable(myMPT);
}
//________________________________________________________________________//
void mRichMaterialList::SetBorosilicate()
{
  const int nEntries1=32;

  G4double PhotonEnergy[nEntries1] =
    { 2.034*eV, 2.068*eV, 2.103*eV, 2.139*eV,     // 610, 600, 590, 580, (nm) 
      2.177*eV, 2.216*eV, 2.256*eV, 2.298*eV,     // 570, 560, 550, 540,
      2.341*eV, 2.386*eV, 2.433*eV, 2.481*eV,     // 530, 520, 510, 500,                                
      2.532*eV, 2.585*eV, 2.640*eV, 2.697*eV,     // 490, 480, 470, 460,
      2.757*eV, 2.820*eV, 2.885*eV, 2.954*eV,     // 450, 440, 430, 420,
      3.026*eV, 3.102*eV, 3.181*eV, 3.265*eV,     // 410, 400, 390, 380,
      3.353*eV, 3.446*eV, 3.545*eV, 3.649*eV,     // 370, 360, 350, 340,
      3.760*eV, 3.877*eV, 4.002*eV, 4.136*eV };   // 330, 320, 310, 300.

  G4double glassRefractiveIndex[nEntries1] =
    { 1.47, 1.47,  1.47, 1.47,  1.47,
      1.47, 1.47,  1.47, 1.47,  1.47,
      1.47, 1.47,  1.47, 1.47,  1.47,
      1.47, 1.47,  1.47, 1.47,  1.47,
      1.47, 1.47,  1.47, 1.47,  1.47,
      1.47, 1.47,  1.47, 1.47,  1.47,
      1.47, 1.47};

  G4double glassAbsorption[nEntries1] =
    {4.25*cm, 4.25*cm, 4.25*cm, 4.25*cm,
     4.25*cm, 4.25*cm, 4.25*cm, 4.25*cm,
     4.25*cm, 4.25*cm, 4.25*cm, 4.25*cm,
     4.25*cm, 4.25*cm, 4.25*cm, 4.25*cm,
     4.25*cm, 4.25*cm, 4.25*cm, 4.25*cm,
     4.25*cm, 00.667*cm, 00.037*cm, 00.333*cm,
     00.001*cm, 00.001*cm, 00.001*cm, 00.001*cm,
     00.001*cm, 00.001*cm, 00.001*cm, 00.001*cm};

  G4MaterialPropertiesTable* myMPT = new G4MaterialPropertiesTable();
  //same photon energy array as aerogel 1
  myMPT->AddProperty("RINDEX", PhotonEnergy, glassRefractiveIndex, nEntries1);
  myMPT->AddProperty("ABSLENGTH", PhotonEnergy, glassAbsorption, nEntries1);

  G4NistManager * nist = G4NistManager::Instance();
  Borosilicate = nist->FindOrBuildMaterial("G4_Pyrex_Glass");
  Borosilicate->SetMaterialPropertiesTable(myMPT);
}
//________________________________________________________________________//
G4Material* mRichMaterialList::GetmRichMaterial(const char* matName)
{
  if (strcmp(matName,"Water")==0) return Water;
  else if (strcmp(matName,"Aluminum")==0) return Aluminum;
  else if (strcmp(matName,"Air")==0) return Air;
  else if (strcmp(matName,"Air_Opt")==0) return Air_Opt;
  else if (strcmp(matName,"Acrylic")==0) return Acrylic;
  else if (strcmp(matName,"Aerogel1")==0) return Aerogel1;
  else if (strcmp(matName,"Aerogel2")==0) return Aerogel2;
  else if (strcmp(matName,"Borosilicate")==0) return Borosilicate;
  else G4cout<<"mRichMaterialList::GetmRichMaterial() ----- ERROR: cannot find material="<<matName<<G4endl;
  return 0;
}
//________________________________________________________________________//
BoxPar::BoxPar() {;}
//________________________________________________________________________//
BoxPar::~BoxPar() {;}
//________________________________________________________________________//
PolyPar::PolyPar() {;}
//________________________________________________________________________//
PolyPar::~PolyPar() {;}
//________________________________________________________________________//
LensPar::LensPar()
{
  n=0;
  f=0;
  diameter=0;
  eff_diameter=0;
  centerThickness=0;
  grooveWidth=0;
}
//________________________________________________________________________//
LensPar::~LensPar() {;}
//________________________________________________________________________//
void LensPar::Set_halfXYZ(G4double halfX, G4double grooveDensity)
{
  halfXYZ[0]=halfX;
  halfXYZ[1]=halfXYZ[0];

  G4double NumberOfGrooves=floor(grooveDensity*(eff_diameter/2.0));
  G4double Rmin1 = (NumberOfGrooves-1)*(grooveWidth) ;
  G4double Rmax1 = (NumberOfGrooves-0)*(grooveWidth) ;
  halfXYZ[2] = (GetSagita(Rmax1)-GetSagita(Rmin1)+centerThickness)/2.0;
}
//________________________________________________________________________//
G4double LensPar::GetSagita(G4double r)
{
  G4double Conic = -1.0;    // original:
  //G4int lens_type = 3; 
  G4int lens_type = 5;
  G4double Curvature;
  G4double Aspher[8] = {0,0,0,0,0,0,0,0};

  if (lens_type == 1) {
    Curvature = 0.00437636761488/mm;
    Aspher[0] = 4.206739256e-05/(mm);
    Aspher[1] = 9.6440152e-10/(mm3);
    Aspher[2] = -1.4884317e-15/(mm2*mm3);
  }
  else  if (lens_type == 2) {
    Curvature = 0.0132/mm;                // r=77mm, f~14cm
    Aspher[0] = 32.0e-05/(mm);
    Aspher[1] = -2.0e-7/(mm3);
    Aspher[2] =  1.2e-13/(mm2*mm3);
  }
  else  if (lens_type == 3) {
    Curvature = 0.0150/mm;         // r=77mm, f~12.5cm
    Aspher[0] = 42.0e-05/(mm);
    Aspher[1] = -3.0e-7/(mm3);
    Aspher[2] =  1.2e-13/(mm2*mm3);
  }
  else if (lens_type == 4) {
    Curvature = 0.0175/mm;         // r=77mm, f~10cm
    Aspher[0] = 72.0e-05/(mm);
    Aspher[1] = -5.0e-7/(mm3);
    Aspher[2] =  1.2e-13/(mm2*mm3);
  }
  else  if (lens_type==5) {
    Curvature=1/(f*(n-1));
  }

  G4double TotAspher = 0.0*mm ;

  for(G4int k=1;k<9;k++){ TotAspher += Aspher[k-1]*std::pow(r,2*k); }

  G4double ArgSqrt = 1.0-(1.0+Conic)*std::pow(Curvature,2)*std::pow(r,2) ; // note conic=-1, so ArgSqrt = 1.0

  if (ArgSqrt < 0.0){
    G4cout << "UltraFresnelLensParameterisation::Sagita: Square Root of <0 !" << G4endl;
  }
  G4double Sagita_value = Curvature*std::pow(r,2)/(1.0+std::sqrt(ArgSqrt)) + TotAspher;

  return Sagita_value ;
}

//________________________________________________________________________//
mRichParameter::mRichParameter(mRichMaterialList* materialList)
{
  int i;
  //----------
  // Constant 
  //----------

  const double myPI=4*atan(1);
  fresnelLens=new LensPar();
  holderBox=new BoxPar();
  hollowVolume=new BoxPar();
  foamHolderBox=new BoxPar();
  foamHolderPoly=new PolyPar();
  aerogel=new BoxPar();
  mirror=new PolyPar();
  glassWindow=new BoxPar();
  sensor=new BoxPar();
  readout=new PolyPar();

  //--------------------Holder box key parameters-----------------------------//
  const G4double BoxDelz = 2.0*mm;   // extract space between components
  const G4double box_thicknessXYZ[4] = {(1./4.)*2.54*cm,(1./4.)*2.54*cm,(1./16.)*2.54*cm,(1./4.)*2.54*cm};

  //------------------------- Aerogel gel key parameters---------------------//
  const G4double foamHolderThicknessXYZ[3]={1.0*cm,1.0*cm,1.0*cm};
  const G4double agel_halfXYZ[3]={5.525*cm,5.525*cm,1.65*cm};

  //------------------------Fresnel lens key parameters----------------------//
  const G4double lens_gap=(2.54/8.0)*cm;    //gap between agel and lens, and between lens and mirror
  
  const G4double lensHalfx=(5.25/2.0)*2.54*cm;
  const G4double grooveDensity=125.0/(2.54*cm);

  fresnelLens->n=1.49;
  fresnelLens->f=6.0*2.54*cm;
  fresnelLens->eff_diameter=15.24*cm;
  fresnelLens->diameter=2.0*sqrt(2.0)*lensHalfx;
  fresnelLens->centerThickness=0.06*2.54*cm;
  fresnelLens->grooveWidth=(G4double) 1.0/grooveDensity;
  fresnelLens->Set_halfXYZ(lensHalfx,grooveDensity);
  //rest of lens parameters are set below

  //---------------------------Photodetector key parameters-------------------//
  const G4double sensorGap=0.05*cm;      //half width of the gap
  const G4double glassWindow_halfXYZ[3]={5.2/2.0*cm,5.2/2.0*cm,0.075*cm};
  const G4double phodet_halfXYZ[3]={2.4*cm,2.4*cm,0.075*cm};

  //--------------------------- mirror key parameters ------------------------//
  const G4double mirrorThickness=0.2*cm;

  //-----------------------Readout Electronics key parameters-----------------//
  const G4double readout_halfz = 0.4*cm;    //redendunt?
  const G4double readoutThickness=0.2*cm;

  //----------
  // calculation
  //----------
  G4double sensor_total_halfx=2*glassWindow_halfXYZ[0]+sensorGap;

  G4double foamHolder_halfXYZ[3];
  foamHolder_halfXYZ[0]=agel_halfXYZ[0]+foamHolderThicknessXYZ[0];
  foamHolder_halfXYZ[1]=foamHolder_halfXYZ[0];
  foamHolder_halfXYZ[2]=foamHolderThicknessXYZ[2]/2.0;

  G4double acrylicBox_halfXYZ[3];
  acrylicBox_halfXYZ[0] = max(max(foamHolder_halfXYZ[0],sensor_total_halfx+readoutThickness), fresnelLens->halfXYZ[0])+0.1\
    *cm + box_thicknessXYZ[0];
  acrylicBox_halfXYZ[1] = acrylicBox_halfXYZ[0];
  acrylicBox_halfXYZ[2] = (BoxDelz+2*foamHolder_halfXYZ[2]+2*agel_halfXYZ[2]
                           +lens_gap+2*fresnelLens->halfXYZ[2]+fresnelLens->f+2*glassWindow_halfXYZ[2]
                           +2*phodet_halfXYZ[2]+(2*readout_halfz+BoxDelz)
                           +box_thicknessXYZ[2]+box_thicknessXYZ[3])/2.0;


  G4double hollow_halfXYZ[3];
  hollow_halfXYZ[0]=acrylicBox_halfXYZ[0]-box_thicknessXYZ[0];
  hollow_halfXYZ[1]=hollow_halfXYZ[0];
  hollow_halfXYZ[2]=(2*acrylicBox_halfXYZ[2]-box_thicknessXYZ[2]-box_thicknessXYZ[3])/2.0;

  G4ThreeVector hollow_posXYZ=G4ThreeVector(0.0*cm,0.0*cm,-acrylicBox_halfXYZ[2]+hollow_halfXYZ[2]+box_thicknessXYZ[2]);

  G4double foamHolder_posz=-hollow_halfXYZ[2]+BoxDelz+foamHolder_halfXYZ[2];
  G4double agel_posz=foamHolder_posz+foamHolder_halfXYZ[2]+agel_halfXYZ[2];
  G4double lens_z=agel_posz+agel_halfXYZ[2]+fresnelLens->halfXYZ[2]+lens_gap;

  G4double glassWindow_z=lens_z-fresnelLens->halfXYZ[2]+fresnelLens->f+glassWindow_halfXYZ[2]; //out of focus. But this makes sense.
  G4double phodet_z=glassWindow_z+glassWindow_halfXYZ[2]+phodet_halfXYZ[2];

  //redendunt:
  G4double readout_z[2];
  readout_z[0]=glassWindow_z-glassWindow_halfXYZ[2];
  readout_z[1]=phodet_z+phodet_halfXYZ[2];

  //----------
  // set holderBox
  //----------
  sprintf(holderBox->name,"HolderBox");
  for (i=0;i<3;i++) holderBox->halfXYZ[i]=acrylicBox_halfXYZ[i];
  holderBox->posXYZ=G4ThreeVector(0*cm,0*cm,0*cm);
  holderBox->material=materialList->GetmRichMaterial("Aluminum");
  holderBox->sensitivity=0;

  holderBox->color=G4Colour(1.0,1.0,1.0);
  holderBox->visibility=true;
  holderBox->wireframe=true;
  holderBox->surface=false;
  
  //----------
  // set HollowVolume
  //----------
  sprintf(hollowVolume->name,"HollowVolume");
  for (i=0;i<3;i++) hollowVolume->halfXYZ[i]=hollow_halfXYZ[i];
  hollowVolume->posXYZ=hollow_posXYZ;
  hollowVolume->material=materialList->GetmRichMaterial("Air_Opt");
  hollowVolume->sensitivity=0;

  hollowVolume->color=G4Colour(1.0,1.0,1.0);
  hollowVolume->visibility=true;
  hollowVolume->wireframe=true;
  hollowVolume->surface=false;

  //----------
  // set FoamHolder_box
  //----------
  sprintf(foamHolderBox->name,"FoamHolder");
  for (i=0;i<3;i++) foamHolderBox->halfXYZ[i]=foamHolder_halfXYZ[i];
  foamHolderBox->posXYZ=G4ThreeVector(0.0*cm,0.0*cm,foamHolder_posz);
  foamHolderBox->material=materialList->GetmRichMaterial("Air_Opt");
  foamHolderBox->sensitivity=0;

  foamHolderBox->color=G4Colour(0.2,0.498,0.369);
  foamHolderBox->visibility=true;
  foamHolderBox->wireframe=true;
  foamHolderBox->surface=false;

  //----------
  // set FoamHolder_polyhedra
  //----------
  sprintf(foamHolderPoly->name,"FoamHolder");
  foamHolderPoly->pos=G4ThreeVector(0,0,0);
  foamHolderPoly->start=45.0*myPI/180.0;
  foamHolderPoly->theta=2*myPI;
  foamHolderPoly->numSide=4;
  foamHolderPoly->num_zLayer=2;

  foamHolderPoly->z[0]=agel_posz-agel_halfXYZ[2];   //front of agel
  foamHolderPoly->z[1]=agel_posz+agel_halfXYZ[2];   //back of sensor plane

  foamHolderPoly->rinner[0]=agel_halfXYZ[0];
  foamHolderPoly->rinner[1]=agel_halfXYZ[0];

  foamHolderPoly->router[0]=foamHolderPoly->rinner[0]+foamHolderThicknessXYZ[0];
  foamHolderPoly->router[1]=foamHolderPoly->rinner[1]+foamHolderThicknessXYZ[0];

  foamHolderPoly->material=materialList->GetmRichMaterial("Air_Opt");
  foamHolderPoly->sensitivity=0;

  foamHolderPoly->color=G4Colour(0.298,0.6,0.471);
  foamHolderPoly->visibility=true;
  foamHolderPoly->wireframe=true;
  foamHolderPoly->surface=false;

  //----------
  // set aerogel
  //----------
  sprintf(aerogel->name,"Aerogel");
  for (i=0;i<3;i++) aerogel->halfXYZ[i]=agel_halfXYZ[i];
  aerogel->posXYZ=G4ThreeVector(0,0,agel_posz);
  aerogel->material=materialList->GetmRichMaterial("Aerogel2");
  //aerogel->material=Air_Opt;
  aerogel->sensitivity=0;

  aerogel->color=G4Colour(1.0,0.65,0.0);
  aerogel->visibility=true;
  aerogel->wireframe=true;
  aerogel->surface=true;

  //----------
  // set Fresnel lens
  //----------
  sprintf(fresnelLens->name,"FresnelLens");
  fresnelLens->pos=G4ThreeVector(0,0,lens_z);
  fresnelLens->material=G4Material::GetMaterial("Acrylic");
  fresnelLens->sensitivity=0;
  fresnelLens->color=G4Colour(0.0,1.0,1.0);
  fresnelLens->visibility=true;
  fresnelLens->wireframe=true;
  fresnelLens->surface=false;

  //----------
  // set mirror
  //----------
  sprintf(mirror->name,"mirror");
  mirror->pos=G4ThreeVector(0,0,0);
  mirror->start=45.0*myPI/180.0;
  mirror->theta=2*myPI;
  mirror->numSide=4;
  mirror->num_zLayer=2;

  mirror->z[0]=lens_z+fresnelLens->halfXYZ[2]+lens_gap;   //back of lens+air gap
  mirror->z[1]=glassWindow_z-glassWindow_halfXYZ[2];   //front of sensor plane

  mirror->rinner[0]=agel_halfXYZ[0];
  mirror->rinner[1]=sensor_total_halfx;

  mirror->router[0]=mirror->rinner[0]+mirrorThickness;
  mirror->router[1]=mirror->rinner[1]+mirrorThickness;

  mirror->material=materialList->GetmRichMaterial("Aluminum");
  mirror->sensitivity=0;

  mirror->color=G4Colour(1.0,1.0,0.0);
  mirror->visibility=true;
  mirror->wireframe=true;

  //----------
  // set glass window
  //----------
  sprintf(glassWindow->name, "glassWindow");
  for (i=0;i<3;i++) glassWindow->halfXYZ[i]=glassWindow_halfXYZ[i];
  glassWindow->posXYZ=G4ThreeVector(glassWindow_halfXYZ[0]+sensorGap,  //the position of the first sensor module.  
                                    glassWindow_halfXYZ[0]+sensorGap,  //the position of other sensor module will
                                    glassWindow_z);                    //be set glass window position in another func.
  glassWindow->material=materialList->GetmRichMaterial("Borosilicate");
  glassWindow->sensitivity=0;

  glassWindow->color=G4Colour(0.101,0.737,0.612);
  glassWindow->visibility=true;
  glassWindow->wireframe=true;
  glassWindow->surface=false;

  //----------
  // set sensor
  //----------
  sprintf(sensor->name, "sensor");
  for (i=0;i<3;i++) sensor->halfXYZ[i]=phodet_halfXYZ[i];
  sensor->posXYZ=G4ThreeVector(0,0,phodet_z);      //temporary. will be set in another func.
  sensor->material=materialList->GetmRichMaterial("Air_Opt");
  sensor->sensitivity=0;

  sensor->color=G4Colour(0.0,0.0,0.63);
  sensor->visibility=true;
  sensor->wireframe=true;
  sensor->surface=true;

  //----------
  // set readout
  //----------
  sprintf(readout->name,"readout");
  readout->pos=G4ThreeVector(0,0,0);
  readout->start=45.0*myPI/180.0;
  readout->theta=2*myPI;
  readout->numSide=4;
  readout->num_zLayer=2;

  readout->z[0]=readout_z[0];
  readout->z[1]=readout_z[1];

  readout->rinner[0]=sensor_total_halfx;
  readout->rinner[1]=readout->rinner[0];

  readout->router[0]=readout->rinner[0]+readoutThickness;
  readout->router[1]=readout->router[0];

  readout->material=materialList->GetmRichMaterial("Aluminum");
  readout->sensitivity=0;

  readout->color=G4Colour(1.0,0.0,0.0);
  readout->visibility=true;
  readout->wireframe=true;
  readout->surface=false;

}
//________________________________________________________________________//
mRichParameter::~mRichParameter(){;}
//________________________________________________________________________//
void mRichParameter::SetPar_glassWindow(G4double x, G4double y)
{
  glassWindow->posXYZ.setX(x);
  glassWindow->posXYZ.setY(y);
}
//________________________________________________________________________//
void mRichParameter::SetPar_sensor(G4double x, G4double y)
{
  sensor->posXYZ.setX(x);
  sensor->posXYZ.setY(y);
}
//________________________________________________________________________//
BoxPar* mRichParameter::GetBoxPar(const char* componentName)
{
  if (strcmp(componentName,"holderBox")==0) return holderBox;
  else if (strcmp(componentName,"hollowVolume")==0) return hollowVolume;
  else if (strcmp(componentName,"foamHolderBox")==0) return foamHolderBox;
  else if (strcmp(componentName,"aerogel")==0) return aerogel;
  else if (strcmp(componentName,"glassWindow")==0) return glassWindow;
  else if (strcmp(componentName,"sensor")==0) return sensor;
  else printf("mRichParameter::GetBoxPar() ----- ERROR: cannot find parameter=%s\n",componentName);

  //return 0;
  return aerogel;
}
//________________________________________________________________________//
LensPar* mRichParameter::GetLensPar(const char* componentName)
{
  if (strcmp(componentName,"fresnelLens")==0) return fresnelLens;
  else printf("mRichParameter::GetLensPar() ----- ERROR: cannot find parameter=%s\n",componentName);

  return 0;
}
//________________________________________________________________________//
PolyPar* mRichParameter::GetPolyPar(const char* componentName)
{
  if (strcmp(componentName,"foamHolderPoly")==0) return foamHolderPoly;
  else if (strcmp(componentName,"mirror")==0) return mirror;
  else if (strcmp(componentName,"readout")==0) return readout;
  else printf("mRichParameter::GetPolyPar() ----- ERROR: cannot find parameter=%s\n",componentName);

  return 0;
}
//________________________________________________________________________//
G4VPhysicalVolume* PHG4mRICHDetector::build_box(BoxPar* par, G4LogicalVolume* motherLV)
{
  G4Box* box = new G4Box(par->name,par->halfXYZ[0],par->halfXYZ[1],par->halfXYZ[2]);
  G4LogicalVolume* log = new G4LogicalVolume(box,par->material,par->name,0,0,0);
  G4VPhysicalVolume* phy=new G4PVPlacement(0,par->posXYZ,log,par->name,motherLV,false,0);

  G4VisAttributes* visAtt = new G4VisAttributes(par->color);
  visAtt->SetVisibility(par->visibility);
  visAtt->SetForceWireframe(par->wireframe);
  visAtt->SetForceSolid(par->surface);
  log->SetVisAttributes(visAtt);

  return phy;
}
//________________________________________________________________________//
G4VPhysicalVolume* PHG4mRICHDetector::build_polyhedra(PolyPar* par, G4LogicalVolume* motherLV)
{
  G4Polyhedra* polyhedra=new G4Polyhedra(par->name, par->start,par->theta, par->numSide,
                                         par->num_zLayer,par->z, par->rinner, par->router);
  G4LogicalVolume* log = new G4LogicalVolume(polyhedra,par->material,par->name,0,0,0);
  G4VPhysicalVolume* phy = new G4PVPlacement(0,par->pos,log,par->name,motherLV,false,0);

  G4VisAttributes* visAtt = new G4VisAttributes(par->color);
  visAtt->SetVisibility(par->visibility);
  visAtt->SetForceWireframe(par->wireframe);
  visAtt->SetForceSolid(par->surface);
  log->SetVisAttributes(visAtt);

  return phy;
}
//________________________________________________________________________//
G4VPhysicalVolume* PHG4mRICHDetector::build_holderBox(mRichParameter* detectorParameter, G4LogicalVolume* motherLV)
{
  G4VPhysicalVolume* holderBox=build_box(detectorParameter->GetBoxPar("holderBox"),motherLV);
  return build_box(detectorParameter->GetBoxPar("hollowVolume"),holderBox->GetLogicalVolume());
}
//________________________________________________________________________//
void PHG4mRICHDetector::build_foamHolder(mRichParameter* detectorParameter,G4LogicalVolume* motherLV)
{
  build_box(detectorParameter->GetBoxPar("foamHolderBox"),motherLV);
  build_polyhedra(detectorParameter->GetPolyPar("foamHolderPoly"),motherLV);
}
//________________________________________________________________________//
void PHG4mRICHDetector::build_aerogel(mRichParameter* detectorParameter,G4VPhysicalVolume* motherPV)
{
  G4VPhysicalVolume* aerogel=build_box(detectorParameter->GetBoxPar("aerogel"),motherPV->GetLogicalVolume());

  const G4int num = 2;
  G4double Ephoton[num] = {2.034*eV, 4.136*eV};

  G4OpticalSurface* OpWaterSurface = new G4OpticalSurface("WaterSurface");
  OpWaterSurface->SetType(dielectric_dielectric);
  OpWaterSurface->SetFinish(ground);
  OpWaterSurface->SetModel(unified);
  new G4LogicalBorderSurface("WaterSurface",aerogel,motherPV,OpWaterSurface);

  G4double RefractiveIndex[num] = {1.35, 1.40};
  G4double SpecularLobe[num]    = {0.3, 0.3};
  G4double SpecularSpike[num]   = {0.2, 0.2};
  G4double Backscatter[num]     = {0.2, 0.2};

  G4MaterialPropertiesTable* myST1 = new G4MaterialPropertiesTable();
  myST1->AddProperty("RINDEX",                Ephoton, RefractiveIndex, num);
  myST1->AddProperty("SPECULARLOBECONSTANT",  Ephoton, SpecularLobe,    num);
  myST1->AddProperty("SPECULARSPIKECONSTANT", Ephoton, SpecularSpike,   num);
  myST1->AddProperty("BACKSCATTERCONSTANT",   Ephoton, Backscatter,     num);

  OpWaterSurface->SetMaterialPropertiesTable(myST1);
}
//________________________________________________________________________//
void PHG4mRICHDetector::build_mirror(mRichParameter* detectorParameter,G4VPhysicalVolume* motherPV)
{
  G4VPhysicalVolume* mirror=build_polyhedra(detectorParameter->GetPolyPar("mirror"),motherPV->GetLogicalVolume());

  //-----------
  //   Optical properties of the interface between the Air and Reflective Surface 
  //   For Mirror, reflectivity is set at 95% and specular reflection is assumed.
  //-----------
  G4OpticalSurface *OpticalAirMirror = new G4OpticalSurface("AirMirrorSurface");
  OpticalAirMirror->SetModel(unified);
  OpticalAirMirror->SetType(dielectric_dielectric);
  OpticalAirMirror->SetFinish(polishedfrontpainted);

  const G4int NUM = 2;
  G4double lambda_min = 200*nm ;
  G4double lambda_max = 700*nm ;

  G4double XX[NUM] = {h_Planck*c_light/lambda_max, h_Planck*c_light/lambda_min} ;
  G4double ICEREFLECTIVITY[NUM]      = { 0.95, 0.95 };

  G4MaterialPropertiesTable *AirMirrorMPT = new G4MaterialPropertiesTable();
  AirMirrorMPT->AddProperty("REFLECTIVITY", XX, ICEREFLECTIVITY,NUM);
  OpticalAirMirror->SetMaterialPropertiesTable(AirMirrorMPT);

  new G4LogicalBorderSurface("Air/Mirror Surface",motherPV,mirror,OpticalAirMirror);
}
//________________________________________________________________________//
void PHG4mRICHDetector::build_sensor(mRichParameter* detectorParameter,G4LogicalVolume* motherLV)
{
  //position of the first sensor module
  G4double last_x=detectorParameter->GetBoxPar("glassWindow")->posXYZ.getX();
  G4double last_y=last_x;

  G4double x,y;

  int i;
  for (i=0;i<4;i++) {
    if (i==0) {
      x=last_x;
      y=last_y;
    }
    else {
      x=-last_y;
      y=last_x;
    }

    detectorParameter->SetPar_glassWindow(x,y);
    build_box(detectorParameter->GetBoxPar("glassWindow"),motherLV);

    detectorParameter->SetPar_sensor(x,y);
    build_box(detectorParameter->GetBoxPar("sensor"),motherLV);

    last_x=x;
    last_y=y;
  }//end of for(i)
}
//________________________________________________________________________//
void PHG4mRICHDetector::build_lens(LensPar* par, G4LogicalVolume* motherLV)
{
  const G4int NumberOfGrooves=floor((par->eff_diameter/2.0)/par->grooveWidth);
  G4Polycone *Groove_poly[NumberOfGrooves];
  G4LogicalVolume *Groove_log[NumberOfGrooves];
  
  G4VisAttributes* SurfaceVisAtt = new G4VisAttributes(G4Colour(0.0,1.0,1.0));
  SurfaceVisAtt->SetVisibility(true);                                                                            
  SurfaceVisAtt->SetForceWireframe(true);
  SurfaceVisAtt->SetForceSolid(true);

  int igroove;
  for (igroove=0;igroove<1000;igroove++) {      //just put a arbitrary large number
    //--------------------------------
    //Grooves' inner and outer radius
    //--------------------------------
    G4double iRmin1 = (igroove+0)*par->grooveWidth;
    G4double iRmax1 = (igroove+1)*par->grooveWidth;
    G4double iRmin2 = iRmin1;
    G4double iRmax2 = iRmin2+0.0001;

    G4double lens_poly_rmin[3] = {iRmin1, iRmin1, iRmin2};
    G4double lens_poly_rmax[3] = {iRmax1, iRmax1, iRmax2};

    if (iRmax1>par->diameter/2.0) break;     //if iRmax1>Lens radius (outside the lens), break

    //--------------------------------
    //phi angle
    //--------------------------------
    G4double phi1;
    G4double phi2;
    G4double deltaPhi;

    if (iRmax1< par->halfXYZ[0]) {   //draw a full circle
      phi1=0;                        //in rad
      deltaPhi=twopi;                //in rad. two pi
    }
    else {
      phi1=acos(par->halfXYZ[0]/iRmax1);   //in rad
      phi2=asin(par->halfXYZ[0]/iRmax1);   //in rad, assume lens is square -> halfy=halfx
      deltaPhi=phi2-phi1;
    }
    //--------------------------------
    //grooves profile
    //--------------------------------
    G4double lens_poly_z[3];
    int numOfLayer;

    if (iRmin1< par->eff_diameter/2.0) {   //if iRmin>=effective radius, dZ=0, i.e. flat
      numOfLayer=3;
      G4double dZ = par->GetSagita(iRmax1) - par->GetSagita(iRmin1);
      lens_poly_z[0]=par->halfXYZ[2];
      lens_poly_z[1]=-par->halfXYZ[2]+dZ;
      lens_poly_z[2]=-par->halfXYZ[2];
    }
    else {
      numOfLayer=2;
      lens_poly_z[0] = par->halfXYZ[2];
      lens_poly_z[1] = par->halfXYZ[2]-par->centerThickness;
      lens_poly_z[2] = 0;
    }

    //--------------------------------
    //build grooves
    //--------------------------------
    int repeat=1;
    if (iRmax1>= par->halfXYZ[0]) { repeat=4; }   //4 edges
    for (int i=0;i<repeat;i++) {
      Groove_poly[i]=new G4Polycone(par->name,phi1,deltaPhi, numOfLayer, lens_poly_z, lens_poly_rmin, lens_poly_rmax);
      Groove_log[i]=new G4LogicalVolume(Groove_poly[i],par->material,par->name,0,0,0);
      new G4PVPlacement(0,par->pos,Groove_log[i],par->name,motherLV,false,0);

      Groove_log[i]->SetVisAttributes(SurfaceVisAtt);
      phi1=phi1+halfpi;   //g4 pre-defined: halfpi=pi/2
    }
  }
  
}
//________________________________________________________________________//
G4LogicalVolume* PHG4mRICHDetector::build_Space(G4LogicalVolume* logicWorld, G4double (&bowlPar)[4])
{
  int i;

  G4Material* Air = G4Material::GetMaterial("G4_AIR");
  G4double r_inner=3*m;
  G4double r_outer=3.7*m;
  G4double phi[2];          //0=min, 1=max
  G4double hin[2],hout[2];    //0=min, 1=max
  G4double cone_l[3];         //0=l_in, 1=l_out, 2=l_out-l_in

  phi[0]=eta2polarAngle(1.9);
  phi[1]=eta2polarAngle(1.1);

  cone_l[0]=0.1*r_outer;
  cone_l[1]=r_outer;
  cone_l[2]=cone_l[1]-cone_l[0];  //=0.9*r_outer

  for (i=0;i<2;i++) {
    hin[i]=cone_l[0]*tan(phi[i]);
    hout[i]=cone_l[1]*tan(phi[i]);
  }
  
  //copy dimension back to Construct()
  bowlPar[0]=r_inner;
  bowlPar[1]=r_outer;
  bowlPar[2]=phi[0];
  bowlPar[3]=phi[1];

  G4VSolid* sphere =new G4Sphere("sphere",
                                       r_inner,r_outer,
                                       0,twopi,
                                       0,twopi);
  
  G4VSolid* cone=new G4Cons("cone",
			    hin[0], hin[1],
			    hout[0],hout[1],
			    cone_l[2]/2.0,0,twopi);

  G4IntersectionSolid* bowl_solid=new G4IntersectionSolid("bowl",sphere,cone,
							  0,G4ThreeVector(0.00*m, 0.00*m, r_outer-cone_l[2]/2.0));
  
  G4LogicalVolume* bowl_log =  new G4LogicalVolume(bowl_solid, Air, "bowl", 0, 0, 0);
  
  G4VisAttributes* visAtt = new G4VisAttributes();
  visAtt->SetVisibility(true);
  visAtt->SetForceSolid(true);
  visAtt->SetColour(G4Color(1.0,0.0,1.0,0.1));    //R,G,B,alpha
  bowl_log->SetVisAttributes(visAtt);
  
  new G4PVPlacement( 0,G4ThreeVector(0, 0, 0),
                     bowl_log, "bowl", logicWorld, 0, false, 1);
  
  return bowl_log;
}
//________________________________________________________________________//
void PHG4mRICHDetector::build_mRICH_wall(G4LogicalVolume*space, G4LogicalVolume*a_mRICH, G4double* bowlPar)
{
  int i,j;
  
  // mRICH half width, height, and length + air gap
  G4double gap=3*cm;         //large gap to avoid overlap. temporary solution
  G4Box* mRICH_box=dynamic_cast<G4Box*>(a_mRICH->GetSolid());
  G4double halfWidth=mRICH_box->GetXHalfLength();// + gap;
  G4double halfHeight=halfWidth;
  G4double halfLength=mRICH_box->GetZHalfLength()+gap;

  // space (bowl shape) dimension
  G4double rinner=bowlPar[0];
  G4double router=bowlPar[1];
  G4double phi_min=bowlPar[2];
  G4double phi_max=bowlPar[3];

  if ((router-rinner)<(2*halfLength)) {
    G4cout<<"PHG4mRICHDetector::build_mRICH_wall ::::::::::: ERROR! Not enough space for mRICH"<<G4endl;
    return;
  }

  //position and name of each copy of mRICH
  G4double x,y,z;
  char name[50];

  G4double c;  // circumference of each ring on the bowl
  int N;       // num. of mRICH on each ring on the bowl

  G4double theta,phi,deltaTheta, deltaPhi;
  deltaPhi=2*atan(halfHeight/rinner);//+pi/180;
  phi=0;

  int count=0;
  i=0;
  for (phi=phi_min+deltaPhi;phi<phi_max-deltaPhi;phi=phi+deltaPhi) {  //add extra space to avoid overlap
    c=twopi*(rinner)*sin(phi);
    N=floor(c/(2*sqrt(2)*halfWidth));
    deltaTheta=2*pi/N;
    
    if (1) {
      printf("-----------------------------------------------------------------------------\n");
      printf("ring %d",i+1);
      printf("(minPhi, maxPhi, delPhi)=(%.2lf,%.2lf,%.2f) deg ",phi_min*180/pi,phi_max*180/pi,deltaPhi*180/pi);
      printf("phi=%.2f deg \n",phi*180/pi);
    }
  
    j=0;
    for (theta=0;theta<=twopi;theta=theta+deltaTheta) {
      x=(rinner+halfLength)*cos(theta)*sin(phi);
      y=(rinner+halfLength)*sin(theta)*sin(phi);
      z=(rinner+halfLength)*cos(phi);
      
      if (0) printf("x,y,z= %.2lf, %0.2lf, %0.2lf\n",x,y,z);
      G4RotationMatrix rot= G4RotationMatrix(); 
      rot.rotateX(phi*(-1)*sin(theta)*180*deg/pi);
      rot.rotateY(phi*cos(theta)*180*deg/pi);
       
      sprintf(name,"mRICH_%d_%d",i+1,j+1);
      
      new G4PVPlacement(G4Transform3D(rot, G4ThreeVector(x, y, z)),
			a_mRICH,name,space,
			0, 0, 1);    //last digit for checking overlapping 
      
      j++;
      count++;
    }
    i++;
    //if (i>0) break;
  }

  printf("-----------------------------------------------------------------------------\n");
  printf("%d detectors are built\n",count);
}
//________________________________________________________________________//
G4double PHG4mRICHDetector::eta2polarAngle(G4double eta)
{
  return 2*atan(exp(-eta));

}
