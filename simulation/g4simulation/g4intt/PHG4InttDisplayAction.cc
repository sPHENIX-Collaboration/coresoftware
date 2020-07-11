#include "PHG4InttDisplayAction.h"

#include <g4main/PHG4DisplayAction.h>  // for PHG4DisplayAction

#include <TSystem.h>

#include <Geant4/G4Colour.hh>          // for G4Colour
#include <Geant4/G4LogicalVolume.hh>
#include <Geant4/G4String.hh>          // for G4String
#include <Geant4/G4VisAttributes.hh>

#include <iostream>
#include <utility>                     // for pair

using namespace std;

PHG4InttDisplayAction::PHG4InttDisplayAction(const std::string &name)
  : PHG4DisplayAction(name)
{
}

PHG4InttDisplayAction::~PHG4InttDisplayAction()
{
  for (auto &it : m_VisAttVec)
  {
    delete it;
  }
  m_VisAttVec.clear();
}

void PHG4InttDisplayAction::ApplyDisplayAction(G4VPhysicalVolume *physvol)
{
  // check if vis attributes exist, if so someone else has set them and we do nothing

  G4Colour colour_air ( 0.0, 0.0, 0.0, 0.0 );
  G4Colour colour_CFRP( 0.4, 0.4, 0.4, 1 );

  for (auto it : m_LogicalVolumeMap)
  {
    G4LogicalVolume *logvol = it.first;
    if (logvol->GetVisAttributes())
    {
      continue;
    }

    G4VisAttributes *visatt = new G4VisAttributes();
    visatt->SetVisibility(true);
    visatt->SetForceSolid(true);
    m_VisAttVec.push_back(visatt);  // for later deletion

    if (it.second == "FPHX")
    {
      visatt->SetColour( G4Colour(1.0, 0.843, 0.0, 0.5) ); // HTML gold
      visatt->SetVisibility( true  );

    }
    else if (it.second == "FPHXGlue" || it.second == "SiGlue")
    {
      visatt->SetColour( G4Colour(0.1, 0.1, 0.1, 0.8) );
      visatt->SetVisibility( true  );

    }
    else if (it.second == "HdiCopper")
    {
      visatt->SetColour( G4Colour(0.7, 0.4, 0, 1) ); // copper color
      visatt->SetVisibility( true );

    }
    else if (it.second == "HdiKapton")
    {
      visatt->SetColour( G4Colour(0.0, 0.590, 1.0, 0.5 ) ); // blue
      visatt->SetVisibility( true );

    }
    else if (it.second == "Ladder" || it.second == "FPHXContainer" || it.second == "FPHXGlueContainer" || it.second == "StaveBox" )
    {
      visatt->SetColour( colour_air );
      //visatt->SetColour( G4Colour(1, 1, 1, 1 ) ); // black for debigging
      visatt->SetForceWireframe( true );
      visatt->SetVisibility( false );
      //visatt->SetVisibility( true );

    }
    else if (it.second == "Rail") // checked, ?
    {
      visatt->SetColour( G4Colour::Cyan());
      visatt->SetVisibility( false );

    }
    else if (it.second == "RohaCell")
    {
      visatt->SetColour( G4Colour(0.9, 0.9, 0.9, 0.5) ); // white
      visatt->SetVisibility( true );

    }
    else if (it.second == "SiActive")
    {
      visatt->SetColour( G4Colour(1.0, 0, 0.0, 0.5) ); // transparent red
      visatt->SetVisibility( true );

    }
    else if (it.second == "SiInActive")
    {
      visatt->SetColour( G4Colour(0, 0, 1, 0.5) ); // transparent blue
      visatt->SetVisibility( true );

    }
    else if (it.second == "StaveCooler")
    {
      visatt->SetColour( colour_CFRP );
      visatt->SetVisibility( true );

    }
    else if (it.second == "StaveCurve")
    {
      visatt->SetColour( colour_CFRP );
      visatt->SetVisibility( true );

    }
    else if (it.second == "StaveGlueBox")
    {
      visatt->SetColour(G4Colour::Cyan());
      visatt->SetVisibility( true );

    }
    else if (it.second == "StavePipe")
    {
      visatt->SetColour( colour_CFRP );
      visatt->SetVisibility( true );

    }
    else if (it.second == "StaveStraightInner")
    {
      visatt->SetColour(G4Colour::Grey());
      visatt->SetVisibility( true );
    }
    else if (it.second == "StaveStraightOuter")
    {
      visatt->SetColour( colour_CFRP );
      visatt->SetVisibility( true );

    }
    else if (it.second == "StaveWater")
    {
      visatt->SetColour(G4Colour::Blue());
      visatt->SetVisibility( true );
    }
    else if (it.second == "Endcap")
    {
      visatt->SetColour(G4Colour::Blue());
      visatt->SetVisibility( false );
    }
    else
    {
      cout << "did not assing color to " << it.first->GetName()
	   << " under " << it.second << endl;
      gSystem->Exit(1);
    }
    logvol->SetVisAttributes(visatt);
  }
  return;
}
