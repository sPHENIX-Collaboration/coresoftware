// $Id: $                                                                                             

/*!
 * \file PHGeom_DSTInspection.C
 * \brief Quick inspection of PHGeoTGeo object in RUN/GEOMETRY node inside a DST file
 * \author Jin Huang <jhuang@bnl.gov>
 * \version $Revision:   $
 * \date $Date: $
 */

#include <fun4all/SubsysReco.h>
#include <fun4all/Fun4AllServer.h>
#include <fun4all/Fun4AllInputManager.h>
#include <fun4all/Fun4AllDummyInputManager.h>
#include <fun4all/Fun4AllOutputManager.h>
#include <fun4all/Fun4AllDstInputManager.h>
#include <fun4all/Fun4AllNoSyncDstInputManager.h>
#include <fun4all/Fun4AllDstOutputManager.h>
#include <phool/recoConsts.h>
#include <phgeom/PHGeomUtility.h>

#include <TEveManager.h>
#include <TGLViewer.h>
#include <TGLUtil.h>
#include <TGLClip.h>
#include <TGeoManager.h>
#include <TROOT.h>
#include <TGeoManager.h>
#include <TEveGeoNode.h>

#include <cassert>
#include <string>

using namespace std;

R__LOAD_LIBRARY(libphgeom.so)
R__LOAD_LIBRARY(libg4dst.so)
R__LOAD_LIBRARY(libfun4all.so)


//! Quick inspection of PHGeoTGeo object in RUN/GEOMETRY node inside a DST file
//! Based on abhisek's display macro
void
PHGeom_DSTInspection(string DST_file_name = "G4sPHENIX.root",
    bool do_clip = true)
{
  TEveManager::Create();

  Fun4AllServer *se = Fun4AllServer::instance();
  se->Verbosity(1);
  recoConsts *rc = recoConsts::instance();
  rc->set_IntFlag("RUNNUMBER", 12345);

  Fun4AllInputManager *hitsin = new Fun4AllDstInputManager("DSTin");
  hitsin->fileopen(DST_file_name);
  se->registerInputManager(hitsin);

  // run one event as example
  se->run(1);

  PHGeomUtility::GetTGeoManager(se->topNode());

  assert(gGeoManager);

  if (!gROOT->GetListOfGeometries()->FindObject(gGeoManager))
    gROOT->GetListOfGeometries()->Add(gGeoManager);
  if (!gROOT->GetListOfBrowsables()->FindObject(gGeoManager))
    gROOT->GetListOfBrowsables()->Add(gGeoManager);
//  gGeoManager->UpdateElements();

  TGeoNode *current = gGeoManager->GetCurrentNode();
  //Alternate drawing
  //current->GetVolume()->Draw("ogl");
  //Print the list of daughters
  //current->PrintCandidates();
  for (int igeom = 0; igeom < current->GetNdaughters(); igeom++)
    {
      TGeoNode *geo_node = (TGeoNode*) current->GetNodes()->UncheckedAt(igeom);
      geo_node->GetVolume()->VisibleDaughters(kFALSE);
      geo_node->GetVolume()->SetTransparency(2);
      //Keep the pipe visible all the time
      if (string(geo_node->GetName()).find("PIPE") != string::npos)
        geo_node->GetVolume()->SetTransparency(0);
    }
  TEveGeoTopNode* eve_node = new TEveGeoTopNode(gGeoManager, current);
  eve_node->SetVisLevel(6);
  gEve->AddGlobalElement(eve_node);
  gEve->FullRedraw3D(kTRUE);

  // EClipType not exported to CINT (see TGLUtil.h):
  // 0 - no clip, 1 - clip plane, 2 - clip box
  TGLViewer *v = gEve->GetDefaultGLViewer();
  if (do_clip)
    {
      v->GetClipSet()->SetClipType( TGLClip::kClipPlane  );
    }
//  v->ColorSet().Background().SetColor(kMagenta + 4);
  v->SetGuideState(TGLUtil::kAxesEdge, kTRUE, kFALSE, 0);
  v->RefreshPadEditor(v);

  v->CurrentCamera().RotateRad(-1.6,0.);
  v->DoDraw();

}

