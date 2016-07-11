// $Id: $                                                                                             

/*!
 * \file PHGeom_DSTInspection.C
 * \brief Quick inspection of PHGeoTGeo object in RUN/GEOMETRY node inside a DST file
 * \author Jin Huang <jhuang@bnl.gov>
 * \version $Revision:   $
 * \date $Date: $
 */

#include <cassert>

//! Quick inspection of PHGeoTGeo object in RUN/GEOMETRY node inside a DST file
//! Based on abhisek's display macro
void
PHGeom_DSTInspection(TString DST_file_name = "DST.root")
{
  TEveManager::Create();

  TFile * _file0 = TFile::Open("DST.root");
  assert(_file0->IsOpen());

  TTree * T1 = (TTree *) _file0->GetObjectChecked("T1", "TTree");
  assert(T1);

  TBranch * tgeom = T1->GetBranch("RUN.GEOMETRY._fGeom");
  assert(tgeom);

  tgeom->GetEntry(0);
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
  v->GetClipSet()->SetClipType(1);
  v->ColorSet().Background().SetColor(kMagenta + 4);
  v->SetGuideState(TGLUtil::kAxesEdge, kTRUE, kFALSE, 0);
  v->RefreshPadEditor(v);
  v->CurrentCamera().RotateRad(-0.5, 0.5);
  v->DoDraw();
}

