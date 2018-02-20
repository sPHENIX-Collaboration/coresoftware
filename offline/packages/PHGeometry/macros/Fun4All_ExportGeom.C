// $Id: $                                                                                             
 
/*!
 * \file Fun4All_ExportGeom.C
 * \brief get geometry from DST file and output TGeoManager ROOT files or other formats of geometry files
 * \author Jin Huang <jhuang@bnl.gov>
 * \version $Revision:   $
 * \date $Date: $
 */

//! \brief get geometry from DST file and output TGeoManager ROOT files or other formats of geometry files
void
Fun4All_ExportGeom(string DST_file_name = "sPHENIX.root_DST.root")
{

  gSystem->Load("libphgeom.so");

  // in case DST contains sPHENIX stuff
  gSystem->Load("libg4calo.so");
  gSystem->Load("libg4vertex.so");
  gSystem->Load("libg4eval.so");

  Fun4AllServer *se = Fun4AllServer::instance();
  se->Verbosity(1);
  recoConsts *rc = recoConsts::instance();
   rc->set_IntFlag("RUNNUMBER", 12345);

  Fun4AllInputManager *hitsin = new Fun4AllDstInputManager("DSTin");
  hitsin->fileopen(DST_file_name);
  se->registerInputManager(hitsin);

  // run one event as example
  se->run(1);

  string output = DST_file_name + "_export.root";
  PHGeomUtility::ExportGeomtry(se->topNode(),DST_file_name + "_export.root");
  cout <<"Done export Geometry to "<<output<<endl;

  se->End();
  delete se;
  cout <<"All done"<<endl;
  gSystem->Exit(0);

}



