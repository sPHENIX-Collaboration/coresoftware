// $Id: $                                                                                             

/*!
 * \file Fun4All_ImportGeom.C
 * \brief Example Fun4All macro part to import an external Geometry to DST node
 * \author Jin Huang <jhuang@bnl.gov>
 * \version $Revision:   $
 * \date $Date: $
 */

void
Fun4All_ImportGeom()
{
  gSystem->Load("libphgeom.so");

  Fun4AllServer *se = Fun4AllServer::instance();
  se->Verbosity(1);
  // just if we set some flags somewhere in this macro
  recoConsts *rc = recoConsts::instance();

  PHGeomFileImport * import = new PHGeomFileImport("./sPHENIX.root");
  se->registerSubsystem(import);

  // dummy input
  Fun4AllInputManager *in = new Fun4AllDummyInputManager("JADE");
  se->registerInputManager(in);

  // output in DST
  Fun4AllDstOutputManager *out = new Fun4AllDstOutputManager("DSTOUT",
      "DST.root");
  se->registerOutputManager(out);

  // run one event as example
  se->run(1);

  se->End();
  std::cout << "All done" << std::endl;
  delete se;
  gSystem->Exit(0);

}

