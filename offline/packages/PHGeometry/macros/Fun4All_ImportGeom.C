// $Id: $                                                                                             

/*!
 * \file Fun4All_ImportGeom.C
 * \brief Example Fun4All macro part to import an external Geometry to DST node
 * \author Jin Huang <jhuang@bnl.gov>
 * \version $Revision:   $
 * \date $Date: $
 */

#include <string>
using namespace std;

//! Read in a Geometry file, and output DST and ROOT TGeo files
void
Fun4All_ImportGeom(const string & geom_file = "./sPHENIX.root")
{
  gSystem->Load("libphgeom.so");

  Fun4AllServer *se = Fun4AllServer::instance();
  se->Verbosity(1);
  // just if we set some flags somewhere in this macro
  recoConsts *rc = recoConsts::instance();

  PHGeomFileImport * import = new PHGeomFileImport(geom_file);
  se->registerSubsystem(import);

  // dummy input
  Fun4AllInputManager *in = new Fun4AllDummyInputManager("JADE");
  se->registerInputManager(in);

  // output in DST
  Fun4AllDstOutputManager *out = new Fun4AllDstOutputManager("DSTOUT",
      geom_file + "_DST.root");
  se->registerOutputManager(out);

  // run one event as example
  se->run(1);

  PHGeomUtility::ExportGeomtry(se->topNode(),geom_file + "_export.root");

  se->End();
  std::cout << "All done" << std::endl;
  delete se;
  gSystem->Exit(0);

}

