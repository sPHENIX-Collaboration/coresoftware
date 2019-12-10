void run_dump(const char *infile, const int evts=100)
{
  gSystem->Load("libg4hough.so");
  gSystem->Load("libphnodedump.so");
  gSystem->Load("libg4vertex_io.so");
  gSystem->Load("libg4bbc_io.so");
  
  Fun4AllServer* se = Fun4AllServer::instance();

  Dumper *dmp = new Dumper();
  gSystem->Exec("mkdir dump");
  dmp->SetOutDir("./dump");

  se->registerSubsystem(dmp);

  Fun4AllInputManager *in = new Fun4AllDstInputManager("DSTin");
  se->registerInputManager(in);
  se->fileopen("DSTin",infile);
  se->run(evts);
  se->End();
  delete se;
}
