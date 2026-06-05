#include <QA.C>

#include <inttcalib/InttCalib.h>

#include <ffamodules/HeadReco.h>
#include <ffamodules/FlagHandler.h>
#include <ffamodules/SyncReco.h>
#include <ffamodules/CDBInterface.h>

#include <ffarawmodules/InttCheck.h>
#include <ffarawmodules/StreamingCheck.h>
#include <ffarawmodules/TpcCheck.h>

#include <fun4all/Fun4AllDstOutputManager.h>
#include <fun4all/Fun4AllInputManager.h>
#include <fun4all/Fun4AllDstInputManager.h>
#include <fun4all/Fun4AllOutputManager.h>
#include <fun4all/Fun4AllServer.h>
#include <fun4all/Fun4AllUtils.h>

#include <fun4allraw/Fun4AllStreamingInputManager.h>
#include <fun4allraw/InputManagerType.h>
#include <fun4allraw/SingleGl1PoolInput.h>
#include <fun4allraw/SingleInttPoolInput.h>
#include <fun4allraw/SingleMicromegasPoolInput.h>
#include <fun4allraw/SingleMvtxPoolInput.h>
#include <fun4allraw/SingleTpcPoolInput.h>
#include <fun4allraw/SingleTpcTimeFrameInput.h>

#include <phool/recoConsts.h>
#include <GlobalVariables.C>

R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libfun4allutils.so)
R__LOAD_LIBRARY(libffamodules.so)
R__LOAD_LIBRARY(libfun4allraw.so)
R__LOAD_LIBRARY(libffarawmodules.so)
R__LOAD_LIBRARY(libcdbobjects.so)
R__LOAD_LIBRARY(libffamodules.so)
R__LOAD_LIBRARY(libPHGarfield.so)

#define Nebdc 24
#define Nserver 2

// Global namespace to assist drawing...
TPolyLine3D *npoly3[48];
TPolyLine3D *spoly3[48];
TPolyLine   *npoly2[48];
TPolyLine   *spoly2[48];
TGeoTube    *tubby;
TCanvas     *canny;
TCanvas     *canny2;

TBox        *boxer1;
TBox        *boxer2;

void TestFieldMap()
{
  recoConsts* rc = recoConsts::instance();

  rc->set_StringFlag("CDB_GLOBALTAG","FieldMapTest");
  rc->set_uint64Flag("TIMESTAMP",1);

  auto cdb = CDBInterface::instance();
  std::string url = cdb->getUrl("FIELDMAP_TRACKING");
  std::cout << "Field map URL:\n" << url << std::endl;

  Fun4AllServer *se = Fun4AllServer::instance();

  Enable::QA  = false;
  Enable::CDB = true;
  
  // Register a whole slew of input managers...
  // NOTE:  This depends upon the requested files being in frog.  Ribbit.
  char nextinput[500];
  char nextfile[500];
  Fun4AllInputManager* in[Nebdc]; 
  for (unsigned int ebdc=0; ebdc<24; ebdc++)
    {
      for (unsigned int server=0; server<2; server++)
	{
	  sprintf(nextinput,"ebdc%02d_%01d",ebdc,server);
	  //sprintf(nextfile,"DST_STREAMING_EVENT_ebdc%02d_%01d_run3line_laser_ana540_nocdbtag_v001-00064890-00000.root",ebdc,server);  // Line Laser
	  sprintf(nextfile,"DST_STREAMING_EVENT_ebdc%02d_%01d_run3auau_ana514_nocdbtag_v001-00075570-00000.root",ebdc,server);          // AuAu Zero Field
	  std::cout << nextfile << " " << nextinput << endl;
	  in[ebdc] = new Fun4AllDstInputManager(nextinput);
	  in[ebdc]->fileopen(nextfile);
	  se->registerInputManager(in[ebdc]);
	}
    }

  // Now register a flag handler because MAAABE it will make the CDB work correctly?
  //SubsysReco *fh = new FlagHandler();
  //se->registerSubsystem(fh);
  
  // Register my analysis module.
  PHGarfield *phg = new PHGarfield();
  se->registerSubsystem(phg);

  se->run(4);

  canny  = new TCanvas("canny","canny",3000,2500);
  canny2 = new TCanvas("canny2","canny2",3000,2500);
  tubby  = new TGeoTube("tubby",20,80,110);

  canny->cd();
  tubby->Draw();

  canny2->cd();
  //gPad->DrawFrame(-150., -10., 150., 10.);
  gPad->DrawFrame(-150., -100., 150., 100.);
  boxer1 = new TBox(-102,20,102,78);
  boxer1->Draw();
  boxer2 = new TBox(-102,-78,102,-20);
  boxer2->Draw("same");
  
  for (int i=0; i<48; i++)
    {
      canny->cd();
      npoly3[i] = phg->ReverseDrift(0,phg->radii[i],102);
      npoly3[i]->SetLineColor(kRed);
      npoly3[i]->SetLineWidth(3);
      npoly3[i]->Draw("same");

      canny2->cd();
      int N = npoly3[i]->GetN();
      float *p =  npoly3[i]->GetP();
      float x[500];
      float y[500];
      float z[500];
      for (int j=0; j<N; j++)
      {
	x[j] = p[3*j+0];
	y[j] = p[3*j+1];
	z[j] = p[3*j+2];
	/*
	std::cout << j
            << "  " << x[j]
            << "  " << y[j]
            << "  " << z[j]
            << std::endl;
	*/
      }
      npoly2[i] = new TPolyLine(N,z,y);
      npoly2[i]->SetLineColor(kRed);
      npoly2[i]->SetLineWidth(3);
      npoly2[i]->Draw("Lsame");
    }

  for (int i=0; i<48; i++)
    {
      canny->cd();
      spoly3[i] = phg->ReverseDrift(0,phg->radii[i],-102);
      spoly3[i]->SetLineColor(kCyan);
      spoly3[i]->SetLineWidth(3);
      spoly3[i]->Draw("same");

      canny2->cd();
      int N = spoly3[i]->GetN();
      float *p =  spoly3[i]->GetP();
      float x[500];
      float y[500];
      float z[500];
      for (int j=0; j<N; j++)
      {
	x[j] = p[3*j+0];
	y[j] = p[3*j+1];
	z[j] = p[3*j+2];
	/*
	std::cout << j
            << "  " << x[j]
            << "  " << y[j]
            << "  " << z[j]
            << std::endl;
	*/
      }
      spoly2[i] = new TPolyLine(N,z,y);
      spoly2[i]->SetLineColor(kCyan);
      spoly2[i]->SetLineWidth(3);
      spoly2[i]->Draw("Lsame");
    }
  //se->dumpHistos("Looker.root");
}
