#ifndef __PMONDISPLAY__
#define __PMONDISPLAY__

#include <string.h>
#include <TPaveLabel.h>
#include <TCanvas.h>

class pmondisplay {   //++CINT

private:
  TCanvas *Pmon;
  TPaveLabel *pltitle; 
  TPaveLabel *pldstatus;
  TPaveLabel *plstatus;
  TPaveLabel *pldevt;
  TPaveLabel *plevt;
  TPaveLabel *pldstream;
  TPaveLabel *plstream;

  void  pmonsetup(const char *objname, const char * title)
    {

      Pmon = new TCanvas(objname, title,-158,36,401,131);
      Pmon->SetFillColor(46);
      Pmon->SetHighLightColor(2);
      Pmon->Range(0,0,1,1);
      Pmon->SetBorderSize(2);
      
      pltitle = new TPaveLabel(0.01 ,0.625 ,0.99 ,1.,title," ");
      pltitle->SetFillColor(42);
      pltitle->SetTextFont(42);
      pltitle->SetTextSize(0.65);
      


      // line 2 Status: xxxxxxx

      pldstatus = new TPaveLabel(0.01 ,0.415 ,0.495 ,0.595,"Status:"," ");
      pldstatus->SetFillColor(42);
      pldstatus->SetTextFont(22);
      pldstatus->SetTextSize(0.5);
      

      plstatus = new TPaveLabel(0.505 ,0.415 ,0.99 ,0.595,"unknown"," ");
      plstatus->SetFillColor(41);
      plstatus->SetTextFont(22);
      plstatus->SetTextSize(0.5);
      



      // line 3 Event: 78987


      pldevt = new TPaveLabel(0.01 ,0.215 ,0.495 ,0.395 ,"Event nr:"," ");
      pldevt->SetFillColor(42);
      pldevt->SetTextFont(22);
      pldevt->SetTextSize(0.5);
      //   plevt->Draw();
      

      plevt = new TPaveLabel(0.505 ,0.215 ,0.99 ,0.395 ,"0"," ");
      plevt->SetFillColor(41);
      plevt->SetTextFont(22);
      plevt->SetTextSize(0.5);
      //   plevt->Draw();
      

      // line 4 Stream: ---/...?.


      pldstream = new TPaveLabel(0.01 ,0.01 ,0.495 ,0.185 ,"Stream: "," ");
      pldstream->SetFillColor(42);
      pldstream->SetTextFont(22);
      pldstream->SetTextSize(0.5);
 
      plstream = new TPaveLabel(0.505 ,0.01 ,0.99 ,0.185 ," "," ");
      plstream->SetFillColor(41);
      plstream->SetTextFont(22);
      plstream->SetTextSize(0.5);
      //  plstream->Draw();
      Pmon->Modified();
      Pmon->Draw();
      Update();
    }


public:


  pmondisplay()
    {
      pmonsetup ("Pmon","Phenix Monitor Status");
    }

  pmondisplay(const char *objname)
    {
      pmonsetup (objname,"Phenix Monitor Status");
    }

  pmondisplay(const char *objname, const char * title)
    {
      pmonsetup (objname,title);
    }

  ~pmondisplay()
    {

      delete pldstream;
      delete pldevt;
      delete pldstatus;

      delete plstream;
      delete plevt;
      delete plstatus;
      delete pltitle;

      delete Pmon;
    }

  void Update()
    {
      pldstatus->Draw();
      pldevt->Draw();
      pldstream->Draw();

      pltitle->Draw();
      plstatus->Draw();
      plevt->Draw();
      plstream->Draw();

      Pmon->Modified();
      Pmon->Update();
      
      //    c1->cd();
    }
  void setStatus(const int i)
    {
      if (i) 
	{
	  plstatus->SetFillColor(3);
	  plstatus->SetLabel("Running");
	}
      else  
	{
	  plstatus->SetFillColor(41);
	  plstatus->SetLabel("Stopped");
	}
    }
  
  void setEvtnr(const int i)
    {
      char s[32];
      sprintf(s,"%d",i);
      plevt->SetLabel(s);
    }
  
  void setStream(const char * stream)
    {
      plstream->SetLabel(stream);
    }
  
};


#endif

  
