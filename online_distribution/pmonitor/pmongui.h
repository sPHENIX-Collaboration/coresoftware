#ifndef __PMONGUI_H__
#define __PMONGUI_H__

#include <TGClient.h>
#include <TGButton.h>
#include <TGLabel.h>


class pmongui : public TGMainFrame {  //++CINT

private:
  TGTextButton *fButton1;
  TGTextButton *fButton2;
  TGTextButton *fButton3;
  TGLabel *statuslabel;
  TGLabel *streamlabel;
  TGLabel *evtnrlabel;


  TGLayoutHints   *fLayout;

public:
  pmongui(const TGWindow *p, UInt_t w, UInt_t h);
  ~pmongui();
  Bool_t ProcessMessage(Long_t msg, Long_t parm1, Long_t parm2);
  int setStatusLabel(const char * status);
  int setStreamLabel(const char * streamname);
  int setEvtnrLabel(const int n);
};

#endif /* __PMONGUI_H__ */
