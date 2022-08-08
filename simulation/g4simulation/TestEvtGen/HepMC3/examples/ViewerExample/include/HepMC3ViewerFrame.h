// -*- C++ -*-
//
// This file is part of HepMC
// Copyright (C) 2014-2019 The HepMC collaboration (see AUTHORS for details)
//
///
/// @class HepMC3ViewerFrame
/// @brief Definition of \b class HepMC3ViewerFrame used for simple GUI viewer
///
#include <TGClient.h>
#include <TBuffer.h>
#include <TGButton.h>
#include <TGFrame.h>
#include "TImage.h"
#include "TCanvas.h"
#include "TGCanvas.h"
#include "TRootEmbeddedCanvas.h"
#include <TGClient.h>
#include <TCanvas.h>
#include <TBuffer.h>
#include <TGButton.h>
#include <TGFrame.h>
#include "TROOT.h"
#include "TImage.h"
#include "TH1S.h"
#include "TGFileDialog.h"

#include "HepMC3/GenEvent.h"
#include "HepMC3/Reader.h"
///
/// @class HepMC3ViewerFrame
/// @brief Definition of \b class HepMC3ViewerFrame
///
class HepMC3ViewerFrame : public TGMainFrame
{
private:
    TGCompositeFrame *fMainFrame; ///< Main frame
    TGCompositeFrame *fButtonFrame;  ///< Button frame
    TGTextButton     *fNextEvent; ///< Button
    TGTextButton *fPreviousEvent; ///< Button
    TGTextButton *fExit; ///< Button
    TGTextButton *fChooseInput; ///< Button
    TGTextButton *fClearEventCache; ///< Button
    TRootEmbeddedCanvas *fEmbEventImageCanvas;   ///< Event canvas
    TRootEmbeddedCanvas *fEmbAnalysisCanvas;     ///< Analysis canvas
    std::shared_ptr<HepMC3::Reader> fReader;                  ///< Reader
    HepMC3::GenEvent *fCurrentEvent;                          ///<Event
    std::vector<HepMC3::GenEvent*> fEventsCache;              ///<Cache of events
    TCanvas* fEventImageCanvas;                               ///< Event canvas
    TCanvas *fAnalysisCanvas;                                 ///<Analysis canvas
    TImage *fGraphImage;                                      ///<Image passed from graphviz
    std::map<std::string, TH1*> fAnalysisH;                   ///< Analysis histograms
    static const size_t m_char_buffer_size=100000;            ///<Size of writer buffer
public:
    void ReadFile(const char* a);                             ///< Open file
    HepMC3ViewerFrame(const TGWindow *p, UInt_t w, UInt_t h); ///< Constructor
    virtual ~HepMC3ViewerFrame();                             ///< Destructor
//Helper functions
//To get image from graphviz
    void DrawEvent(); ///< Draw evemt
//To do extra analysiz of the event
    void DoAnalysis(); ///< Do analysis
    // slots
    void NextEvent();          ///< slot
    void PreviousEvent();///< slot
    void ClearEventCache();///< slot
    void ChooseInput();///< slot
//   ClassDef(HepMC3ViewerFrame, 0)
};
