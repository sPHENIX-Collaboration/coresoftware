#include "HepMC3/ReaderRootTree.h"
#include "HepMC3/ReaderFactory.h"
#include "HepMC3/GenEvent.h"
#include "HepMC3/GenParticle.h"
#include "HepMC3/GenVertex.h"
#include "HepMC3ViewerFrame.h"

/* Older graphviz versions can have conflictiong declarations  of memcmp/strcmp function
 * This can break compilation with -pedantic. Uncomenting line below can fix it.
 */
// #define _PACKAGE_ast 1

#include <graphviz/gvc.h>
#define CONSERVATION_TOLERANCE 1e-5

static  char*  create_image_from_dot(char* m_buffer)
{
    GVC_t * gvc=gvContext();
    Agraph_t *g= agmemread(m_buffer);
    gvLayout(gvc,g,"dot");

    int err;
    char *data;
    unsigned int length;

    if (!g)
        return NULL;
    err = gvRenderData(gvc, g, "png", &data, &length);
    if (err)
        return NULL;
    data = (char*)realloc(data, length + 1);
    delete g;
    delete gvc;
    return data;
}

static bool show_as_parton(HepMC3::ConstGenParticlePtr p )
{
    const int pd=std::abs(p->pid());
    bool parton=false;

    if (pd==81||pd==82||pd<25) parton=true;
    if (
        (pd/1000==1||pd/1000==2||pd/1000==3||pd/1000==4||pd/1000==5)
        &&(pd%1000/100==1||pd%1000/100==2||pd%1000/100==3||pd%1000/100==4)
        &&(pd%100==1||pd%100==3)
    )
        parton=true;
    if (p->status()==4)  parton=true;
    return parton;
}

static char*  write_event_to_dot(char* used_cursor,const HepMC3::GenEvent &evt,int used_style=1)
{
    used_cursor += sprintf(used_cursor, "digraph graphname%d {\n",evt.event_number());
    used_cursor += sprintf(used_cursor, "v0[label=\"Machine\"];\n");
    for(auto v: evt.vertices() )
    {
        if (used_style!=0)
        {
            if (used_style==1) //paint decay and fragmentation vertices in green
            {
                if (v->status()==2) used_cursor += sprintf(used_cursor, "node [color=\"green\"];\n");
                else  used_cursor += sprintf(used_cursor, "node [color=\"black\"];\n");
            }
        }
        HepMC3::FourVector in=HepMC3::FourVector(0,0,0,0);
        HepMC3::FourVector out=HepMC3::FourVector(0,0,0,0);
        double energy=0;
        for(auto p1: v->particles_in()  ) {
            in+=p1->momentum();
            energy+=std::abs(p1->momentum().e());
        }
        for(auto p2: v->particles_out() ) {
            out+=p2->momentum();
            energy+=std::abs(p2->momentum().e());
        }
        HepMC3::FourVector momviolation(0,0,0,0);
        momviolation+=in;
        momviolation-=out;
        double energyviolation=std::sqrt(momviolation.length2()  +momviolation.e()*momviolation.e()       );
        bool violation=false;
        if (energyviolation>CONSERVATION_TOLERANCE*energy) violation=true;

        if(violation)
        {
            used_cursor += sprintf(used_cursor, "node [shape=rectangle];\n");
            used_cursor += sprintf(used_cursor, "v%d [label=\"%d\nd=%4.2f\"];\n", -v->id(),v->id(),energyviolation);
        }
        else
        {
            used_cursor += sprintf(used_cursor, "node [shape=ellipse];\n");
            used_cursor += sprintf(used_cursor, "v%d[label=\"%d\"];\n", -v->id(),v->id());
        }

        used_cursor += sprintf(used_cursor, "node [shape=ellipse];\n");
    }
    for(auto p: evt.beams() )
    {
        if (!p->end_vertex()) continue;
        used_cursor += sprintf(used_cursor, "node [shape=point];\n");
        used_cursor += sprintf(used_cursor, "v0 -> v%d [label=\"%d(%d)\"];\n", -p->end_vertex()->id(),p->id(),p->pid());
    }

    for(auto v: evt.vertices() )
    {

        for(auto p: v->particles_out() )
        {
            {
                if (used_style!=0)
                {
                    if (used_style==1) //paint suspected partons and 81/82 in red
                    {
                        if (show_as_parton(p)&&p->status()!=1) used_cursor += sprintf(used_cursor, "edge [color=\"red\"];\n");
                        else        used_cursor +=sprintf(used_cursor, "edge [color=\"black\"];\n");
                    }
                }
                if (!p->end_vertex())
                {
                    used_cursor += sprintf(used_cursor, "node [shape=point];\n");
                    used_cursor += sprintf(used_cursor, "v%d -> o%d [label=\"%d(%d)\"];\n", -v->id(),p->id(),p->id(),p->pid());
                    continue;
                }
                else
                    used_cursor += sprintf(used_cursor, "v%d -> v%d [label=\"%d(%d)\"];\n", -v->id(),-p->end_vertex()->id(),p->id(),p->pid());
            }
        }
    }
    used_cursor += sprintf(used_cursor, "labelloc=\"t\";\nlabel=\"Event %d; Vertices %lu; Particles %lu;\";\n", evt.event_number(), evt.vertices().size(), evt.particles().size());
    used_cursor += sprintf(used_cursor,"}\n\n");

    return used_cursor;
}


void HepMC3ViewerFrame::DrawEvent()
{
    char* m_buffer = new char[m_char_buffer_size]();
    char* m_cursor=m_buffer;
    m_cursor=write_event_to_dot(m_cursor,*(fCurrentEvent));
    char *buf=create_image_from_dot(m_buffer);
    fEmbEventImageCanvas->MapSubwindows();

    if(!fEventImageCanvas)  fEventImageCanvas=new TCanvas("fEmbEventImageCanvas","fEmbEventImageCanvas",1024,768);

    fEventImageCanvas->cd();
    fEventImageCanvas->Clear();
    double d=0.60;

    fGraphImage = TImage::Create();
    fGraphImage->SetName("Event");
    fGraphImage->SetImageBuffer(&buf, TImage::kPng);

    fGraphImage->SetConstRatio(kFALSE);

    TPad *p1 = new TPad("i1", "i1", 0.05, 0.05, 0.05+d*fGraphImage->GetWidth()/fGraphImage->GetHeight(), 0.95);
    p1->Draw();
    p1->cd();

    fGraphImage->Draw("xxx");
    delete [] m_buffer;
    gPad->Update();
    DoAnalysis();
}

void HepMC3ViewerFrame::DoAnalysis()
{
    fEmbAnalysisCanvas->MapSubwindows();
    fAnalysisCanvas->cd();
    fAnalysisCanvas->Clear();
    for (auto h: fAnalysisH) h.second->Delete();
    fAnalysisH.clear();

    /*   */
    TH1S* particles1= new TH1S();
    fAnalysisH["particles1"]=particles1;
    particles1->SetTitle("Flavour: all particles; PDG ID; Number of particles");
    particles1->SetFillColor(kBlue);
    for(auto p: fCurrentEvent->particles() )
        particles1->Fill((std::to_string(p->pid())).c_str(),1.0);
    particles1->LabelsOption(">","X");
    /*   */
    TH1S* particles2= new TH1S();
    fAnalysisH["particles2"]=particles2;
    particles2->SetTitle("Flavour: particles with status 1; PDG ID; Number of particles");
    particles2->SetFillColor(kBlue);
    for(auto p: fCurrentEvent->particles() )
        if(p->status()==1) particles2->Fill((std::to_string(p->pid())).c_str(),1.0);
    particles2->LabelsOption(">","X");
    /*   */
    std::vector<double> masses;
    for(auto p: fCurrentEvent->particles() )
        if(show_as_parton(p)) masses.push_back(p->momentum().m());
    TH1D* particles3= new TH1D("particles3","Mass:  parton particles; Mass, GeV; Number of particles",masses.size(),
                               0,*std::max_element(masses.begin(),masses.end()));
    fAnalysisH["particles3"]=particles3;
    particles3->SetFillColor(kBlue);
    for(auto m: masses) particles3->Fill(m);


    fAnalysisCanvas->cd();
    TPad *p1 = new TPad("i1", "i1", 0.00, 0.75, 1.0, 1.0);
    p1->Draw();
    p1->cd();
    particles1->Draw();
    fAnalysisCanvas->cd();
    TPad *p2 = new TPad("i2", "i2", 0.00, 0.50, 1.0, 0.75);
    p2->Draw();
    p2->cd();
    particles2->Draw();
    fAnalysisCanvas->cd();
    TPad *p3 = new TPad("i3", "i3", 0.00, 0.25, 1.0, 0.50);
    p3->Draw();
    p3->cd();
    particles3->Draw();

    gPad->Update();
}

void HepMC3ViewerFrame::ClearEventCache()
{
    fEventsCache.clear();
    fCurrentEvent=nullptr;
}

void HepMC3ViewerFrame::PreviousEvent()
{
    auto pos=find(fEventsCache.begin(),fEventsCache.end(),fCurrentEvent);
    if (pos==fEventsCache.begin()) return;
    pos--;
    fCurrentEvent=*(pos);
    if (pos==fEventsCache.end()) printf("This event was not found in the cache.  Cache size is %zu\n",fEventsCache.size());
    DrawEvent();
}

void HepMC3ViewerFrame::ReadFile(const char* a) {
    fReader=HepMC3::deduce_reader(a);
}

void HepMC3ViewerFrame::NextEvent()
{
    if (fCurrentEvent==nullptr||fEventsCache.back()==fCurrentEvent)
    {
        HepMC3::GenEvent* evt1=new HepMC3::GenEvent(HepMC3::Units::GEV,HepMC3::Units::MM);
        bool ok=fReader->read_event(*(evt1));
        ok=(ok&&!fReader->failed());
        if (ok)
        {
            fEventsCache.push_back(evt1);
            fCurrentEvent=evt1;
        }
        else return;
    }
    else
    {
        auto pos=find(fEventsCache.begin(),fEventsCache.end(),fCurrentEvent);
        pos++;
        fCurrentEvent=*(pos);
    }
    DrawEvent();
}
void HepMC3ViewerFrame::ChooseInput()
{
    static const char *FileType[] = {"All", "*.*","HepMC", "*.hepmc*","LHEF", "*.lhe*","ROOT", "*.root", 0, 0 };
    static TString dir("./");
    TGFileInfo fi;
    fi.fFileTypes = FileType;
    fi.fIniDir = StrDup(dir);
    new TGFileDialog(gClient->GetRoot(), this, kFDOpen, &fi);
    if (fReader) fReader->close();
    fReader=HepMC3::deduce_reader(fi.fFilename);
}

HepMC3ViewerFrame::HepMC3ViewerFrame(const TGWindow *p, UInt_t w, UInt_t h) :
    TGMainFrame(p, w, h)
{
    fMainFrame = new TGCompositeFrame(this, 1350, 500, kHorizontalFrame|kFixedWidth);
    fButtonFrame = new TGCompositeFrame(fMainFrame, 150, 200, kFixedWidth);

    fEmbEventImageCanvas =new TRootEmbeddedCanvas("MainCanvaslegent", fMainFrame, 850, 500);

    fEmbAnalysisCanvas =new TRootEmbeddedCanvas("EmbAnalysisCanvaslegend", fMainFrame, 350, 500);


    fMainFrame->AddFrame(fEmbEventImageCanvas,new TGLayoutHints(kLHintsTop | kLHintsExpandX| kLHintsExpandY, 1, 1, 2, 2));
    fMainFrame->AddFrame(fEmbAnalysisCanvas,new TGLayoutHints(kLHintsTop | kFixedWidth| kLHintsExpandY, 1, 1, 2, 2));
    fMainFrame->AddFrame(fButtonFrame,new TGLayoutHints(kLHintsTop, 1, 1, 2, 2));


    fChooseInput = new TGTextButton(fButtonFrame, "&Choose input");
    fChooseInput->Connect("Clicked()", "HepMC3ViewerFrame", this, "ChooseInput()");
    fChooseInput->SetToolTipText("Click to choose file");
    fButtonFrame->AddFrame(fChooseInput, new TGLayoutHints(kLHintsTop | kLHintsExpandX, 1, 1, 2, 2));



    fNextEvent = new TGTextButton(fButtonFrame, "&Next event");
    fNextEvent->Connect("Clicked()", "HepMC3ViewerFrame", this, "NextEvent()");
    fNextEvent->SetToolTipText("Click to display next event");
    fButtonFrame->AddFrame(fNextEvent, new TGLayoutHints(kLHintsExpandX|kLHintsLeft, 1, 1, 2, 2));


    fPreviousEvent = new TGTextButton(fButtonFrame, "&Previous event");
    fPreviousEvent->Connect("Clicked()", "HepMC3ViewerFrame", this, "PreviousEvent()");
    fPreviousEvent->SetToolTipText("Click to display previous event");
    fButtonFrame->AddFrame(fPreviousEvent, new TGLayoutHints( kLHintsExpandX|kLHintsLeft, 1, 1, 2, 2));


    fClearEventCache = new TGTextButton(fButtonFrame, "&Clear event cache");
    fClearEventCache->Connect("Clicked()", "HepMC3ViewerFrame", this, "ClearEventCache()");
    fClearEventCache->SetToolTipText("Click to clear event cache ");
    fButtonFrame->AddFrame(fClearEventCache, new TGLayoutHints( kLHintsExpandX|kLHintsLeft, 1, 1, 2, 2));

    fExit = new TGTextButton(fButtonFrame, "&Exit ","gApplication->Terminate(0)");
    fExit->SetToolTipText("Click to exit");
    fButtonFrame->AddFrame(fExit, new TGLayoutHints( kLHintsExpandX|kLHintsLeft,1,1,2,2));

    AddFrame(fMainFrame, new TGLayoutHints(kLHintsTop |kLHintsExpandX| kLHintsExpandY, 1, 1, 2, 2));

    SetWindowName("Event viewer");
    MapSubwindows();
    Resize(GetDefaultSize());
    MapWindow();

    fReader=nullptr;
    fEventImageCanvas=fEmbEventImageCanvas->GetCanvas();
    fAnalysisCanvas=fEmbAnalysisCanvas->GetCanvas();
    fCurrentEvent=nullptr;
    fGraphImage = TImage::Create();
}
HepMC3ViewerFrame::~HepMC3ViewerFrame()
{
    fMainFrame->Cleanup();
    fReader->close();
    Cleanup();
}
