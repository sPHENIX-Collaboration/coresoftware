#include <TApplication.h>
#include "HepMC3ViewerFrame.h"
int main(int argc, char **argv)
{
    if (argc>2)
    {
        fprintf(stderr, "%s: only one optional argument is supported: the name of file to open.\n", argv[0]);
        return 1;
    }
    TApplication theApp("App", &argc, argv);

    if (gROOT->IsBatch())
    {
        fprintf(stderr, "%s: cannot run in batch mode\n", argv[0]);
        return 1;
    }
    HepMC3ViewerFrame *G=new HepMC3ViewerFrame(gClient->GetRoot(), 350, 80);
    if (theApp.Argc()>1) G->ReadFile(theApp.Argv()[1]);
    theApp.Run();
}