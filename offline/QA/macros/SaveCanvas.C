// $Id: $

/*!
 * \file SaveCanvas.C
 * \brief Save canvas as png, eps, cxx, root, etc.
 * \author Jin Huang <jhuang@bnl.gov>
 * \version $Revision:   $
 * \date $Date: $
 */

#ifndef SaveCanvas_C

#define SaveCanvas_C

#include <TList.h>
#include <TClass.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TString.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TStyle.h>
#include <TObject.h>
#include <TH1F.h>

#include <iostream>

using namespace std;

//! Service function to SaveCanvas()
void
SavePad(TPad * p)
{
  if (!p)
    return;

  TList * l = p->GetListOfPrimitives();
//  l->Print();

  TIter next(l);
  TObject *obj = NULL;
  while ((obj = next()))
    {

      if (obj->IsA()->GetBaseClassOffset(TClass::GetClass("TPad")) >= 0)
        {
          if ((TPad *) obj != p)
            SavePad((TPad *) obj);
        }
      else if (obj->IsA()->GetBaseClassOffset(TClass::GetClass("TH1")) >= 0)
        {
          cout << "Save TH1 " << obj->GetName() << endl;
          obj->Clone()->Write(obj->GetName(), TObject::kOverwrite);
        }
      else if (obj->IsA()->GetBaseClassOffset(TClass::GetClass("TF1")) >= 0)
        {
          cout << "Save TF1 " << obj->GetName() << endl;
          obj->Clone()->Write(obj->GetName(), TObject::kOverwrite);
        }
      else if (obj->IsA()->GetBaseClassOffset(TClass::GetClass("TGraph")) >= 0)
        {
          cout << "Save TGraph " << obj->GetName() << endl;
          obj->Clone()->Write(obj->GetName(), TObject::kOverwrite);
        }
    }
}

//! Save canvas to multiple formats
/*!
 *  @param[in] c    pointer to the canvas
 *  @param[in] name Base of the file name. The default is the name of the cavas
 *  @param[in] bEPS true = save .eps and .pdf format too.
 */
void
SaveCanvas(TCanvas * c, TString name = "", Bool_t bEPS = kTRUE)
{
  if (name.Length() == 0)
    name = c->GetName();

  c->Print(name + ".png");

  TDirectory * oldd = gDirectory;

  TString rootfilename;

  c->Print(rootfilename = name + ".root");

  TFile f(rootfilename, "update");

  SavePad(c);

  f.Close();

  oldd->cd();

  if (bEPS)
    {
//      c->Print(name + ".pdf");

      float x = 20;
      float y = 20;
      gStyle->GetPaperSize(x, y);

      gStyle->SetPaperSize(c->GetWindowWidth() / 72 * 2.54,
          c->GetWindowHeight() / 72 * 2.54);
//      c->Print(name + ".eps");
      c->Print(name + ".svg");
      gSystem->Exec("rsvg-convert -f pdf -o "+name + ".pdf " + name + ".svg");
      gSystem->Exec("rm -fv " +  name + ".svg");

      gStyle->SetPaperSize(x, y);
    }
  //      c->Print(name+".C");
}

//! example to use this SaveCanvas()
/*!
 *  Output:
 *  The canvas data will be saved to RootFileName.root, as well as
 *  RootFileName.png for presentation and RootFileName.eps for Latex
 *
 *  How to use RootFileName.root
 *  open RootFileName.root with root.
 *  It contains the canvas object "CanvasTest", which can be redraw with CanvasTest -> Draw()
 *  It also contains a copy of the histograms and graphs data for use in root again, e.g.
 *  root [3]  h1->GetBinContent(30)
 *
 */
void
example_save_canvas()
{

  TCanvas *c1 = new TCanvas("CanvasTest", "CanvasTest", 800, 900);

  TH1F * h1 = new TH1F("h1", "histo from a gaussian", 100, -3, 3);
  h1->FillRandom("gaus", 10000);

  h1->Draw();

  // single call to save c1 to file RootFileName.*
  SaveCanvas(c1, "RootFileName");

}

#endif

