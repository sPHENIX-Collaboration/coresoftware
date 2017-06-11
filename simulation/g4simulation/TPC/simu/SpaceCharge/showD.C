void showD() {
  gStyle->SetOptStat(0);
  //TFile *file = new TFile( "sPHENIX20.root" );
  TFile *file = new TFile( "TPCCAGE_20_78_211_2.root" );

  TH3F *a = (TH3F*) file->Get("mapDeltaR");
  a->GetYaxis()->SetRange(1,1);
  TH2D *a1 = (TH2D*) a->Project3D("xz");
  a->GetYaxis()->SetRange(1,-1);
  int aa = 0.5*a->GetZaxis()->GetNbins();
  a->GetZaxis()->SetRange(aa,aa);
  TH2D *a2 = (TH2D*) a->Project3D("xy");

  TH3F *b = (TH3F*) file->Get("mapRDeltaPHI");
  b->GetYaxis()->SetRange(1,1);
  TH2D *b1 = (TH2D*) b->Project3D("xz");
  b->GetYaxis()->SetRange(1,-1);
  int bb = 0.5*b->GetZaxis()->GetNbins();
  b->GetZaxis()->SetRange(bb,bb);
  TH2D *b2 = (TH2D*) b->Project3D("xy");

  float rmax = 1.2*a1->GetYaxis()->GetBinCenter( a1->GetYaxis()->GetNbins() );
  TH2D* dummy_his = new TH2D("dummy", "histo title;x [cm];y [cm]", 100, -rmax, rmax, 100, -rmax, rmax);

  a1->SetTitle( a->GetTitle() );
  a2->SetTitle( a->GetTitle() );
  b1->SetTitle( b->GetTitle() );
  b2->SetTitle( b->GetTitle() );

  TCanvas *main = new TCanvas();
  main->Divide(2,1);
  main->cd(1);
  a1->Draw("colz");
  main->cd(2);
  dummy_his->SetTitle( a->GetTitle() );
  dummy_his->Draw("colz");
  a2->Draw("colz pol same");

  TCanvas *main2 = new TCanvas();
  main2->Divide(2,1);
  main2->cd(1);
  b1->Draw("colz");
  main2->cd(2);
  dummy_his->SetTitle( b->GetTitle() );
  dummy_his->Draw("colz");
  b2->Draw("colz pol same");
}
