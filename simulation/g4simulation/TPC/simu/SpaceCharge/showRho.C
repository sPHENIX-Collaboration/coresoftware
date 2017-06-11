void showRho() {
  gStyle->SetOptStat(0);
  TFile *file = new TFile( "TPCCAGE_20_78_211_0.root" );
  TH3F *a = (TH3F*) file->Get("rho");
  a->GetYaxis()->SetRange(1,1);
  TH2D *a1 = (TH2D*) a->Project3D("xz");
  a->GetYaxis()->SetRange(1,-1);
  int bb = 0.5*a->GetZaxis()->GetNbins();
  a->GetZaxis()->SetRange(bb,bb);
  TH2D *a2 = (TH2D*) a->Project3D("xy");
  float rmax = 1.2*a1->GetYaxis()->GetBinCenter( a1->GetYaxis()->GetNbins() );
  TH2D* dummy_his = new TH2D("dummy", "histo title;x [cm];y [cm]", 100, -rmax, rmax, 100, -rmax, rmax);
  a1->SetTitle( a->GetTitle() );
  a2->SetTitle( a->GetTitle() );
  dummy_his->SetTitle( a->GetTitle() );
  a2->RebinY(200);
  a2->Scale(1/200.);

  TCanvas *main = new TCanvas();
  main->Divide(2,1);
  main->cd(1);
  a1->Draw("colz");
  main->cd(2);
  dummy_his->Draw("colz");
  a2->Draw("colz pol same");
}
