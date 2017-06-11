void showDD() {
  gStyle->SetOptStat(0);
  TFile *file0 = new TFile( "sPHENIX20.root" );
  TFile *file1 = new TFile( "TPCCAGE_20_78_211_2.root" );
  TFile *file2 = new TFile( "TPCCAGE_30_78_211_2.root" );

  TH3F *a0 = (TH3F*) file0->Get("mapDeltaR");
  TH3F *a1 = (TH3F*) file1->Get("mapDeltaR");
  TH3F *a2 = (TH3F*) file2->Get("mapDeltaR");

  int nb0 = a0->GetZaxis()->FindBin( -0.25 );
  int nb1 = a1->GetZaxis()->FindBin( -0.25 );
  int nb2 = a2->GetZaxis()->FindBin( -0.25 );

  TH1D *a1_0 = (TH1D*) a0->ProjectionX("x0",1,1,nb0,nb0);
  TH1D *a1_1 = (TH1D*) a1->ProjectionX("x1",1,1,nb1,nb1);
  TH1D *a1_1c = (TH1D*) a1_1->Clone("r_ro");
  TH1D *a1_2 = (TH1D*) a2->ProjectionX("x2",1,1,nb2,nb2);

  a1_1->SetLineColor(kBlue-3);
  a1_2->SetLineColor(kRed-3);
  a1_2->SetLineStyle(2);
  a1_1c->SetLineColor(kBlue-3);
  a1_1c->SetFillColor(kOrange-3);
  a1_1c->SetFillStyle(3001);
  a1_1c->GetXaxis()->SetRangeUser(31.5,75.0);

  new TCanvas();
  TH2D* dummy_his = new TH2D("dummy", "Spacecharge Distortions in TPC;R [cm];#delta R [cm]", 100, 19.5, 78.5,100,-0.8,+2.5);
  dummy_his->Draw();
  //a1_0->Draw("same");
  a1_1c->Draw("][same");
  a1_1->Draw("same");
  a1_2->Draw("same");

  TLegend *leg = new TLegend(0.5,0.5,0.9,0.9, "Ne CF_{4} iC_{4}H_{10} @ 50 kHz IBF0.3\%");
  leg->SetFillColor(kWhite);
  leg->AddEntry(a1_2,"Inner Wall Rad at 30 cm");
  leg->AddEntry(a1_1,"Inner Wall Rad at 20 cm");
  leg->AddEntry(a1_1c,"Instrumented in sPHENIX");
  leg->Draw();

  TH3F *b0 = (TH3F*) file0->Get("mapRDeltaPHI");
  TH3F *b1 = (TH3F*) file1->Get("mapRDeltaPHI");
  TH3F *b2 = (TH3F*) file2->Get("mapRDeltaPHI");

  TH1D *b1_0 = (TH1D*) b0->ProjectionX("z0",1,1,nb0,nb0);
  TH1D *b1_1 = (TH1D*) b1->ProjectionX("z1",1,1,nb1,nb1);
  TH1D *b1_1c = (TH1D*) b1_1->Clone("rp_ro");
  TH1D *b1_2 = (TH1D*) b2->ProjectionX("z2",1,1,nb2,nb2);

  b1_1->SetLineColor(kBlue-3);
  b1_2->SetLineColor(kRed-3);
  b1_2->SetLineStyle(2);
  b1_1c->SetLineColor(kBlue-3);
  b1_1c->SetFillColor(kOrange-3);
  b1_1c->SetFillStyle(3001);
  b1_1c->GetXaxis()->SetRangeUser(31.5,75.0);

  new TCanvas();
  TH2D* dummy_his2 = new TH2D("dummy2", "Spacecharge Distortions in TPC;R [cm];R #delta #varphi [cm]", 100, 19.5, 78.5,100,-1.4,+5.2);
  dummy_his2->Draw();
  //b1_0->Draw("same");
  b1_1c->Draw("][same");
  b1_1->Draw("same");
  b1_2->Draw("same");

  leg->Draw();
}
