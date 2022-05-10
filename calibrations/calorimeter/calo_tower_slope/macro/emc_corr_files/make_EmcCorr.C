
{

   TFile f("inv_emcal_corr_sec12bands.root","recreate");
   TTree t1("emc_corr_tree","");
   Int_t towid;
   Float_t corr;
   t1.Branch("corr",&corr,"corr/F");
   t1.Branch("towid",&towid,"towid/I");

   for (Int_t i=0;i< 96; i++) {

     int k = i/8;
     //     int ll = (6-abs(6-k));
     int ll = k%6;

     for (int j=0; j < 256; j++) {
       
       towid = i*1000 + j;
       corr = 1.0/(0.9+ll*0.2);
       //       corr = 0.9+ll*0.2;
       //corr = 0.3000000;
       if (i%8==0&& j%128==0) 
	 cout << "towid" << towid <<  " corr  " << corr  << endl;
       t1.Fill();

     }
   }
  
   f.Write();
   //   f.Close();
 

}
