/**
 * Using 16-bit short integers to store 32-bit floating-point values.
 * Joint work with Massachusetts Institute of Technology (MIT) and Institute for Information Industry (III)
 * Questions please contact fishyu@iii.org.tw
 * Date: May 22, 2021
 */

#include "compressor.h"

#include <TFile.h>
#include <TTree.h>

#include <iostream>
#include <fstream>
#include <string>
#include <utility>

using namespace std;

void doCompression(ofstream &myfile, TString filename);
Float_t computeSd(Float_t* avg, Int_t n_entries, TTree* t, Float_t* gen_);

void output_before_vs_after(
 TString filename, 
 vector<UShort_t>* order, 
 vector<Float_t>* dict, 
 Int_t n_entries, 
 TTree* t, 
 Float_t* gen_
);

uint64_t timeSinceEpochMillisec()
{
   using namespace std::chrono;
   return duration_cast<milliseconds>(system_clock::now().time_since_epoch()).count();
}

void compress_clu_res_float32(std::string &filename)
{
 ofstream logFile;
 logFile.open ("log.csv");
 logFile << "file,points,variable,dltSD,bits,dltdltSD/dltSD,avg,LB,UB,msec" << endl;
 doCompression(logFile, filename);
 logFile.close();
}

void doCompression(ofstream &myfile, TString filename)
{
 TFile *f = new TFile(filename,"read");
 TTree *t = (TTree*)f->Get("ntp_trkres");
 Float_t genv1_;
 Float_t genv2_;
 //get branches
 t->SetBranchAddress("v1",&genv1_);
 t->SetBranchAddress("v2",&genv2_);
 Int_t n_entries = (Int_t)t->GetEntries();

 // start compress v1
 Float_t avg;
 Float_t dltSD = computeSd(&avg, n_entries, t, &genv1_);
 for (Int_t numBits = 3; numBits < 12; ++numBits) {
   vector<UShort_t> v1_order; // order, to replace v1 and save to the new ROOT file
   vector<Float_t> v1_dict; // dictionay, to save to the new ROOT file
   vector<size_t> v1_cnt; // how many data points in the dictionary
   uint64_t start = timeSinceEpochMillisec();
   Float_t dltdltSD = approx(&v1_order, &v1_dict, &v1_cnt, n_entries, t, &genv1_, (size_t) pow(2, numBits));
   myfile 
   << filename << "," 
   << v1_order.size() << ",dphi," 
   << dltSD << "," 
   << numBits << "," 
   << dltdltSD / dltSD << "," 
   << avg << "," 
   << avg - 4 * dltSD << "," 
   << avg + 4 * dltSD << "," 
   << timeSinceEpochMillisec() - start << endl;
   output_before_vs_after(filename + "_dltPHI_" + numBits + "bits.csv", &v1_order, &v1_dict, n_entries, t, &genv1_);
 }
 // end compress v1
 
 // start compress v2
 dltSD = computeSd(&avg, n_entries, t, &genv2_);
 for (Int_t numBits = 3; numBits < 12; ++numBits) {
   vector<UShort_t> v2_order; // order, to replace v2 and save to the new ROOT file
   vector<Float_t> v2_dict; // dictionay, to save to the new ROOT file
   vector<size_t> v2_cnt; // how many data points in the dictionary
   uint64_t start = timeSinceEpochMillisec();
   Float_t dltdltSD = approx(&v2_order, &v2_dict, &v2_cnt, n_entries, t, &genv2_, (size_t) pow(2, numBits));
   myfile 
     << filename << "," 
     << v2_order.size() << ",dz," 
     << dltSD << "," 
     << numBits << "," 
     << dltdltSD / dltSD << "," 
     << avg << "," 
     << avg - 4 * dltSD << "," 
     << avg + 4 * dltSD << "," 
     << timeSinceEpochMillisec() - start << endl;
     output_before_vs_after(filename + "_dltZ_" + numBits + "bits.csv", &v2_order, &v2_dict, n_entries, t, &genv2_);
  }
}

Float_t computeSd(Float_t* avg, Int_t n_entries, TTree* t, Float_t* gen_)
{
  Double_t sumSquared = 0;
  Double_t sum = 0;
  for (Int_t j = 0 ; j < n_entries; j++){
    t->GetEntry(j);
    
    sumSquared += (Double_t) *gen_ * (Double_t) *gen_;
    sum += (Double_t) *gen_;
  }
  
  Double_t _avg = sum / (Double_t) n_entries;
  *avg = _avg;

  return sqrt((sumSquared / (Double_t) n_entries) - _avg * _avg);
}

void output_before_vs_after(TString filename, vector<UShort_t>* order, vector<Float_t>* dict, Int_t n_entries, TTree* t, Float_t* gen_) 
{
  ofstream myfile;
  myfile.open(filename);
  myfile << "before,after" << endl;
  for (Int_t j = 0 ; j < n_entries; j++){
    t->GetEntry(j);
	myfile << *gen_ << "," << dict->at(order->at(j)) << endl;
  }
  myfile.close();
}

