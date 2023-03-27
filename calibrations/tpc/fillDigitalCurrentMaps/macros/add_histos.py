#!/usr/bin/env python3
from ROOT import TCanvas, TFile, TH1F,TH2F,TH3F, TF1, TGraph, gROOT
import os
import re

gROOT.SetBatch(True)

dirName = './Files/'
bXs = [1508071, 3016509, 4524020, 6032112, 7540028, 9048092, 10556072, 12064371, 13572143, 15080178, 16588072, 18096105]
h_names = ['_h_DC_E',
            '_h_DC_SC',
            '_h_DC_SC_XY',
            '_h_R',
            '_h_SC_ibf',
            '_h_hits']

for bX in bXs:
    print(bX)
    name = 'hist_G4Hits_sHijing_0-12fm_bX{}'.format(bX)
    outputName = './Files/Summary_ADC_NoW_hist_AA_event_0_bX{}.root'.format(bX)
    files = [f for f in os.listdir(dirName) if re.match(name, f)]
    n=0
    histos = []
    for file in files:
        h=0
        f = TFile.Open(dirName+file)
        for h_name in h_names:
            newName=h_name+'_{}'.format(n)
            if n==0:
                newName=h_name
            hist = f.Get(h_name).Clone()
            hist.SetDirectory(0)
            if n==0:
                histos.append(hist)
            if n>0:
                histos[h].Add(hist) 
            h+=1
        n+=1

    outfile = TFile(outputName, "RECREATE")
    for hist in histos:
        hist.Sumw2(False)
        hist.Write()
    outfile.Write()
    outfile.Close()
#print(files)
#
#