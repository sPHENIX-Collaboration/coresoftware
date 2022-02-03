#!/usr/bin/env python3
from ROOT import TCanvas, TFile, TH1F,TH2F,TH3F, TF1, TGraph, gROOT
import os
import re
import glob

gROOT.SetBatch(True)

dirName = '/sphenix/user/shulga/Work/IBF/DistortionMap/Files/'
bXs = [1508071, 3016509, 4524020, 6032112, 7540028, 9048092, 10556072, 12064371, 13572143, 15080178, 16588072, 18096105]
#bXs = [18096105]
h_names = ['_h_hits','_h_R','_h_DC_E']
for i in range(30):
    h_names.append('_h_SC_ibf_{}'.format(i))
    h_names.append('_h_SC_prim_{}'.format(i))


for ib,bX in enumerate(bXs):
    print(bX)
    name = 'mdc2_ADCBins_UseFieldMaps_hist_G4Hits_sHijing_0-12fm_bX{}*'.format(bX)
    outputName = '/sphenix/user/shulga/Work/IBF/DistortionMap/Files/Summary_hist_mdc2_UseFieldMaps_AA_event_{}_bX{}.root'.format(ib,bX)

    filePattern = dirName+name
    files = sorted(glob.glob(filePattern))
    #print(files)
    #n=0
    histos = []
    for n,file in enumerate(files):
        #h=0
        f = TFile.Open(file)
        for h,h_name in enumerate(h_names):
            newName=h_name+'_{}'.format(n)
            if n==0:
                newName=h_name
            hist = f.Get(h_name).Clone(h_name)
            if n==0:
                histos.append(hist)
            if n>0:
                histos[h].Add(hist) 
            hist.SetDirectory(0)
            #h+=1
        #n+=1

    outfile = TFile(outputName, "RECREATE")
    for hist in histos:
        hist.Sumw2(False)
        hist.Write()
    outfile.Write()
    outfile.Close()
#print(files)
#
#