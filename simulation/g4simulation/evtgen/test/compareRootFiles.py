#!/usr/bin/env python

########################################################################
# Copyright 1998-2020 CERN for the benefit of the EvtGen authors       #
#                                                                      #
# This file is part of EvtGen.                                         #
#                                                                      #
# EvtGen is free software: you can redistribute it and/or modify       #
# it under the terms of the GNU General Public License as published by #
# the Free Software Foundation, either version 3 of the License, or    #
# (at your option) any later version.                                  #
#                                                                      #
# EvtGen is distributed in the hope that it will be useful,            #
# but WITHOUT ANY WARRANTY; without even the implied warranty of       #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the        #
# GNU General Public License for more details.                         #
#                                                                      #
# You should have received a copy of the GNU General Public License    #
# along with EvtGen.  If not, see <https://www.gnu.org/licenses/>.     #
########################################################################

import sys
import ROOT

def compare1D(can, hist1, hist2):
  ''' Function to overlay two 1D histograms'''
  hist1.SetLineColor(ROOT.kBlue)
  hist2.SetLineColor(ROOT.kRed)
  hist1.Draw()
  hist2.Draw('same')
  can.Update()

def compare2D(can, hist1, hist2, pdf):
  ''' This function plots difference of 2D histograms. Be aware that in
case of identical histograms result will be empty plot.'''
  histDiff = hist1.Clone(hist1.GetName()+'diff')
  histDiff.Add(hist1, hist2, 1, -1)
  histDiff.SetMinimum(-50)
  histDiff.SetMaximum(50)
  histDiff.Draw('colz')
  can.Update()
  canvas.Print(pdf)

if __name__ == '__main__':
  file1 = ROOT.TFile.Open(sys.argv[1])
  file2 = ROOT.TFile.Open(sys.argv[2])
  pdfName = sys.argv[3]
  
  histNames = []
  histTypes = []
  
  for key in file1.GetListOfKeys():
    histNames.append(key.GetName())
    histTypes.append(key.GetClassName())
  
  print(histNames)
  print(histTypes)
  
  canvas = ROOT.TCanvas()
  canvas.Draw()
  canvas.Print(pdfName+'[')
  
  for name, typ in zip(histNames, histTypes):
    hh1 = file1.Get(name)
    hh2 = file2.Get(name)
    if 'TH1' in typ:
      compare1D(canvas, hh1, hh2)
      canvas.Print(pdfName)
    elif 'TH2' in typ:
      compare2D(canvas, hh1, hh2, pdfName)
    print(hh1.GetName())
    hh1.Chi2Test(hh2, "P")
  
  canvas.Print(pdfName+']')
  file1.Close()
  file2.Close()
  
