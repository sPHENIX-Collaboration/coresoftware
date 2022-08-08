#!/usr/bin/env python

from __future__ import division, print_function
import sys
import re


class PPupdate:
  """Code to update particle properties to input from PDG."""
  pdgProperties = dict()
  hbar = 6.58212E-25
#  hbar = 6.58211928E-25
  def readPDGfile(self, filename):
    ff = open(filename)
    for line in ff:
      if line[0] == '*':  # Comment
        continue
      if len(line) < 107:
        continue
      ids = re.split(' +',line[0:32].strip())
      mass = float( line[32:51].strip() )
      width = float( '0'+line[70:88].strip() )
      name = line[107:len(line)]
      for i in ids:
        self.pdgProperties[int(i)] = [ mass, width ]
    ff.close()

  def updateProperties(self, inputfile, outputfile):
    ffin = open( inputfile, 'r' )
    ffout = open( outputfile, 'w' )
    for line in ffin:
      fields = re.split( ' +', line.strip() )
      if not (fields[0]=='add' and fields[2]=='Particle'):
        ffout.write(line)
        continue
      idKey = int(fields[4])
      if idKey not in self.pdgProperties.keys():
        idKey *= -1
      if idKey not in self.pdgProperties.keys():
        ffout.write(line) 
        continue
      mass = self.pdgProperties[idKey][0]
      width = self.pdgProperties[idKey][1]
      oldLifetime = float( fields[10] )
      oldWidth = float( fields[6] )
      lifetime = 0
      if idKey == 4212:
        print(1, width, oldWidth, lifetime, oldLifetime)
      if width != 0 and width < 1e-5 and oldLifetime > 2e-15:
        lifetime = self.hbar/width
        if oldWidth < 2e-20:
          width = 0
      if idKey == 4212:
        print(2, width, oldWidth, lifetime, oldLifetime)
      if width == 0 and oldLifetime > 2e-15:
        lifetime = oldLifetime
      if idKey == 4212:
        print(3, width, oldWidth, lifetime, oldLifetime)
      if oldLifetime > 2e-18 and (( lifetime - oldLifetime )/oldLifetime < 1e-5):
        lifetime = oldLifetime
      if idKey == 4212:
        print(4, width, oldWidth, lifetime, oldLifetime)
      if width == 0 or width < 2e-24 or oldLifetime > 2e-15:
        width = oldWidth
      if idKey == 4212:
        print(5, width, oldWidth)
#      print(idKey, mass, fields[3], width, lifetime, fields[10])
      massStr = '%11.7e' % mass
      widthStr = '%11.7e' % width
      lifetimeStr = '%11.7e' % lifetime
      ffout.write( '%s%s%s%s%s%s%s%s%s%s%s%s\n' % 
                                   ( '{:<5}'.format(fields[0]), 
                                     '{:<2}'.format(fields[1]) ,
                                     '{:<10}'.format(fields[2]),
                                     '{:<22}'.format(fields[3]),
                                     '{:>14}'.format(fields[4]),
                                     '{:>15}'.format(massStr),
                                     '{:>15}'.format(widthStr),
                                     '{:>15}'.format(fields[7]),
                                     '{:>6}'.format(fields[8]),
                                     '{:>6}'.format(fields[9]),
                                     '{:>15}'.format(lifetimeStr),
                                     '{:>11}'.format(fields[11])
                                   ))
    ffin.close()
    ffout.close() 

updater = PPupdate()
updater.readPDGfile(sys.argv[1])
updater.updateProperties(sys.argv[2], sys.argv[3])

