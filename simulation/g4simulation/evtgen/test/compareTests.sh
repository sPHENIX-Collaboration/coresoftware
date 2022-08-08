#!/usr/bin/bash

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

oldDir=$1
newDir=$2
pdfDir=$3

# First we do all root files
for ff in ${oldDir}/*root
do
  fileName=$(basename $ff)
  baseFileName=`echo ${fileName}|cut -d'.' -f1`
  ./compareRootFiles.py ${oldDir}/${fileName} ${newDir}/${fileName} ${pdfDir}/${baseFileName}.pdf
done

# Now compare all dat files
for ff in ${oldDir}/*dat
do
  fileName=$(basename ${ff})
  echo $fileName `diff ${oldDir}/${fileName} ${newDir}/${fileName}|wc -l `
done

