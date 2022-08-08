
########################################################################
# Copyright 1998-2021 CERN for the benefit of the EvtGen authors       #
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

# - CMake buildscript for Doxygen documentation subdirectory

# - Doxygen documentation
# Find an install of Doxygen (we're not too fussy about the version) using CMake's find_package command
find_package(Doxygen REQUIRED)

# Configure+Copy Doxygen template file from current source dir to current binary dir
# '@ONLY' means replace any instances of '@VARIABLE@' in the input file with
# the value of the CMake variable named 'VARIABLE' in the output file
configure_file(Doxyfile.in Doxyfile @ONLY)
configure_file(mainpage.dox.in mainpage.dox @ONLY)
configure_file(../AUTHORS AUTHORS.md @ONLY)

# Run Doxygen
add_custom_command(
  # What this command generates
  OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/html/index.html"
  # Actually run doxygen
  COMMAND ${DOXYGEN_EXECUTABLE}
  # ... in the same directory as we generated the Doxyfile...
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  # Output requires regeneration when
  # i   - we modify the Doxygen configuration
  # ii  - the sources being documented change
  # iii - any auxillary inputs change
  DEPENDS Doxyfile.in mainpage.dox.in ${PROJECT_SOURCE_DIR}/AUTHORS ${PROJECT_SOURCE_DIR}/README.md ${PROJECT_SOURCE_DIR}/Pythia8_README.md ${PROJECT_SOURCE_DIR}/Tauola_README.md ${PROJECT_SOURCE_DIR}/History.md ${PROJECT_SOURCE_DIR}/EvtGen ${PROJECT_SOURCE_DIR}/EvtGenBase ${PROJECT_SOURCE_DIR}/EvtGenModels ${PROJECT_SOURCE_DIR}/EvtGenExternal ${PROJECT_SOURCE_DIR}/src
  # Log that the command is run
  COMMENT "Building Doxygen for ${PROJECT_NAME}"
  )

# However, adding a command doesn't change anything in the build
# system - need to create a target that depends on the output of the
# command.
add_custom_target(
  # Name of the target - e.g. in Makefiles we can run 'make doc'
  doc
  # Add this target to the default target - i.e. documentation is
  # always built when building the whole project
  ALL
  # Make it depend on the output of the custom command
  DEPENDS "${CMAKE_CURRENT_BINARY_DIR}/html/index.html"
  )

# Install the documentation
install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/html DESTINATION ${CMAKE_INSTALL_DOCDIR})

