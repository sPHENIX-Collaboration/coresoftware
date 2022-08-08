
/***********************************************************************
* Copyright 1998-2020 CERN for the benefit of the EvtGen authors       *
*                                                                      *
* This file is part of EvtGen.                                         *
*                                                                      *
* EvtGen is free software: you can redistribute it and/or modify       *
* it under the terms of the GNU General Public License as published by *
* the Free Software Foundation, either version 3 of the License, or    *
* (at your option) any later version.                                  *
*                                                                      *
* EvtGen is distributed in the hope that it will be useful,            *
* but WITHOUT ANY WARRANTY; without even the implied warranty of       *
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the        *
* GNU General Public License for more details.                         *
*                                                                      *
* You should have received a copy of the GNU General Public License    *
* along with EvtGen.  If not, see <https://www.gnu.org/licenses/>.     *
***********************************************************************/

#include "EvtGenExternal/EvtPythia6CommandConverter.hh"

#include "EvtGenBase/EvtReport.hh"

#include <iostream>
#include <stdlib.h>

using std::endl;

std::vector<std::string> convertPythia6Command( Command command )
{
    std::string module = command["MODULE"];
    std::string param = command["PARAM"];
    std::string value = command["VALUE"];
    std::vector<std::string> commandStrings;
    if ( module == "MSTJ" ) {
        switch ( atoi( param.c_str() ) ) {
            //1,2,3
            case 11:
                switch ( atoi( value.c_str() ) ) {
                    case 3:
                        commandStrings.push_back( "StringZ:usePetersonC = on" );
                        commandStrings.push_back( "StringZ:usePetersonB = on" );
                        commandStrings.push_back( "StringZ:usePetersonH = on" );
                        break;
                    case 1:
                    case 4:
                        EvtGenReport( EVTGEN_ERROR, "EvtGen" )
                            << "Pythia6 parameter: MSTJ(11)=" << value
                            << " is only implicitly supported." << endl;
                        EvtGenReport( EVTGEN_ERROR, "EvtGen" )
                            << "Please use MSTJ(11)=5 and ensure PARJ(46) and PARJ(47) are both set appropriately."
                            << endl;
                        ::abort();
                    case 5:
                        commandStrings.push_back( "StringZ:usePetersonC = off" );
                        commandStrings.push_back( "StringZ:usePetersonB = off" );
                        commandStrings.push_back( "StringZ:usePetersonH = off" );
                        break;
                    default:
                        EvtGenReport( EVTGEN_ERROR, "EvtGen" )
                            << "Pythia6 parameter: MSTJ(11)=" << value
                            << " is not currently supported." << endl;
                        EvtGenReport( EVTGEN_ERROR, "EvtGen" )
                            << "Please use MSTJ(11)=3 or MSTJ(11)=5." << endl;
                        ::abort();
                }
                break;
            case 12:
                switch ( atoi( value.c_str() ) ) {
                    case 2:
                        commandStrings.push_back(
                            "StringFlav:suppressLeadingB = off" );
                        break;
                    case 3:
                        commandStrings.push_back(
                            "StringFlav:suppressLeadingB = on" );
                        break;
                    default:
                        EvtGenReport( EVTGEN_ERROR, "EvtGen" )
                            << "Pythia6 parameter: MSTJ(12)=" << value
                            << " is not currently supported." << endl;
                        EvtGenReport( EVTGEN_ERROR, "EvtGen" )
                            << "Please use MSTJ(12)=2 or MSTJ(12)=3." << endl;
                        ::abort();
                }
                break;
            //13-19,21-24
            case 26:
                switch ( atoi( value.c_str() ) ) {
                    case 0:
                        commandStrings.push_back( "ParticleDecays:mixB = off" );
                        break;
                    case 1:
                    case 2:
                        commandStrings.push_back( "ParticleDecays:mixB = on" );
                        break;
                }
                break;
            //28,38-50
            //51 Inclusion of BE effects - TODO
            case 52:
                switch ( atoi( value.c_str() ) ) {
                    case 9:
                        commandStrings.push_back( "BoseEinstein:Eta = on" );
                        [[fallthrough]];
                    case 7:
                        commandStrings.push_back( "BoseEinstein:Kaon = on" );
                        [[fallthrough]];
                    case 3:
                        commandStrings.push_back( "BoseEinstein:Pion = on" );
                        break;
                    default:
                        EvtGenReport( EVTGEN_ERROR, "EvtGen" )
                            << "Pythia6 parameter: MSTJ(52)=" << value
                            << " is not allowed." << endl;
                        EvtGenReport( EVTGEN_ERROR, "EvtGen" )
                            << "Please select 3,7 or 9." << endl;
                        ::abort();
                }
                break;
            //53-57,91-93,101-121
            default:
                EvtGenReport( EVTGEN_WARNING, "EvtGen" )
                    << "Pythia6 parameter: " << module << "(" << param
                    << ") is not currently supported and will be ignored."
                    << endl;
                EvtGenReport( EVTGEN_WARNING, "EvtGen" )
                    << "A similar Pythia8 parameter may be available." << endl;
        }
    } else if ( module == "PARJ" ) {
        switch ( atoi( param.c_str() ) ) {
            case 1:
                commandStrings.push_back( "StringFlav:probQQtoQ = " + value );
                break;
            case 2:
                commandStrings.push_back( "StringFlav:probStoUD = " + value );
                break;
            case 3:
                commandStrings.push_back( "StringFlav:probSQtoQQ = " + value );
                break;
            case 4:
                commandStrings.push_back( "StringFlav:probQQ1toQQ0 = " + value );
                break;
            case 5:
                commandStrings.push_back( "StringFlav:popcornRate = " + value );
                break;
            case 6:
                commandStrings.push_back( "StringFlav:popcornSpair = " + value );
                break;
            case 7:
                commandStrings.push_back( "StringFlav:popcornSmeson = " + value );
                break;
            //8-10 Advanced popcorn model - can't find these in Pythia8 (unsupported?)
            case 11:
                commandStrings.push_back( "StringFlav:mesonUDvector = " + value );
                break;
            case 12:
                commandStrings.push_back( "StringFlav:mesonSvector = " + value );
                break;
            case 13:
                commandStrings.push_back( "StringFlav:mesonCvector = " + value );
                commandStrings.push_back( "StringFlav:mesonBvector = " + value );
                break;
            case 14:
                commandStrings.push_back( "StringFlav:mesonUDL1S0J1 = " + value );
                commandStrings.push_back( "StringFlav:mesonSL1S0J1 = " + value );
                commandStrings.push_back( "StringFlav:mesonCL1S0J1 = " + value );
                commandStrings.push_back( "StringFlav:mesonBL1S0J1 = " + value );
                break;
            case 15:
                commandStrings.push_back( "StringFlav:mesonUDL1S1J0 = " + value );
                commandStrings.push_back( "StringFlav:mesonSL1S1J0 = " + value );
                commandStrings.push_back( "StringFlav:mesonCL1S1J0 = " + value );
                commandStrings.push_back( "StringFlav:mesonBL1S1J0 = " + value );
                break;
            case 16:
                commandStrings.push_back( "StringFlav:mesonUDL1S1J1 = " + value );
                commandStrings.push_back( "StringFlav:mesonSL1S1J1 = " + value );
                commandStrings.push_back( "StringFlav:mesonCL1S1J1 = " + value );
                commandStrings.push_back( "StringFlav:mesonBL1S1J1 = " + value );
                break;
            case 17:
                commandStrings.push_back( "StringFlav:mesonUDL1S1J2 = " + value );
                commandStrings.push_back( "StringFlav:mesonSL1S1J2 = " + value );
                commandStrings.push_back( "StringFlav:mesonCL1S1J2 = " + value );
                commandStrings.push_back( "StringFlav:mesonBL1S1J2 = " + value );
                break;
            case 18:
                commandStrings.push_back( "StringFlav:decupletSup = " + value );
                break;
            case 19:
                commandStrings.push_back( "StringFlav:lightLeadingBSup = " +
                                          value );
                commandStrings.push_back( "StringFlav:heavyLeadingBSup = " +
                                          value );
                break;
            //21-24 Gaussian PT distributions for primary hadrons - can't find these in Pythia8
            case 25:
                commandStrings.push_back( "StringFlav:etaSup = " + value );
                break;
            case 26:
                commandStrings.push_back( "StringFlav:etaPrimeSup = " + value );
                break;
            //31,32
            case 33:
                commandStrings.push_back( "StringFragmentation:stopMass = " +
                                          value );
                break;
            //34 Stop mass for MSTJ(11)=2 - can't find MSTJ(11)=2 analogue in Pythia 8 so leaving this out too
            //36
            case 37:
                commandStrings.push_back( "StringFragmentation:stopSmear = " +
                                          value );
                break;
            //39,40
            case 41:
                commandStrings.push_back( "StringZ:aLund = " + value );
                break;
            case 42:
                commandStrings.push_back( "StringZ:bLund = " + value );
                break;
            //43,44
            case 45:
                commandStrings.push_back( "StringZ:aExtraDiquark = " + value );
                break;
            case 46:
                commandStrings.push_back( "StringZ:rFactC = " + value );
                break;
            case 47:
                commandStrings.push_back( "StringZ:rFactB = " + value );
                break;
            //48,49,50,51-55,59,61-63,64,65,66,71,72,73,74
            case 76:
                commandStrings.push_back( "ParticleDecays:xBdMix = " + value );
                break;
            case 77:
                commandStrings.push_back( "ParticleDecays:xBsMix = " + value );
                break;
            //80-90 Time-like parton showers - can't find these in Pythia8
            case 91:
                commandStrings.push_back( "BoseEinstein:widthSep = " + value );
                break;
            case 92:
                commandStrings.push_back( "BoseEinstein:lambda = " + value );
                break;
            case 93:
                commandStrings.push_back( "BoseEinstein:QRef = " + value );
                break;
            //94-96 Further BE parameters - can't find these in Pythia8
            //121-171 parameters for ee event generation - can't find these in Pythia8
            //180-195 Various coupling constants & parameters related to couplings - can't find these in Pythia8
            default:
                EvtGenReport( EVTGEN_WARNING, "EvtGen" )
                    << "Pythia6 parameter: " << module << "(" << param
                    << ") is not currently supported and will be ignored."
                    << endl;
                EvtGenReport( EVTGEN_WARNING, "EvtGen" )
                    << "A similar Pythia8 parameter may be available." << endl;
        }
    } else {
        EvtGenReport( EVTGEN_WARNING, "EvtGen" )
            << "Pythia6 parameter: " << module << "(" << param
            << ") is not currently supported and will be ignored." << endl;
        EvtGenReport( EVTGEN_WARNING, "EvtGen" )
            << "A similar Pythia8 parameter may be available." << endl;
    }
    return commandStrings;
}
