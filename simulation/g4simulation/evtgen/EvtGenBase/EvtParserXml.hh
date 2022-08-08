
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

#ifndef EVTPARSERXML_HH
#define EVTPARSERXML_HH

#include <fstream>
#include <string>
#include <vector>

class EvtParserXml final {
  public:
    bool open( std::string filename );
    bool close();

    bool readNextTag();

    std::string getTagTitle() { return _tagTitle; }
    std::string getParentTagTitle();
    int getLineNumber() { return _lineNo; }
    bool isTagInline() { return _inLineTag; }

    std::string readAttribute( std::string attribute,
                               std::string defaultValue = "" );
    bool readAttributeBool( std::string attribute, bool defaultValue = false );
    int readAttributeInt( std::string attribute, int defaultValue = -1 );
    double readAttributeDouble( std::string attribute, double defaultValue = -1. );

  private:
    std::ifstream _fin;
    std::string _line;
    int _lineNo = 0;

    std::string _tag;
    std::string _tagTitle;
    bool _inLineTag;
    std::vector<std::string> _tagTree;

    bool processTagTree();

    bool expandEnvVars( std::string& str );
    bool isAlphaNum( char c );
};

#endif
