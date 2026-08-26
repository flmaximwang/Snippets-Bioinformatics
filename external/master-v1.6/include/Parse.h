/*
----------------------------------------------------------------------------
This file is part of MASTER.

MASTER is free software: you can redistribute it and/or modify it under
the terms of the GNU Lesser General Public License as published by the
Free Software Foundation, either version 3 of the License, or (at your
option) any later version.

MASTER is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details.

You should have received a copy of the GNU Lesser General Public License
along with MASTER.  If not, see <http://www.gnu.org/licenses/>.

Copyright (C) 2014 Jianfu Zhou, Gevorg Grigoryan
----------------------------------------------------------------------------
*/

#if !defined(PARSE_H)
#define PARSE_H

#include "Common.h"

class Parse {
	public:
		static void extractNumRes(FILE* &, int &, const string &);
		static void extractPdb(FILE* &, vector<string> &, const int &, const string &);
		static void extractSeq(FILE* &, const int &, char* &, const string &);
		static void outBrk(fstream & ofs, const string & pdbfname, const vector<string> & legalAA = vector<string>());
		static void outNumRes(fstream &, const int &);
		static void outPdb(FILE* &, const int &, const string &, const string &);
		static void outSeq(FILE* &, const int &, fstream &, const string &);
		static vector<string> parsePdb(Structure & sys, const string & pdbfname, const vector<string> & legalAA = vector<string>());
		static vector<string> parsePdb(Structure & sys, const Structure & S, const vector<string> & legalAA = vector<string>());
};

#endif
