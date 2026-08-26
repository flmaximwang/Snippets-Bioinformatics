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

#if !defined(PARSEOPTIONS_H)
#define PARSEOPTIONS_H

#include "Common.h"

class ParseOptions
{
	public:
		ParseOptions() {}
		~ParseOptions() {}

		string getBrkFile() const { return _brkfname; }
		string getNumResFile() const { return _nrfname; }
		vector<string>& getPdbFiles() { return _pdbfnames; }
		vector<string>& getPdsFiles() { return _pdsfnames; }
		string getSeqFile() const { return _seqfname; }

		void setBrkFile(const string & fn) { _brkfname = fn; }
		void setNumResFile(const string & fn) { _nrfname = fn; }
		void setPdbFile(const string &);
		void setPdbFiles(const string &);
		void setPdsFile(const string &);
		void setPdsFiles(const string &);
		void setSeqFile(const string & fn) { _seqfname = fn; }

		void usage();

	protected:
		string _brkfname;
		string _nrfname;
		vector<string> _pdbfnames;
		vector<string> _pdsfnames;
		string _seqfname;
};

#endif
