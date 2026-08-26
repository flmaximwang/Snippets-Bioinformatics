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

#if !defined(CREATEOPTIONS_H)
#define CREATEOPTIONS_H

#include "Common.h"

class CreateOptions
{
	public:
		CreateOptions();
		~CreateOptions() {}

		bool getBinary() const { return _bin; }
		double getContRmsdThresh() const { return _contRmsdThresh; }
		double getDistCut() const { return _dcut; }
		double getDistStep() const { return _dstep; }
		string getFileExt() const { return _ext; }
		vector<string>& getLegalAA() { return _legalAA; }
		bool getNonRed() const { return _nr; }
		bool getNonRedPdb() const { return _nrPdb; }
		double getOptRmsdThresh() const { return _optRmsdThresh; }
		vector<string>& getPdbFiles() { return _pdbfnames; }
		vector<string>& getPdsFiles() { return _pdsfnames; }
		string getPdsType() const { return _type; }
		double getPhiStep() const { return _phistep; }
		vector<string>& getPostPdbFiles() { return _opdbfnames; }
		double getPsiStep() const { return _psistep; }		
		double getSeqIdenThresh() const { return _seqIdenThresh; }
		string getWordTer() const { return _ter; }
		string getWordSep() const { return _sep; }
				
		void setContRmsdThresh(const char*);
		void setDistCut(const char*);
		void setDistStep(const char*);
		void setLegalAA(const bool &);
		void setNonRed(const bool & f) { _nr = f; }
		void setNonRedPdb(const bool & f) { _nrPdb = f; }
		void setOptRmsdThresh(const char*);
		void setPdbFile(const string &);
		void setPdbFiles(const string &);
		void setPdsFile(const string & fn, const bool & chk = true);
		void setPdsFiles(const string &);
		void setPdsType(const string &);
		void setPhiStep(const char*);
		void setPostPdbFile(const string &);
		void setPostPdbFiles(const string &);
		void setPsiStep(const char*);	
		void setSeqIdenThresh(const char*);

		void usage();

	protected:
		bool _bin;
		double _contRmsdThresh; // --gRMSD
		double _dcut;
		double _dstep;
		string _ext;
		vector<string> _legalAA;
		bool _nr;
		bool _nrPdb;
		vector<string> _opdbfnames;
		double _optRmsdThresh; // --lRMSD
		vector<string> _pdbfnames;
		vector<string> _pdsfnames;
		double _phistep;
		double _psistep;
		string _sep;
		double _seqIdenThresh; // --seqID
		string _ter;
		string _type;
};

#endif
