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

#if !defined(TARGETSTRUCT_H)
#define TARGETSTRUCT_H

#include "Common.h"
#include "ProteinStruct.h"

class TargetStruct : public ProteinStruct
{	
	public:
		TargetStruct() : ProteinStruct() {}
		TargetStruct(Structure & S, const string & pdbfname = "target", const bool & clean = false, const bool & nr = false, const double & seqIdenThresh = 0.9, const double & optRmsdThresh = 1.0, const double & contRmsdThresh = 2.0, const bool & nrPdb = false, const double & phistep = 10.0, const double & psistep = 10.0, const double & dcut = 25.0, const double & dstep = 5.0);
		~TargetStruct() {}

		void fillTargetStruct(Structure & S, const string & pdbfname, const bool & clean = false, const bool & nr = false, const double & seqIdenThresh = 0.9, const double & optRmsdThresh = 1.0, const double & contRmsdThresh = 2.0, const bool & nrPdb = false, const double & phistep = 10.0, const double & psistep = 10.0, const double & dcut = 25.0, const double & dstep = 5.0);

		// function members for createPDS only
		TargetStruct(const string &, const vector < string > &, const bool &, const double &, const double &, const double &, const bool &, const double &, const double &, const double &, const double &);

		void checkTargetStructFile(FILE* &, const string &, const CreateOptions &);
		void setDihedDistr(const double &, const double &);
		void setDistDistr(const double &, const double &);
		void setTargetStruct(const double &, const double &, const double &, const double &);
		void writeDihedDistr(fstream &, const CreateOptions &);
		void writeDistDistr(fstream &, const CreateOptions &);
		void writeTargetStruct(fstream &, const CreateOptions &);
		void writeTargetStructFile(fstream &, const CreateOptions &);
		int writeTargetStructFileHeader(fstream &, const CreateOptions &);

		// function members for master only
		TargetStruct(const string &);

		void initDihedDistr();
		void initDistDistr();
		void readDihedDistr(vector<int> &, const double &, const double &, const double &, const double &);
		void readDistDistr(vector<int> &, const int &, const double &, const double &);
		void readTargetStruct(const string &);
	
	protected:
		double _dcut;
		vector<vector<vector<int> > > _diheddistr;
		vector<vector<vector<int> > > _distdistr;
		double _dstep;
		double _phistep;
		double _psistep;
		vector<vector<bool> > _readdiheddistr;
		vector<vector<bool> > _readdistdistr;

		// data members for master only
		int _os_diheddistr;
		int _os_distdistr;
};

#endif
