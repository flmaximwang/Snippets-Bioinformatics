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

#if !defined(QUERYSTRUCT_H)
#define QUERYSTRUCT_H

#include "CenResGrp.h"
#include "Common.h"
#include "ProteinStruct.h"

class QueryStruct : public ProteinStruct
{
	public:
		QueryStruct() : ProteinStruct() {}
		QueryStruct(Structure & S, const string & pdbfname = "query", const bool & clean = false, const double & dcut = 25.0);
		~QueryStruct() {}

		void fillQueryStruct(Structure & S, const string & pdbfname, const bool & clean = false, const double & dcut = 25.0);
		CenResGrp& getCenResGrp() { return _crgrp; }
		double getRmsdCut() const { return _rthresh; }
		int getRmsdMode() const { return _rmode; }
		double getRmsdTune() const { return _tune; }
		void initCenResGrp();
		void setSearchConstraints(const double & rthresh, const vector<vector<int> > & gaplen = vector<vector<int> >(), const int & rmode = 0, const double & tune = 0.5);

		// function members for createPDS only
		QueryStruct(const string &, const vector<string> &, const double &);

		void checkQueryStructFile(FILE* &, const string &, const CreateOptions &);
		void setCenRes(const double &);
		void setCenResDihed();
		void setNumSeg();
		void setQueryStruct(const double &);
		void setResBefBrk();
		void writeCenRes(fstream &, const CreateOptions &);
		void writeCenResDihed(fstream &, const CreateOptions &);
		void writeQueryStruct(fstream &, const CreateOptions &);
		void writeQueryStructFile(fstream &, const CreateOptions &);
		int writeQueryStructFileHeader(fstream &, const CreateOptions &);
		void writeResBefBrk(fstream &, const CreateOptions &);

		// function members for master only
		QueryStruct(const string &);

		int getCenRes(const int & i) const { return _cenres[i]; }
		vector<int>& getCenRes() { return _cenres; }
		vector<vector<double> >& getCenResDihed() { return _crdihed; }
		int getNumSeg() const { return _numseg; }
		int getResBefBrk(const int & i) const { return _befbrk[i]; }
		vector<int>& getResBefBrk() { return _befbrk; }
		void readCenRes();
		void readCenResDihed();
		void readNumSeg();
		void readQueryStruct(const string &);
		void readResBefBrk();

	protected:
		vector<int> _befbrk;
		vector<int> _cenres;
		vector<vector<double> > _crdihed;
		CenResGrp _crgrp;
		int _numseg;
		int _rmode;
		double _rthresh;
		double _tune;
};

#endif
