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

#if !defined(PROTEINSTRUCT_H)
#define PROTEINSTRUCT_H

#include "Common.h"
#include "CreateOptions.h"
#include "Redundancy.h"
#include "Parse.h"

class ProteinStruct
{	
	public:
		ProteinStruct() { _ifp = NULL; }
		~ProteinStruct() {}

		vector<string> legalResidueNames(const bool & clean = false);
		void parsePdb(Structure & S0, const string & pdbfname, const vector<string> & legalAA, const bool & nr = false, const double & seqIdenThresh = 0.9, const double & optRmsdThresh = 1.0, const double & contRmsdThresh = 2.0, const bool & nrPdb = false);

		// function members for createPDS only
		ProteinStruct(const string & pdbfname, const vector<string> & legalAA, const bool & nr = false, const double & seqIdenThresh = 0.9, const double & optRmsdThresh = 1.0, const double & contRmsdThresh = 2.0, const bool & nrPdb = false);

		void checkBBCoords(FILE* &, const string &, const CreateOptions &);
		void checkCACoords(FILE* &, const string &, const CreateOptions &);
		void checkProteinStructFile(FILE* &, const string &, const CreateOptions &);
		Structure& getProteinSys() { return _sys; }
		void parsePdb(const string & pdbfname, const vector<string> & legalAA, const bool & nr = false, const double & seqIdenThresh = 0.9, const double & optRmsdThresh = 1.0, const double & contRmsdThresh = 2.0, const bool & nrPdb = false);
		void setBBCoords();
		void setCACoords();
		void setDihed();
		void setDist();
		void setFileName(const string & pdbfname) { _fname = pdbfname; }
		void setNumRes() { _numres = _sys.positionSize(); }
		void setPdbInfo();
		void setProteinStruct();
		void setSeq();
		void writeBBCoords(fstream &, const CreateOptions &);
		void writeCACoords(fstream &, const CreateOptions &);
//			void writeDihedralAngles(fstream &, const CreateOptions &);
//			void writeDistance(fstream &, const CreateOptions &);
		void writePdbInfo(fstream &, const CreateOptions &);
		void writeProteinStruct(fstream &, const CreateOptions &);
		void writeProteinStructFile(fstream &, const CreateOptions &);
		int writeProteinStructFileHeader(fstream &, const CreateOptions &, const int &);
		void writeSeq(fstream &, const CreateOptions &);

		// function members for master only
		ProteinStruct(const string &);

		double calcDistBB(const int &, const int &);
		double calcDistCA(const int &, const int &);
		void closeProteinStructFile();
		bool fullBB();
		vector<vector<double> >& getBBCoords() { return _bbcoor; }
		vector<vector<double> >& getCACoords() { return _cacoor; }
		string getFileName() { return _fname; }
		bool getFullBB(const int & i) const { return _fullbb[i]; }
		vector<bool>& getFullBB() { return _fullbb; }
		int getNumRes() const { return _numres; }
		vector<string>& getPdbInfo() { return _pdbinfo; }	
		void initPdbInfo();
		void readBBCoords();
		void readCACoords();
		void readNumRes();
		void readPdbInfo();
		void readPdbInfo(const vector<int> &);
		void readProteinStruct(const string &);
		void readSeq(char* &);

	protected:
		vector<vector<double> > _bbcoor;
		vector<vector<double> > _cacoor;
		string _fname;
		vector<bool> _fullbb;
		int _numres;
		vector<string> _pdbinfo;

		// data members for createPDS only
		vector<vector<double> > _dihed;
		map<pair<int, int>, double> _dist;
		AtomPointerVector _pselca;
		vector<string> _seq;
		Structure _sys;

		// data members for master only
		FILE* _ifp;
};

#endif
