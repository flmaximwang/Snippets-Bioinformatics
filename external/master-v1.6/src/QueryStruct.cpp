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

#include "QueryStruct.h"
#include "Search.h"

QueryStruct::QueryStruct(Structure & S, const string & pdbfname, const bool & clean, const double & dcut)
{
	fillQueryStruct(S, pdbfname, clean, dcut);
}

void QueryStruct::fillQueryStruct(Structure & S, const string & pdbfname, const bool & clean, const double & dcut)
{
	parsePdb(S, pdbfname, legalResidueNames(clean));
	setQueryStruct(dcut);
	initCenResGrp();
}

void QueryStruct::initCenResGrp()
{
	int i;

	_crgrp.initCenResGrp(_numseg);
	for (i = 0; i < _numseg; i++)
	{
		_crgrp[i].setResIdx(_cenres[i]);
		_crgrp[i].setSegIdx(i);
		_crgrp[i].setBegRes(_befbrk[i] + 1);
		_crgrp[i].setEndRes(_befbrk[i + 1]);
		_crgrp[i].setSegLen(_befbrk[i + 1] - _befbrk[i]);
		_crgrp[i].setPhi(_crdihed[i][0]);
		_crgrp[i].setPsi(_crdihed[i][1]);
	}
}

void QueryStruct::setSearchConstraints(const double & rthresh, const vector < vector < int > > & gaplen, 
										   const int & rmode, const double & tune)
{
	int i;

	_rthresh = rthresh;
	_rmode = rmode;
	_tune = tune;

	for (i = 0; i < _numseg; i++)
	{
		// assume that all the other segments have perfect matches
		// NOTE: BB and CA share the same RMSD cutoff, and factor 4 (i.e., number of BB atoms) makes no difference
		_crgrp[i].setSelfRmsdCut(rmsdCut(rmode, rthresh, _numres, 0.0, 0, _crgrp[i].getSegLen(), tune));
		_crgrp[i].setRank(0.0);
		if ((i > 0) && (gaplen.size() > 0))
		{
			_crgrp[i].setGapLen(gaplen[i - 1]);
			if (_crgrp[i].getGapLen().size() > 0)
			{
				_crgrp[i - 1].setRank(float(i - 1 - _numseg));
				_crgrp[i].setRank(float(i - _numseg));
			}
		}
	}
}

// function members for createPDS only
QueryStruct::QueryStruct(const string & pdbfname, const vector<string> & legalAA, const double & dcut)
{
	parsePdb(pdbfname, legalAA);
	setQueryStruct(dcut);
}

void QueryStruct::checkQueryStructFile(FILE* & ifp, const string & qsfname,const CreateOptions & copts)
{
	checkProteinStructFile(ifp, qsfname, copts);
}

void QueryStruct::setCenRes(const double & dcut)
{
	_cenres.clear();

	int i, j, k, crind;
	double radius, maxdist, d;
	for (i = 0; i < _numseg; i++)
	{
		radius = MAX_DOUBLE;
		crind = -1;
		for (j = _befbrk[i] + 1; j <= _befbrk[i + 1]; j++)
		{
			maxdist = 0.0;
			for (k = _befbrk[i] + 1; k <= _befbrk[i + 1]; k++)
			{
				d = _pselca[j]->distance(*_pselca[k]);
				if (d > maxdist)
				{
					maxdist = d;
				}						
			}
			if (maxdist < radius)
			{
				radius = maxdist;
				crind = j;
			}
		}
		ASSERT(crind > -1, "could not find a center residue");
		_cenres.push_back(crind);
	}	

	// distance between any pair of center residues should be <= distance cutoff
	// otherwise, not a good query structure
	for (i = 0; i < (_cenres.size() - 1); i++)
	{
		for (j = i + 1; j < _cenres.size(); j++)
		{
			ASSERT(_pselca[_cenres[i]]->distance(*_pselca[_cenres[j]]) <= dcut, "distance between certain pair of center residues > distance cutoff, and %s is not a good query", _fname.c_str());
		}
	}
}

void QueryStruct::setCenResDihed()
{
	_crdihed.clear();

	int i;

	for (i = 0; i < _cenres.size(); i++)
	{
		_crdihed.push_back(_dihed[_cenres[i]]);
	}
}

void QueryStruct::setNumSeg()
{
	_numseg = _befbrk.size() - 1;
}

void QueryStruct::setQueryStruct(const double & dcut)
{
	setProteinStruct();
	setResBefBrk();
	setNumSeg();
	setCenRes(dcut);
	setCenResDihed();

#if defined(DEBUG_QS)
	cout << "residues before breaks: ";
	for (int i = 0; i < _befbrk.size(); i++)
	{
		cout << _befbrk[i] << " ";
	}
	cout << endl;
	cout << "center residues: ";
	for (int i = 0; i < _cenres.size(); i++)
	{
		cout << _cenres[i] << " ";
	}
	cout << endl;
	cout << "number of segments: " << _numseg << endl;
	exit(-1);
#endif	
}

void QueryStruct::setResBefBrk() {
	_befbrk.clear();

	_befbrk.push_back(-1);
	vector<Residue*> residues = _sys.getResidues();
	for (int i = 0; i < (residues.size() - 1); i++) {
		if (hasBreak(*(residues[i]), *(residues[i + 1]))) {
			_befbrk.push_back(i);
		}
	}
	_befbrk.push_back(residues.size() - 1);
}

void QueryStruct::writeCenRes(fstream & ofs,const CreateOptions & copts)
{
	int i;
	if (!copts.getBinary())
	{
		writeString(ofs, "CENRES\n", copts.getBinary());
	}
	for (i = 0; i < _cenres.size(); i++)
	{
		writeDatum(ofs, int(_cenres[i]), copts.getBinary());
		writeString(ofs, copts.getWordTer(), copts.getBinary());
	}
	if (!copts.getBinary())
	{
		writeString(ofs, "END\n", copts.getBinary());
	}
}

void QueryStruct::writeCenResDihed(fstream & ofs,const CreateOptions & copts)
{
	int i;
	if (!copts.getBinary())
	{
		writeString(ofs, "CENTERRESIDUEDIHEDRALANGLES\n", copts.getBinary());
	}
	for (i = 0; i < _crdihed.size(); i++)
	{
		writeDatum(ofs, double(_crdihed[i][PHI_IDX]), copts.getBinary());
		writeString(ofs, copts.getWordSep(), copts.getBinary());
		writeDatum(ofs, double(_crdihed[i][PSI_IDX]), copts.getBinary());
		writeString(ofs, copts.getWordTer(), copts.getBinary());
	}
	if (!copts.getBinary())
	{
		writeString(ofs, "END\n", copts.getBinary());
	}
}

void QueryStruct::writeQueryStruct(fstream & ofs, const CreateOptions & copts)
{
	writeProteinStruct(ofs, copts);
	writeResBefBrk(ofs, copts);
	writeCenRes(ofs, copts);
	writeCenResDihed(ofs, copts);
}

void QueryStruct::writeQueryStructFile(fstream & ofs, const CreateOptions & copts)
{
	int countbytes = writeQueryStructFileHeader(ofs, copts);
	writeQueryStruct(ofs, copts);

	if (copts.getBinary())
	{
		streampos begin, end;
		ofs.seekp(0, ios::beg);
		begin = ofs.tellp();
		ofs.seekp(0, ios::end);
		end = ofs.tellp();
		ASSERT(countbytes == (end - begin), "size of query structure file not consistent");
	}

	cout << "Query structure for " << _fname << " created." << endl;
}

int QueryStruct::writeQueryStructFileHeader(fstream & ofs, const CreateOptions & copts)
{
	int countbytes = writeProteinStructFileHeader(ofs, copts, NUM_SEC_QS);

	if (!copts.getBinary())
	{
		writeDatum(ofs, int(_numseg), copts.getBinary());
		writeString(ofs, copts.getWordTer(), copts.getBinary());
		return -1;
	}

	// offset of residues before breaks section		
	writeDatum(ofs, int(countbytes), copts.getBinary());
	countbytes += sizeof(int) * _befbrk.size();
	// offset of center residues section
	writeDatum(ofs, int(countbytes), copts.getBinary());
	countbytes += sizeof(int) * _cenres.size();
	// offset of center residue dihedral angles section
	writeDatum(ofs, int(countbytes), copts.getBinary());
	countbytes += sizeof(double) * 2 * _cenres.size();
	// number of segments
	writeDatum(ofs, int(_numseg), copts.getBinary());	

	return countbytes;
}

void QueryStruct::writeResBefBrk(fstream & ofs,const CreateOptions & copts)
{
	int i;
	if (!copts.getBinary())
	{
		writeString(ofs, "RESIDUESBEFOREBREAKS\n", copts.getBinary());
	}
	for (i = 0; i < _befbrk.size(); i++)
	{
		writeDatum(ofs, int(_befbrk[i]), copts.getBinary());
		writeString(ofs, copts.getWordTer(), copts.getBinary());
	}
	if (!copts.getBinary())
	{
		writeString(ofs, "END\n", copts.getBinary());
	}
}

// function members for master only
QueryStruct::QueryStruct(const string & qsfname)
{
	readQueryStruct(qsfname);
}

void QueryStruct::readCenRes()
{
	int offset;
	ASSERT(fseek(_ifp, sizeof(int) * SEC_CENRES, SEEK_SET) == 0, "could not seek offset of center residues section in file %s", _fname.c_str());
	ASSERT(fread((void*) &offset, sizeof(int), 1, _ifp) == 1, "could not read offset of center residues section in file %s", _fname.c_str());
	ASSERT(fseek(_ifp, offset, SEEK_SET) == 0, "could not seek center residues section in file %s", _fname.c_str());	
	int* cr = (int*) malloc(_numseg * sizeof(int));
	ASSERT(fread((void*) cr, sizeof(int), _numseg, _ifp) == _numseg, "could not read center residues section in file %s", _fname.c_str());
	_cenres.assign(cr, cr + _numseg);
	free((void*) cr);
}

void QueryStruct::readCenResDihed()
{
	_crdihed.clear();
	int i, offset;
	vector<double> dihed;
	ASSERT(fseek(_ifp, sizeof(int) * SEC_CRDIHED, SEEK_SET) == 0, "could not seek offset of center residue dihedral angles section in file %s", _fname.c_str());
	ASSERT(fread((void*) &offset, sizeof(int), 1, _ifp) == 1, "could not read offset of center residue dihedral angles section in file %s", _fname.c_str());
	ASSERT(fseek(_ifp, offset, SEEK_SET) == 0, "could not seek center residue dihedral angles section in file %s", _fname.c_str());
	double* phipsi = (double*) malloc(_numseg * 2 * sizeof(double));
	ASSERT(fread((void*) phipsi, sizeof(double), (_numseg * 2), _ifp) == (_numseg * 2), "could not read center residue dihedral angles section in file %s", _fname.c_str());
	for (i = 0; i < (_numseg * 2); i += 2)
	{
		dihed.assign(phipsi + i, phipsi + i + 2);
		_crdihed.push_back(dihed);
	}
	free((void*) phipsi);
}

void QueryStruct::readNumSeg()
{
	ASSERT(fseek(_ifp, sizeof(int) * SEC_NUMSEG, SEEK_SET) == 0, "could not seek number of segments in file %s", _fname.c_str());	
	ASSERT(fread((void*) &_numseg, sizeof(int), 1, _ifp) == 1, "could not read number of segments in file %s", _fname.c_str());
}

void QueryStruct::readQueryStruct(const string & qsfname)
{
	readProteinStruct(qsfname);
	readNumSeg();
	readCenRes();
	readCenResDihed();
	readResBefBrk();

#if defined(DEBUG_QS)
	int i;
	cout << "number of segments: " << _numseg << "\n";	
	cout << "center residues: ";
	for (i = 0; i < _cenres.size(); i++)
	{
		cout << " " << _cenres[i];
	}
	cout << "\n";
	cout << "center residue dihedral angles: ";
	for (i = 0; i < _crdihed.size(); i++)
	{
		printf(" (%f,%f)", _crdihed[i][0], _crdihed[i][1]);
	}
	cout << "\n";
	cout << "residues before breaks: ";
	for (i = 0; i < _befbrk.size(); i++)
	{
		cout << " " << _befbrk[i];
	}
	cout << "\n";
	exit(-1);
#endif
}

void QueryStruct::readResBefBrk()
{
	int offset;		
	ASSERT(fseek(_ifp, sizeof(int) * SEC_BEFBRK, SEEK_SET) == 0, "could not seek offset of residues before breaks section in file %s", _fname.c_str());
	ASSERT(fread((void*) &offset, sizeof(int), 1, _ifp) == 1, "could not read offset of residues before breaks section in file %s", _fname.c_str());
	ASSERT(fseek(_ifp, offset, SEEK_SET) == 0, "could not seek residues before breaks section in file %s", _fname.c_str());
	int* rbb = (int*) malloc((_numseg + 1) * sizeof(int));
	ASSERT(fread((void*) rbb, sizeof(int), (_numseg + 1), _ifp) == (_numseg + 1), "could not read indices before breaks section in file %s", _fname.c_str());	
	_befbrk.assign(rbb, rbb + _numseg + 1);
	free((void*) rbb);
}

