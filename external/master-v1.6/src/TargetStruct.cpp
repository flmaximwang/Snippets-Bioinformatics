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

#include "TargetStruct.h"

TargetStruct::TargetStruct(Structure & S, const string & pdbfname, const bool & clean, const bool & nr, const double & seqIdenThresh, 
							  const double & optRmsdThresh, const double & contRmsdThresh, const bool & nrPdb, const double & phistep, 
							  const double & psistep, const double & dcut, const double & dstep)
{
	fillTargetStruct(S, pdbfname, clean, nr, seqIdenThresh, optRmsdThresh, contRmsdThresh, nrPdb, 
					 phistep, psistep, dcut, dstep);
}

void TargetStruct::fillTargetStruct(Structure & S, const string & pdbfname, const bool & clean, const bool & nr, const double & seqIdenThresh, 
							  const double & optRmsdThresh, const double & contRmsdThresh, const bool & nrPdb, const double & phistep, 
							  const double & psistep, const double & dcut, const double & dstep)
{
	parsePdb(S, pdbfname, legalResidueNames(clean), nr, seqIdenThresh, optRmsdThresh, contRmsdThresh, nrPdb);
	setTargetStruct(phistep, psistep, dcut, dstep);
}

// function members for createPDS only
TargetStruct::TargetStruct(const string & pdbfname, const vector < string > & legalAA, const bool & nr, const double & seqIdenThresh, const double & optRmsdThresh, const double & contRmsdThresh, const bool & nrPdb, const double & phistep, const double & psistep, const double & dcut, const double & dstep)
{
	parsePdb(pdbfname, legalAA, nr, seqIdenThresh, optRmsdThresh, contRmsdThresh, nrPdb);
	setTargetStruct(phistep, psistep, dcut, dstep);
}

void TargetStruct::checkTargetStructFile(FILE * & ifp,const string & tsfname,const CreateOptions & copts)
{
	checkProteinStructFile(ifp, tsfname, copts);
}

void TargetStruct::setDihedDistr(const double & phistep, const double & psistep)
{
	_diheddistr.clear();
	_readdiheddistr.clear();
	_phistep = phistep;
	_psistep = psistep;

	int i, j, ri, numphibin, numpsibin;
	double phi, psi;

	// (0, step], (step, 2step], (2step, 3step], ...
	numphibin = int(ceil(360.0 / phistep)) + 1; // 1 for IMPOSSIBLE_ANGLE
	numpsibin = int(ceil(360.0 / psistep)) + 1; // 1 for IMPOSSIBLE_ANGLE
	_diheddistr.assign(numphibin, vector<vector<int> >(numpsibin, vector<int>()));
	for (ri = 0; ri < _dihed.size(); ri++)
	{
		phi = _dihed[ri][PHI_IDX];
		psi = _dihed[ri][PSI_IDX];
		if (phi >= IMPOSSIBLE_ANGLE)
		{
			i = numphibin - 1;
		}
		else
		{
			i = int(ceil((phi + 180.0) / phistep)) - 1;
		}		
		if (psi >= IMPOSSIBLE_ANGLE)
		{
			j = numpsibin - 1;
		}
		else
		{
			j = int(ceil((psi + 180.0) / psistep)) - 1;
		}		
		_diheddistr[i][j].push_back(ri);
	}

	_readdiheddistr.assign(numphibin, vector<bool>(numpsibin, true));
}

void TargetStruct::setDistDistr(const double & dcut, const double & dstep)
{
	_distdistr.clear();
	_readdistdistr.clear();
	_dcut = dcut;
	_dstep = dstep;

	int i, j, numbin, bi;
	double d;
			
	// (0, dstep], (dstep, 2dstep], (2dstep, 3dstep], ...
	numbin = int(ceil(dcut / dstep));
	_distdistr.assign(_numres, vector<vector<int> >(numbin, vector<int>()));
	for (i = 0; i < (_numres - 1); i++)
	{
		for (j = i + 1; j < _numres; j++)
		{
			d = _dist[make_pair(i, j)];
			if (d > dcut)
			{
				continue;
			}
			if (d <= 0.0)
			{
				bi = 0;
			}
			else
			{
				bi = int(ceil(d / dstep)) - 1;
			}
			_distdistr[i][bi].push_back(j);
			_distdistr[j][bi].push_back(i);
		}
	}

	_readdistdistr.assign(_numres, vector<bool>(numbin, true));
}

void TargetStruct::setTargetStruct(const double & phistep, const double & psistep, const double & dcut, const double & dstep)
{
	setProteinStruct();
	setDihedDistr(phistep, psistep);
	setDistDistr(dcut, dstep);
}

void TargetStruct::writeDihedDistr(fstream & ofs,const CreateOptions & copts)
{
	int i, j, k, countbytes;
	if (copts.getBinary())
	{
		writeDatum(ofs, double(_phistep), copts.getBinary());
		writeDatum(ofs, double(_psistep), copts.getBinary());
		countbytes = sizeof(double) * 2 + sizeof(int) * _diheddistr.size() * _diheddistr[0].size(); // 2 for phi and psi steps
		for (i = 0; i < _diheddistr.size(); i++)
		{
			for (j = 0; j < _diheddistr[i].size(); j++)
			{
				if (_diheddistr[i][j].size() == 0)
				{
					writeDatum(ofs, int(_diheddistr[i][j].size()), copts.getBinary());
				}
				else
				{
					writeDatum(ofs, int(countbytes), copts.getBinary());
					countbytes += sizeof(int) * (1 + _diheddistr[i][j].size()); // 1 for number of elements
				}
			}
		}
		for (i = 0; i < _diheddistr.size(); i++)
		{
			for (j = 0; j < _diheddistr[i].size(); j++)
			{
				if (_diheddistr[i][j].size() > 0)
				{
					writeDatum(ofs, int(_diheddistr[i][j].size()), copts.getBinary());
					for (k = 0; k < _diheddistr[i][j].size(); k++)
					{
						writeDatum(ofs, int(_diheddistr[i][j][k]), copts.getBinary());
					}
				}
			}
		}
	}
	else
	{
		writeString(ofs, "DIHEDRALANGLESDISTRIBUTION\n", copts.getBinary());
		writeDatum(ofs, double(_phistep), copts.getBinary());
		writeString(ofs, copts.getWordSep(), copts.getBinary());
		writeDatum(ofs, double(_psistep), copts.getBinary());
		writeString(ofs, copts.getWordTer(), copts.getBinary());
		for (i = 0; i < _diheddistr.size(); i++)
		{
			for (j = 0; j < _diheddistr[i].size(); j++)
			{
				writeDatum(ofs, int(_diheddistr[i][j].size()), copts.getBinary());
				writeString(ofs, copts.getWordTer(), copts.getBinary());
				for (k = 0; k < _diheddistr[i][j].size(); k++)
				{
					writeDatum(ofs, int(_diheddistr[i][j][k]), copts.getBinary());
					writeString(ofs, copts.getWordTer(), copts.getBinary());
				}
			}
		}
		writeString(ofs, "END\n", copts.getBinary());
	}
}

void TargetStruct::writeDistDistr(fstream & ofs,const CreateOptions & copts)
{
	int i, j, k, countbytes;
	if (copts.getBinary())
	{
		writeDatum(ofs, double(_dcut), copts.getBinary());
		writeDatum(ofs, double(_dstep), copts.getBinary());
		countbytes = sizeof(double) * 2 + sizeof(int) * _distdistr.size() * _distdistr[0].size(); // 2 for distance cutoff and step
		for (i = 0; i < _distdistr.size(); i++)
		{
			for (j = 0; j < _distdistr[i].size(); j++)
			{
				if (_distdistr[i][j].size() == 0)
				{
					writeDatum(ofs, int(_distdistr[i][j].size()), copts.getBinary());
				}
				else
				{
					writeDatum(ofs, int(countbytes), copts.getBinary());
					countbytes += sizeof(int) * (1 + _distdistr[i][j].size()); // 1 for number of elements
				}
			}
		}
		for (i = 0; i < _distdistr.size(); i++)
		{
			for (j = 0; j < _distdistr[i].size(); j++)
			{
				if (_distdistr[i][j].size() > 0)
				{
					writeDatum(ofs, int(_distdistr[i][j].size()), copts.getBinary());
					for (k = 0; k < _distdistr[i][j].size(); k++)
					{
						writeDatum(ofs, int(_distdistr[i][j][k]), copts.getBinary());
					}
				}
			}
		}
	}
	else
	{
		writeString(ofs, "DISTANCEDISTRIBUTION\n", copts.getBinary());
		writeDatum(ofs, double(_dcut), copts.getBinary());
		writeString(ofs, copts.getWordSep(), copts.getBinary());
		writeDatum(ofs, double(_dstep), copts.getBinary());
		writeString(ofs, copts.getWordTer(), copts.getBinary());
		for (i = 0; i < _distdistr.size(); i++)
		{
			for (j = 0; j < _distdistr[i].size(); j++)
			{
				writeDatum(ofs, int(_distdistr[i][j].size()), copts.getBinary());
				writeString(ofs, copts.getWordTer(), copts.getBinary());
				for (k = 0; k < _distdistr[i][j].size(); k++)
				{
					writeDatum(ofs, int(_distdistr[i][j][k]), copts.getBinary());
					writeString(ofs, copts.getWordTer(), copts.getBinary());
				}
			}
		}
		writeString(ofs, "END\n", copts.getBinary());
	}
}

void TargetStruct::writeTargetStruct(fstream & ofs, const CreateOptions & copts)
{
	writeProteinStruct(ofs, copts);
	writeDihedDistr(ofs, copts);
	writeDistDistr(ofs, copts);
}

void TargetStruct::writeTargetStructFile(fstream & ofs, const CreateOptions & copts)
{
	int countbytes = writeTargetStructFileHeader(ofs, copts);
	writeTargetStruct(ofs, copts);

	if (copts.getBinary())
	{
		streampos begin, end;
		ofs.seekp(0, ios::beg);
		begin = ofs.tellp();
		ofs.seekp(0, ios::end);
		end = ofs.tellp();
		ASSERT(countbytes == (end - begin), "size of target structure file not consistent");
	}

	cout << "Target structure for " << _fname << " created." << endl;	
}

int TargetStruct::writeTargetStructFileHeader(fstream & ofs, const CreateOptions & copts)
{
	int countbytes = writeProteinStructFileHeader(ofs, copts, NUM_SEC_TS);

	if (!copts.getBinary())
	{
		return -1;
	}

	int i, j;
	// offset of dihedral angles distribution section
	writeDatum(ofs, int(countbytes), copts.getBinary());
	countbytes += sizeof(double) * 2 + sizeof(int) * _diheddistr.size() * _diheddistr[0].size(); // 2 for phi and psi steps
	for (i = 0; i < _diheddistr.size(); i++)
	{
		for (j = 0; j < _diheddistr[i].size(); j++)
		{
			if (_diheddistr[i][j].size() > 0)
			{
				countbytes += sizeof(int) * (1 + _diheddistr[i][j].size()); // 1 for number of elements
			}
		}
	}
	// offset of distance distribution section
	writeDatum(ofs, int(countbytes), copts.getBinary());
	countbytes += sizeof(double) * 2 + sizeof(int) * _distdistr.size() * _distdistr[0].size(); // 2 for distance cutoff and step
	for (i = 0; i < _distdistr.size(); i++)
	{
		for (j = 0; j < _distdistr[i].size(); j++)
		{
			if (_distdistr[i][j].size() > 0)
			{
				countbytes += sizeof(int) * (1 + _distdistr[i][j].size()); // 1 for number of elements
			}
		}
	}
	
	return countbytes;
}

// function members for master only
TargetStruct::TargetStruct(const string & tsfname)
{
	readTargetStruct(tsfname);
}

void TargetStruct::initDihedDistr()
{
	ASSERT(fseek(_ifp, sizeof(int) * SEC_DIHEDDISTR, SEEK_SET) == 0, "could not seek offset of dihedral angles distribution section in file %s", _fname.c_str());
	ASSERT(fread((void*) &_os_diheddistr, sizeof(int), 1, _ifp) == 1, "could not read offset of dihedral angles distribution section in file %s", _fname.c_str());
	ASSERT(fseek(_ifp, _os_diheddistr, SEEK_SET) == 0, "could not seek dihedral angles distribution section in file %s", _fname.c_str());
	ASSERT(fread((void*) &_phistep, sizeof(double), 1, _ifp) == 1, "could not read phi step in file %s", _fname.c_str());
	ASSERT(fread((void*) &_psistep, sizeof(double), 1, _ifp) == 1, "could not read psi step in file %s", _fname.c_str());
	int numphibin = int(ceil(360.0 / _phistep)) + 1; // 1 for IMPOSSIBLE_ANGLE
	int numpsibin = int(ceil(360.0 / _psistep)) + 1; // 1 for IMPOSSIBLE_ANGLE
	_diheddistr.assign(numphibin, vector<vector<int> >(numpsibin, vector<int>()));
	_readdiheddistr.assign(numphibin, vector<bool>(numpsibin, false));
}

void TargetStruct::initDistDistr()
{
	ASSERT(fseek(_ifp, sizeof(int) * SEC_DISTDISTR, SEEK_SET) == 0, "could not seek offset of distance distribution section in file %s", _fname.c_str());
	ASSERT(fread((void*) &_os_distdistr, sizeof(int), 1, _ifp) == 1, "could not read offset of distance distribution section in file %s", _fname.c_str());
	ASSERT(fseek(_ifp, _os_distdistr, SEEK_SET) == 0, "could not seek distance distribution section in file %s", _fname.c_str());
	ASSERT(fread((void*) &_dcut, sizeof(double), 1, _ifp) == 1, "could not read distance cutoff in file %s", _fname.c_str());
	ASSERT(fread((void*) &_dstep, sizeof(double), 1, _ifp) == 1, "could not read distance step in file %s", _fname.c_str());
	int numbin = int(ceil(_dcut / _dstep));
	_distdistr.assign(_numres, vector<vector<int> >(numbin, vector<int>()));
	_readdistdistr.assign(_numres, vector<bool>(numbin, false));
}

void TargetStruct::readDihedDistr(vector<int> & mcand, const double & phi, const double & psi, const double & phieps, const double & psieps)
{
	mcand.clear();

	int i;
	if ((phieps >= 180.0) && (psieps >= 180.0))
	{
		for (i = 0; i < _numres; i++)
		{
			mcand.push_back(i);
		}
		return;
	}

	int numphibin = int(ceil(360.0 / _phistep)) + 1; // 1 for IMPOSSIBLE_ANGLE
	int numpsibin = int(ceil(360.0 / _psistep)) + 1; // 1 for IMPOSSIBLE_ANGLE
	vector<int> phibi, psibi; // bin indices
	int j;
	if (phi >= IMPOSSIBLE_ANGLE)
	{
		for (i = 0; i < numphibin; i++)
		{
			phibi.push_back(i);
		}
	}
	else
	{
		int bphibin = int(ceil((phi - phieps + 180.0) / _phistep)) - 1;
		int ephibin = int(ceil((phi + phieps + 180.0) / _phistep)) - 1;
		if (ephibin - bphibin >= numphibin - 2)
		{
			for (i = 0; i < numphibin; i++)
			{
				phibi.push_back(i);
			}
		}
		else
		{
			for (i = 0; (bphibin + i) <= ephibin; i++)
			{
				j = bphibin + i;
				while ((j < 0) || (j > (numphibin - 2)))
				{
					if (j < 0)
					{
						j += (numphibin - 1);
					}
					if (j > (numphibin - 2))
					{
						j -= (numphibin - 1);
					}
				}
				phibi.push_back(j);
			}
			sort(phibi.begin(), phibi.end());
			phibi.push_back(numphibin - 1);
		}
	}	
	if (psi >= IMPOSSIBLE_ANGLE)
	{
		for (i = 0; i < numpsibin; i++)
		{
			psibi.push_back(i);
		}
	}
	else
	{
		int bpsibin = int(ceil((psi - psieps + 180.0) / _psistep)) - 1;
		int epsibin = int(ceil((psi + psieps + 180.0) / _psistep)) - 1;
		if (epsibin - bpsibin >= numpsibin - 2)
		{
			for (i = 0; i < numpsibin; i++)
			{
				psibi.push_back(i);
			}
		}
		else
		{
			for (i = 0; (bpsibin + i) <= epsibin; i++)
			{
				j = bpsibin + i;
				while ((j < 0) || (j > (numpsibin - 2)))
				{
					if (j < 0)
					{
						j += (numpsibin - 1);
					}
					if (j > (numpsibin - 2))
					{
						j -= (numpsibin - 1);
					}
				}
				psibi.push_back(j);
			}
			sort(psibi.begin(), psibi.end());
			psibi.push_back(numpsibin - 1);
		}
	}	

	if ((phibi.size() >= numphibin) && (psibi.size() >= numpsibin))
	{
		for (i = 0; i < _numres; i++)
		{
			mcand.push_back(i);
		}
		return;
	}

	bool f;
	int* mc = (int*) malloc(_numres * sizeof(int));
	int numelem, pos;
	vector<int> offset;

	for (i = 0; i < phibi.size(); i++)
	{
		f = false;
		offset.assign(psibi.size(), 0);

		for (j = 0; j < psibi.size(); j++)
		{
			if (_readdiheddistr[phibi[i]][psibi[j]])
			{
				if (_diheddistr[phibi[i]][psibi[j]].size() > 0)
				{
					mcand.insert(mcand.end(), _diheddistr[phibi[i]][psibi[j]].begin(), _diheddistr[phibi[i]][psibi[j]].end());
				}
				continue;
			}

			pos = _os_diheddistr + sizeof(double) * 2 + sizeof(int) * (phibi[i] * numpsibin + psibi[j]);
			if (pos != ftell(_ifp))
			{
				ASSERT(fseek(_ifp, pos, SEEK_SET) == 0, "could not seek offset of dihedral angles distribution");
			}
			ASSERT(fread((void*) &offset[j], sizeof(int), 1, _ifp) == 1, "could not read offset of dihedral angles distribution");
			if (0 == offset[j])
			{
				_readdiheddistr[phibi[i]][psibi[j]] = true;
			}
			else if (!f)
			{
				f = true;
			}
		}

		if (f)
		{
			for (j = 0; j < psibi.size(); j++)
			{
				if (0 == offset[j])
				{
					continue;
				}

				pos = _os_diheddistr + offset[j];
				if (pos != ftell(_ifp))
				{
					ASSERT(fseek(_ifp, pos, SEEK_SET) == 0, "could not seek dihedral angles distribution");
				}
				ASSERT(fread((void*) &numelem, sizeof(int), 1, _ifp) == 1, "could not read number of elements in dihedral angles distribution");
				ASSERT(fread((void*) mc, sizeof(int), numelem, _ifp) == numelem, "could not read dihedral angles distribution");
				_diheddistr[phibi[i]][psibi[j]].insert(_diheddistr[phibi[i]][psibi[j]].end(), mc, mc + numelem);
				_readdiheddistr[phibi[i]][psibi[j]] = true;
				mcand.insert(mcand.end(), mc, mc + numelem);
			}
		}
	}
	free((void*) mc);
}

void TargetStruct::readDistDistr(vector<int> & mcand, const int & ri, const double & d, const double & deps)
{
	mcand.clear();
	if (d - deps > _dcut)
	{
		return;
	}

	int bbin = int(ceil(((d - deps < _dstep) ? _dstep : (d - deps)) / _dstep)) - 1;
	int ebin = int(ceil(((d + deps > _dcut) ? _dcut : (d + deps)) / _dstep)) - 1;
	bool f = false;
	int i, pos;
	int numbin = int(ceil(_dcut / _dstep));
	vector<int> offset(ebin - bbin + 1, 0);

	for (i = 0; (bbin + i) <= ebin; i++)
	{
		if (_readdistdistr[ri][bbin + i])
		{
			if (_distdistr[ri][bbin + i].size() > 0)
			{
				mcand.insert(mcand.end(), _distdistr[ri][bbin + i].begin(), _distdistr[ri][bbin + i].end());
			}
			continue;
		}

		pos = _os_distdistr + sizeof(double) * 2 + sizeof(int) * (ri * numbin + bbin + i);
		if (pos != ftell(_ifp))
		{
			ASSERT(fseek(_ifp, pos, SEEK_SET) == 0, "could not seek offset of distance distribution");
		}
		ASSERT(fread((void*) &offset[i], sizeof(int), 1, _ifp) == 1, "could not read offset of distance distribution");
		if (0 == offset[i])
		{
			_readdistdistr[ri][bbin + i] = true;
		}
		else if (!f)
		{
			f = true;
		}
	}

	if (f)
	{
		int* mc = (int*) malloc(_numres * sizeof(int));
		int numelem;
		for (i = 0; (bbin + i) <= ebin; i++)
		{
			if (0 == offset[i])
			{
				continue;
			}

			pos = _os_distdistr + offset[i];
			if (pos != ftell(_ifp))
			{
				ASSERT(fseek(_ifp, pos, SEEK_SET) == 0, "could not seek distance distribution");
			}
			ASSERT(fread((void*) &numelem, sizeof(int), 1, _ifp) == 1, "could not read number of elements in distance distribution");
			ASSERT(fread((void*) mc, sizeof(int), numelem, _ifp) == numelem, "could not read distance distribution");
			_distdistr[ri][bbin + i].insert(_distdistr[ri][bbin + i].end(), mc, mc + numelem);
			_readdistdistr[ri][bbin + i] = true;
			mcand.insert(mcand.end(), mc, mc + numelem);
		}
		free((void*) mc);
	}
}

void TargetStruct::readTargetStruct(const string & tsfname)
{
	readProteinStruct(tsfname);
	initDihedDistr();
	initDistDistr();
}

