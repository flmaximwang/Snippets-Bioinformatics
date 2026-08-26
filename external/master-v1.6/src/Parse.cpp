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

#include "Parse.h"

void Parse::extractNumRes(FILE* & ifp, int & numres, const string & pdsfname)
{
	numres = 0;
	ASSERT(fseek(ifp, sizeof(int) * SEC_NUMRES, SEEK_SET) == 0, "could not seek number of residues in file %s", pdsfname.c_str());	
	ASSERT(fread((void*) &numres, sizeof(int), 1, ifp) == 1, "could not read number of residues in file %s", pdsfname.c_str());
	ASSERT(numres > 0, "invalid number of residues in file %s", pdsfname.c_str());
}

void Parse::extractPdb(FILE* & ifp, vector<string> & pdbinfo, const int & numres, const string & pdsfname)
{
	pdbinfo.clear();
	
	int i, offset, os_pdbinfo, pos, reslen, maxreslen = 10000; // guess as to what the max possible residue string will be
	char* resinfo = (char*) malloc((maxreslen + 1) * sizeof(char)); // 1 for '\0'
	
	ASSERT(fseek(ifp, sizeof(int) * SEC_PDBINFO, SEEK_SET) == 0, "could not seek offset of PDB information section in file %s", pdsfname.c_str());
	ASSERT(fread((void*) &os_pdbinfo, sizeof(int), 1, ifp) == 1, "could not read offset of PDB information section in file %s", pdsfname.c_str());
	ASSERT(fseek(ifp, os_pdbinfo, SEEK_SET) == 0, "could not seek PDB information section in file %s", pdsfname.c_str());
	ASSERT(fread((void*) &offset, sizeof(int), 1, ifp) == 1, "could not read offset of PDB information in file %s", pdsfname.c_str());
	pos = os_pdbinfo + offset;
	ASSERT(fseek(ifp, pos, SEEK_SET) == 0, "could not seek PDB information in file %s", pdsfname.c_str());
	for (i = 0; i < numres; i++)
	{
		ASSERT(fread((void*) &reslen, sizeof(int), 1, ifp) == 1, "could not read length of PDB information in file %s", pdsfname.c_str());
		if (reslen > maxreslen)
		{
			maxreslen = reslen;
			free((void*) resinfo);
			resinfo = (char*) malloc((maxreslen + 1) * sizeof(char)); // 1 for '\0'
		}
		ASSERT(fread((void*) resinfo, sizeof(char), reslen, ifp) == reslen, "could not read PDB information in file %s", pdsfname.c_str());
		resinfo[reslen] = '\0';
		pdbinfo.push_back(string(resinfo));
	}
	free((void*) resinfo);
}

void Parse::extractSeq(FILE* & ifp, const int & numres, char* & seq, const string & pdsfname)
{
	int offset;
	const int seqlen = LEN_AA_CODE * numres;
	seq = new char [seqlen];

	ASSERT(fseek(ifp, sizeof(int) * SEC_SEQ, SEEK_SET) == 0, "could not seek offset of sequence section in file %s", pdsfname.c_str());
	ASSERT(fread((void*) &offset, sizeof(int), 1, ifp) == 1, "could not read offset of sequence section in file %s", pdsfname.c_str());
	ASSERT(fseek(ifp, offset, SEEK_SET) == 0, "could not seek sequence section in file %s", pdsfname.c_str());
	ASSERT(fread((void*) seq, sizeof(char), seqlen, ifp) == seqlen, "could not read sequence section in file %s", pdsfname.c_str());
}

void Parse::outBrk(fstream & ofs, const string & pdbfname, const vector<string> & legalAA)
{
	vector<int> brk;
	int i;
	Structure sys;

	Parse::parsePdb(sys, pdbfname, legalAA);
	findBrk(sys, brk);

	if (ofs.is_open())
	{
		for (i = 0; i < brk.size(); i++)
		{
			ofs << brk[i] << " ";
		}
		ofs << "\n";
	}
	else
	{
		for (i = 0; i < brk.size(); i++)
		{
			cout << brk[i] << " ";
		}
		cout << "\n";
	}
}

void Parse::outNumRes(fstream & ofs, const int & numres)
{
	if (ofs.is_open())
	{
		ofs << numres << "\n";
	}
	else
	{
		cout << numres << "\n";
	}
}

void Parse::outPdb(FILE* & ifp, const int & numres, const string & pdbfname, const string & pdsfname)
{
	int i;
	fstream ofs;
	vector<string> pdbinfo;
	
	extractPdb(ifp, pdbinfo, numres, pdsfname);
	openFileCPP(ofs, pdbfname, (ios::out | ios::trunc));
	for (i = 0; i < pdbinfo.size(); i++)
	{
		ofs << pdbinfo[i];
	}		
	ofs.close();
}

void Parse::outSeq(FILE* & ifp, const int & numres, fstream & ofs, const string & pdsfname)
{
	string aa;
	int i;
	char* seq;
	const int seqlen = LEN_AA_CODE * numres;

	extractSeq(ifp, numres, seq, pdsfname);

	for (i = 0; i < seqlen; i += LEN_AA_CODE)
	{
		aa.append(seq + i, LEN_AA_CODE);
		aa.append(" ");
	}
	aa[aa.size() - 1] = '\n';
	if (ofs.is_open())
	{
		ofs << aa;
	}
	else
	{
		cout << aa;
	}
}

vector<string> Parse::parsePdb(Structure & sys, const Structure & S, const vector<string> & legalAA)
{
	sys.reset();

	vector<string> acceptAA;
	Chain* currChain;

	for (int i = 0; i < S.chainSize(); i++)
	{
		currChain = NULL;

		Chain& chain = S[i];
		for (int j = 0; j < chain.residueSize(); j++)
		{
			Residue& r = chain[j];

			if (!legalAA.empty())
			{
				if (!isProtein(legalAA, r.getName()))
				{
					warning("residue '%s,%d%c' has been skipped due to unidentified residue name '%s'",
						r.getChainID().c_str(), r.getNum(), r.getIcode(), r.getName().c_str());
					continue;
				}

				// If skipping only those residues missing CA atom, when it comes to sorting a list of matches by full-backbone RMSD, 
				// which were originally sorted by CA RMSD, MASTER will get an error when processing matches missing the other backbone atoms.
				// It is important to keep the consistency between CA RMSD based and full-backbone RMSD based searches.
				if (!hasFullBackbone(r))
				{
					warning("residue '%s,%d%c' (name '%s') has been skipped due to missing backbone atom(s)",
						r.getChainID().c_str(), r.getNum(), r.getIcode(), r.getName().c_str());
					continue;
				}
			}
			else
			{
				if (!hasFullBackbone(r))
				{
					error("residue '%s,%d%c' (name '%s') is missing backbone atom(s)",
						r.getChainID().c_str(), r.getNum(), r.getIcode(), r.getName().c_str());
				}
			}

			if (find(acceptAA.begin(), acceptAA.end(), r.getName()) == acceptAA.end())
			{
				acceptAA.push_back(r.getName());
			}

			if (NULL == currChain)
			{
				currChain = new Chain(chain.getID(), chain.getSegID());
			}
			Residue* res = new Residue(r);
			currChain->appendResidue(res);
		}

		if (NULL != currChain)
		{
			sys.appendChain(currChain);
		}
	}

	return acceptAA;
}

vector<string> Parse::parsePdb(Structure & sys, const string & pdbfname, const vector<string> & legalAA)
{
	Structure S(pdbfname);
	return Parse::parsePdb(sys, S, legalAA);
}

