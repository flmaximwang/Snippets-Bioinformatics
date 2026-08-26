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

#include "ProteinStruct.h"

vector<string> ProteinStruct::legalResidueNames(const bool & clean)
{
	vector<string> legalAA = vector<string>();

	if (clean)
	{
		legalAA.push_back("ALA");
		legalAA.push_back("ARG");
		legalAA.push_back("ASN");
		legalAA.push_back("ASP");
		legalAA.push_back("CYS");
		legalAA.push_back("GLN");
		legalAA.push_back("GLU");
		legalAA.push_back("GLY");
		legalAA.push_back("HIS"); legalAA.push_back("HSC"); legalAA.push_back("HSD"); legalAA.push_back("HSE"); legalAA.push_back("HSP");
		legalAA.push_back("ILE");
		legalAA.push_back("LEU");
		legalAA.push_back("LYS");
		legalAA.push_back("MET"); legalAA.push_back("MSE");
		legalAA.push_back("PHE");
		legalAA.push_back("PRO");
		legalAA.push_back("SER");
		legalAA.push_back("THR");
		legalAA.push_back("TRP");
		legalAA.push_back("TYR");
		legalAA.push_back("VAL");

		legalAA.push_back("CSO");
		legalAA.push_back("HIP");
		legalAA.push_back("PTR");
		legalAA.push_back("SEC");
		legalAA.push_back("SEP");
		legalAA.push_back("TPO");
	}

	return legalAA;
}

void ProteinStruct::parsePdb(Structure & S0, const string & pdbfname, const vector<string> & legalAA, const bool & nr, const double & seqIdenThresh, const double & optRmsdThresh, const double & contRmsdThresh, const bool & nrPdb)
{
	_sys.reset();

	setFileName(pdbfname);

	Structure S;
	vector<string> acceptAA = Parse::parsePdb(S, S0, legalAA);

	if (nr) {
		cout << "Removing redundant chains...\n";

		map<Chain*, AtomPointerVector> chainCa = getChainCa(S);
		map<Chain*, vector<double> > chainCaCentroid = getChainCaCentroid(chainCa);
		map<Chain*, vector<Chain*> > chainContacts = getChainContacts(S);
		map<Chain*, int> chainGrpIdx;
		vector<vector<Chain*> > chainGrps;
		vector<vector<Chain*> > chainGrpsRefined;
		map<Chain*, bool> chainKept;
		map<Chain*, vector<string> > chainSeq = getChainSeq(S);
		map<Chain*, bool> chainRemoved;
		for (int i = 0; i < S.chainSize(); i++)
		{
			chainRemoved[&(S.getChain(i))] = false; // True: can be removed safely
													// False: can't be removed
		}
		Chain* chn;
		vector<CorrChains> coChains;
		vector<Chain*> curChainSet, mostOvlpChainSet;
		bool f;
		int i, j, k, minChainIdx, minGrpIdx, minNumChains, ttlNumChains;
		map<Chain*, vector<Chain*> >::iterator it, itA, itB;
		vector<double> t(3), ct(3);
		vector<vector<double> > u(3, vector<double>(3)), cu(3, vector<double>(3));

		for (itA = chainContacts.begin(); itA != chainContacts.end(); ++itA) {
			if (chainGrpIdx.find(itA->first) == chainGrpIdx.end()) {
				chainGrps.push_back(vector<Chain*>());
				chainGrps[chainGrps.size() - 1].push_back(itA->first);
				chainGrpIdx[itA->first] = chainGrps.size() - 1;
			}

			it = itA;
			++it;
			for (itB = it; itB != chainContacts.end(); ++itB) {
				if ((chainGrpIdx.find(itA->first) != chainGrpIdx.end())
					&& (chainGrpIdx.find(itB->first) != chainGrpIdx.end())
					&& (chainGrpIdx[itA->first] == chainGrpIdx[itB->first])) {
					continue;
				}

				if (alignChains(acceptAA, chainSeq[itA->first], chainSeq[itB->first], seqIdenThresh, chainCa[itA->first], chainCa[itB->first], 1, optRmsdThresh, ct, cu, contRmsdThresh, t, u)) {
					// two central chains are redundant
					coChains = getCorrChains(itA->second, itB->second, chainCaCentroid, ct, cu);
					f = true;
					for (i = 0; i < coChains.size(); i++) {
						if (!alignChains(acceptAA, chainSeq[coChains[i].chnA], chainSeq[coChains[i].chnB], seqIdenThresh, chainCa[coChains[i].chnA], chainCa[coChains[i].chnB], 0, optRmsdThresh, t, u, contRmsdThresh, ct, cu)) {
							// two corresponding neighbor chains are not redundant
							f = false;
							break;
						}
					}
					if (f) {
						// chain A plus its neighbors and chain B plus its neighbors are redundant
						if (itA->second.size() < itB->second.size()) {
							// chain B plus its neighbors contains chain A plus its neighbors
							// therefore safely remove chain A and its neighbors	
							chainRemoved[itA->first] = true;
						} else if (itA->second.size() > itB->second.size()) {
							// chain A plus its neighbors contains chain B plus its neighbors
							// therefore safely remove chain B and its neighbors
							chainRemoved[itB->first] = true;
						} else if ((chainGrpIdx.find(itA->first) != chainGrpIdx.end()) && (chainGrpIdx.find(itB->first) != chainGrpIdx.end())) {
							// chain A plus its neighbors == chain B plus its neighbors
							// A and B are in different chain groups
							// merge the two chain groups
							mergeChainGroups(itA->first, itB->first, chainGrpIdx, chainGrps);
						} else if ((chainGrpIdx.find(itA->first) != chainGrpIdx.end()) && (chainGrpIdx.find(itB->first) == chainGrpIdx.end())) {
							// chain A plus its neighbors == chain B plus its neighbors
							// A is in a chain group, while B is not in any chain group
							// add chain B to chain A's group
							chainGrps[chainGrpIdx[itA->first]].push_back(itB->first);
							chainGrpIdx[itB->first] = chainGrpIdx[itA->first];
						} else if ((chainGrpIdx.find(itA->first) == chainGrpIdx.end()) && (chainGrpIdx.find(itB->first) != chainGrpIdx.end())) {
							// chain A plus its neighbors == chain B plus its neighbors
							// B is in a chain group, while A is not in any chain group
							// add chain A to chain B's group
							chainGrps[chainGrpIdx[itB->first]].push_back(itA->first);
							chainGrpIdx[itA->first] = chainGrpIdx[itB->first];
						} else if ((chainGrpIdx.find(itA->first) == chainGrpIdx.end()) && (chainGrpIdx.find(itB->first) == chainGrpIdx.end())) {
							// chain A plus its neighbors == chain B plus its neighbors
							// Neither A nor B are in any chain group
							// assign a new group for chain A and B
							chainGrps.push_back(vector<Chain*>());
							chainGrps[chainGrps.size() - 1].push_back(itA->first);
							chainGrps[chainGrps.size() - 1].push_back(itB->first);
							chainGrpIdx[itA->first] = chainGrps.size() - 1;
							chainGrpIdx[itB->first] = chainGrps.size() - 1;
						}
					}
				}
				// 1. two central chains are not redundant
				// 2. two corresponding neighbor chains are not redundant
				// 3. chain A plus its neighbors contains chain B plus its neighbors, or vice versa
				if (chainGrpIdx.find(itB->first) == chainGrpIdx.end()) {
					chainGrps.push_back(vector<Chain*>());
					chainGrps[chainGrps.size() - 1].push_back(itB->first);
					chainGrpIdx[itB->first] = chainGrps.size() - 1;
				}
			}
		}

		for (i = 0; i < chainGrps.size(); i++) {
			if (chainGrps[i].size() <= 0) {
				continue;
			}

			f = false;
			for (j = 0; j < chainGrps[i].size(); j++) {
				if (chainRemoved[chainGrps[i][j]]) {
					f = true;
					break;
				}
			}
			if (f) {
				for (j = 0; j < chainGrps[i].size(); j++) {
					chainRemoved[chainGrps[i][j]] = true;

#if defined(DEBUG_NR)
					chn = chainGrps[i][j];
					Structure ac;
					ac.appendChain(chn);
					for (k = 0; k < chainContacts[chn].size(); k++) {
						ac.appendChain(chainContacts[chn][k]);
					}
					ac.writePDB(fileName(pdbfname, false, true) + "." + MstUtils::toString(i + 1) + "." + MstUtils::toString(j + 1) + ".r.pdb");
#endif				
				}
				chainGrps[i].clear();
			}
		}

		minNumChains = S.chainSize() + 1;
		minGrpIdx = -1;
		for (i = 0; i < chainGrps.size(); i++) {
			if (chainGrps[i].size() <= 0) {
				continue;
			}

#if defined(DEBUG_NR)
			for (j = 0; j < chainGrps[i].size(); j++) {
				chn = chainGrps[i][j];
				Structure ac;
				ac.appendChain(chn);
				for (k = 0; k < chainContacts[chn].size(); k++) {
					ac.appendChain(chainContacts[chn][k]);
				}
				ac.writePDB(fileName(pdbfname, false, true) + "." + MstUtils::toString(i + 1) + "." + MstUtils::toString(j + 1) + ".k.pdb");
			}
#endif

			chainGrpsRefined.push_back(chainGrps[i]);
			if (chainGrpsRefined[chainGrpsRefined.size() - 1].size() < minNumChains) {
				minNumChains = chainGrpsRefined[chainGrpsRefined.size() - 1].size();
				minGrpIdx = chainGrpsRefined.size() - 1;
			}
		}
		ASSERT((minNumChains > 0) && (minNumChains <= S.chainSize()) 
			&& (minGrpIdx >= 0) && (minGrpIdx < chainGrpsRefined.size()), 
			"could not find group with minimum number of chains");

#if defined(DEBUG_NR)
		printf("group %d contains %d chains, which is minimum.\n", minGrpIdx + 1, minNumChains);
		for (i = 0; i < chainGrpsRefined.size(); i++) {
			printf("group %d contains %d chains.\n", i + 1, chainGrpsRefined[i].size());
		}
		exit(-1);
#endif

		minNumChains = S.chainSize() + 1;
		minChainIdx = -1;
		for (i = 0; i < chainGrpsRefined[minGrpIdx].size(); i++) {
			ttlNumChains = findMostOvlpChainSet(minGrpIdx, i, chainGrpsRefined, chainContacts, curChainSet);
			if (ttlNumChains < minNumChains)
			{
				minNumChains = ttlNumChains;
				minChainIdx = i;
				mostOvlpChainSet = curChainSet;
			}
		}
		ASSERT((minNumChains > 0) && (minNumChains <= S.chainSize()) 
			&& (minChainIdx >= 0) && (minChainIdx < chainGrpsRefined[minGrpIdx].size()), 
			"could not find most overlapping chain set");

		for (i = 0; i < mostOvlpChainSet.size(); i++) {
			chn = mostOvlpChainSet[i];

			if (nrPdb) {
				Structure ac;
				ac.appendChain(chn);
				for (j = 0; j < chainContacts[chn].size(); j++) {
					ac.appendChain(chainContacts[chn][j]);
				}
				ac.writePDB(fileName(pdbfname, false, true) + "." + MstUtils::toString(i + 1) + ".f.pdb");
			}

			chainKept[chn] = true;
			for (j = 0; j < chainContacts[chn].size(); j++) {
				chainKept[chainContacts[chn][j]] = true;
			}
		}
		for (i = 0; i < S.chainSize(); i++) {
			if (chainKept.find(&(S.getChain(i))) == chainKept.end()) {
				continue;
			}
			Chain* currChn = new Chain(S[i]);
			_sys.appendChain(currChn);
		}
	} else {
		_sys = S;
	}

	vector<Atom*> _sysAtoms = _sys.getAtoms();
	_pselca.resize(0);
	for (int i = 0; i < _sysAtoms.size(); i++) {
		if (_sysAtoms[i]->isNamed("CA")) _pselca.push_back(_sysAtoms[i]);
	}
	ASSERT(_pselca.size() == _sys.positionSize(), "numbers of CA atoms and residues not consistent in file %s", pdbfname.c_str());
	ASSERT(_pselca.size() > 0, "no valid residues in file %s", pdbfname.c_str());

#if defined(DEBUG_PARSEPDB)
	cout << "file name: " << _fname << endl;
	cout << "full atoms:" << endl;
	cout << (AtomPointerVector) _sys.getAtoms() << endl;
	cout << "CA atoms:" << endl;
	cout << _pselca << endl;
	exit(-1);
#endif
}

// function members for createPDS only
ProteinStruct::ProteinStruct(const string & pdbfname, const vector<string> & legalAA, const bool & nr, const double & seqIdenThresh, const double & optRmsdThresh, const double & contRmsdThresh, const bool & nrPdb)
{
	parsePdb(pdbfname, legalAA, nr, seqIdenThresh, optRmsdThresh, contRmsdThresh, nrPdb);
	setProteinStruct();
}

void ProteinStruct::checkBBCoords(FILE * & ifp, const string & psfname, const CreateOptions & copts)
{
	int i, j, numres, coordlen, offset, pos;
	double *xyz = NULL;
	vector<vector<double> > bbcoor;

	if (copts.getBinary())
	{
		ASSERT(fseek(ifp, sizeof(int) * SEC_NUMRES, SEEK_SET) == 0, "could not seek number of residues in file %s", psfname.c_str());
		ASSERT(fread((void*) &numres, sizeof(int), 1, ifp) == 1, "could not read number of residues in file %s", psfname.c_str());
		ASSERT(numres == _numres, "number of residues not consistent in file %s", psfname.c_str());

		bbcoor.assign(numres * NUM_BBA, vector<double>());
		coordlen = numres * NUM_BBA * NUM_COORDS;
		xyz = (double*) malloc(coordlen * sizeof(double));
		ASSERT(fseek(ifp, sizeof(int) * SEC_BBCOOR, SEEK_SET) == 0, "could not seek offset of BB coords section in file %s", psfname.c_str());
		ASSERT(fread((void*) &offset, sizeof(int), 1, ifp) == 1, "could not read offset of BB coords section in file %s", psfname.c_str());
		ASSERT(fseek(ifp, offset, SEEK_SET) == 0, "could not seek BB coords section in file %s", psfname.c_str());
		ASSERT(fread((void*) xyz, sizeof(double), coordlen, ifp) == coordlen, "could not read BB coords section in file %s", psfname.c_str());
		for (i = 0; i < numres; i++)
		{
			for (j = 0; j < NUM_BBA; j++)
			{
				pos = i * NUM_BBA * NUM_COORDS + j * NUM_COORDS;
				bbcoor[i * NUM_BBA + j].assign(xyz + pos, xyz + pos + NUM_COORDS);
			}
		}

		ASSERT(bbcoor.size() == _bbcoor.size(), "number of BB coords not consistent in file %s", psfname.c_str());
		for (i = 0; i < bbcoor.size(); i++)
		{
			ASSERT((_bbcoor[i][0] == bbcoor[i][0]) && (_bbcoor[i][1] == bbcoor[i][1]) && (_bbcoor[i][2] == bbcoor[i][2]), "BB coords section not consistent in file %s", psfname.c_str());
		}

		free((void*) xyz);
	}
	else
	{
		error("could not check the correctness of BB coords in text file %s", psfname.c_str());
	}
	
	cout << "BB coords section in file " << psfname << " checked." << endl;
}

void ProteinStruct::checkCACoords(FILE * & ifp, const string & psfname, const CreateOptions & copts)
{
	int i, numres, coordlen, offset, pos;
	double *xyz = NULL;
	vector<vector<double> > cacoor;

	if (copts.getBinary())
	{
		ASSERT(fseek(ifp, sizeof(int) * SEC_NUMRES, SEEK_SET) == 0, "could not seek number of residues in file %s", psfname.c_str());
		ASSERT(fread((void*) &numres, sizeof(int), 1, ifp) == 1, "could not read number of residues in file %s", psfname.c_str());
		ASSERT(numres == _numres, "number of residues not consistent in file %s", psfname.c_str());

		cacoor.assign(numres, vector<double>());
		coordlen = numres * NUM_COORDS;
		xyz = (double*) malloc(coordlen * sizeof(double));
		ASSERT(fseek(ifp, sizeof(int) * SEC_CACOOR, SEEK_SET) == 0, "could not seek offset of CA coords section in file %s", psfname.c_str());
		ASSERT(fread((void*) &offset, sizeof(int), 1, ifp) == 1, "could not read offset of CA coords section in file %s", psfname.c_str());
		ASSERT(fseek(ifp, offset, SEEK_SET) == 0, "could not seek CA coords section in file %s", psfname.c_str());
		ASSERT(fread((void*) xyz, sizeof(double), coordlen, ifp) == coordlen, "could not read CA coords section in file %s", psfname.c_str());
		for (i = 0; i < numres; i++)
		{
			pos = i * NUM_COORDS;
			cacoor[i].assign(xyz + pos, xyz + pos + NUM_COORDS);
		}

		ASSERT(cacoor.size() == _cacoor.size(), "number of CA coords not consistent in file %s", psfname.c_str());
		for (i = 0; i < cacoor.size(); i++)
		{
			ASSERT((_cacoor[i][0] == cacoor[i][0]) && (_cacoor[i][1] == cacoor[i][1]) && (_cacoor[i][2] == cacoor[i][2]), "CA coords section not consistent in file %s", psfname.c_str());
		}

		free((void*) xyz);
	}
	else
	{
		error("could not check the correctness of CA coords in text file %s", psfname.c_str());
	}
	
	cout << "CA coords section in file " << psfname << " checked." << endl;
}

void ProteinStruct::checkProteinStructFile(FILE* & ifp, const string & psfname, const CreateOptions & copts)
{
	checkBBCoords(ifp, psfname, copts);
	checkCACoords(ifp, psfname, copts);
}

void ProteinStruct::parsePdb(const string & pdbfname, const vector < string > & legalAA, const bool & nr, const double & seqIdenThresh, const double & optRmsdThresh, const double & contRmsdThresh, const bool & nrPdb)
{
	Structure S0(pdbfname);
	parsePdb(S0, pdbfname, legalAA, nr, seqIdenThresh, optRmsdThresh, contRmsdThresh, nrPdb);
}

void ProteinStruct::setBBCoords() {
	_bbcoor.assign(NUM_BBA * _sys.positionSize(), vector<double>(NUM_COORDS));
	_fullbb.assign(_sys.residueSize(), true);

	Atom* a;
	vector<Residue*> residues = _sys.getResidues();
	double x, y, z;

	for (int i = 0; i < residues.size(); i++) {
		Residue& r = *(residues[i]);

		// find the four backbone atoms: N, CA, C, and O (some could be named several ways)
		for (int k = 0; k < 4; k++) {
			vector<string> atomNames;
			int atomIndex;
			bool required = false;

			switch (k) {
				case 0:
					atomNames.push_back("N"); atomNames.push_back("NT");
					atomIndex = N_IDX;
					break;

				case 1:
					atomNames.push_back("CA");
					atomIndex = CA_IDX;
					required = true;
					break;

				case 2:
					atomNames.push_back("C");
					atomIndex = C_IDX;
					break;

				case 3:
					atomNames.push_back("O"); atomNames.push_back("OT1"); atomNames.push_back("OT2"); atomNames.push_back("OXT");
					atomIndex = O_IDX;
					break;

				default:
					MstUtils::error("something weird happened", "ProteinStruct::setBBCoords");
			}

			for (int m = 0; m < atomNames.size(); m++) {
				a = NULL;
				a = r.findAtom(atomNames[m], false);
				if (a != NULL) break;
			}

			if (a == NULL) {
				if (required) {
					MstUtils::error("CA atom should have existed", "ProteinStruct::setBBCoords");
				}

				x = y = z = IMPOSSIBLE_COORD;
				_fullbb[i] = false;
			} else {
				x = a->getX(); y = a->getY(); z = a->getZ();
			}

			_bbcoor[NUM_BBA * i + atomIndex][0] = x;
			_bbcoor[NUM_BBA * i + atomIndex][1] = y;
			_bbcoor[NUM_BBA * i + atomIndex][2] = z;
		}
	}

#if defined(DEBUG_BB)
	for (int i = 0; i < residues.size(); i++) {
		cout << endl;
		for (int j = 0; j < residues[i]->atomSize(); j++) {
			cout << (*residues[i])[j].pdbLine() << endl;
		}
		printf("N\t%.3f\t%.3f\t%.3f\n", _bbcoor[4 * i + 0][0], _bbcoor[4 * i + 0][1], _bbcoor[4 * i + 0][2]);
		printf("CA\t%.3f\t%.3f\t%.3f\n", _bbcoor[4 * i + 1][0], _bbcoor[4 * i + 1][1], _bbcoor[4 * i + 1][2]);
		printf("C\t%.3f\t%.3f\t%.3f\n", _bbcoor[4 * i + 2][0], _bbcoor[4 * i + 2][1], _bbcoor[4 * i + 2][2]);
		printf("O\t%.3f\t%.3f\t%.3f\n", _bbcoor[4 * i + 3][0], _bbcoor[4 * i + 3][1], _bbcoor[4 * i + 3][2]);
	}
	exit(-1);
#endif
}

void ProteinStruct::setCACoords()
{
	_cacoor.assign(_pselca.size(), vector<double>(NUM_COORDS));
	for (int i = 0; i < _pselca.size(); i++)
	{
		_cacoor[i][0] = _pselca[i]->getX();
		_cacoor[i][1] = _pselca[i]->getY();
		_cacoor[i][2] = _pselca[i]->getZ();
	}

#if defined(DEBUG_CA)
	for (int i = 0; i < _cacoor.size(); i++)
	{
		cout << _pselca[i]->pdbLine() << endl;
		printf("%f %f %f\n", _cacoor[i][0], _cacoor[i][1], _cacoor[i][2]);
	}
	exit(-1);
#endif
}

void ProteinStruct::setDihed() {
	_dihed.clear();

	int i;
	double phi, psi;
	vector<Residue*> residues = _sys.getResidues();

	for (i = 0; i < residues.size(); i++) {
		phi = IMPOSSIBLE_ANGLE;
		psi = IMPOSSIBLE_ANGLE;

		if (i > 0) {
			Residue& r0 = *(residues[i - 1]);
			Residue& r1 = *(residues[i]);
			if ((!hasBreak(r0, r1)) && (r0.atomExists("C")) && (r1.atomExists("CA")) && (r1.atomExists("C"))) {
				if (r1.atomExists("N")) {
					phi = CartesianGeometry::dihedral(r0.findAtom("C", true), r1.findAtom("N", true), r1.findAtom("CA", true), r1.findAtom("C", true));
				} else if (r1.atomExists("NT")) {
					phi = CartesianGeometry::dihedral(r0.findAtom("C", true), r1.findAtom("NT", true), r1.findAtom("CA", true), r1.findAtom("C", true));
				}
			}
		}
		if (i < (_sys.positionSize() - 1)) {
			Residue& r1 = *(residues[i]);
			Residue& r2 = *(residues[i + 1]);
			if ((!hasBreak(r1, r2)) && (r1.atomExists("CA")) && (r1.atomExists("C"))) {
				if (r1.atomExists("N")) {
					if (r2.atomExists("N")) {
						psi = CartesianGeometry::dihedral(r1.findAtom("N", true), r1.findAtom("CA", true), r1.findAtom("C", true), r2.findAtom("N", true));
					} else {
						if (r2.atomExists("NT")) {
							psi = CartesianGeometry::dihedral(r1.findAtom("N", true), r1.findAtom("CA", true), r1.findAtom("C", true), r2.findAtom("NT", true));
						}
					}
				} else if (r1.atomExists("NT")) {
					if (r2.atomExists("N")) {
						psi = CartesianGeometry::dihedral(r1.findAtom("NT", true), r1.findAtom("CA", true), r1.findAtom("C", true), r2.findAtom("N", true));
					} else if (r2.atomExists("NT")) {
						psi = CartesianGeometry::dihedral(r1.findAtom("NT", true), r1.findAtom("CA", true), r1.findAtom("C", true), r2.findAtom("NT", true));
					}
				}
			}
		}

		phi = (phi == -180.0 ? 180.0 : phi);
		psi = (psi == -180.0 ? 180.0 : psi);
		_dihed.push_back(vector<double>());
		_dihed[_dihed.size() - 1].push_back(phi);
		_dihed[_dihed.size() - 1].push_back(psi);
	}
}

void ProteinStruct::setDist()
{
	_dist.clear();

	int i, j;
	for (i = 0; i < (_pselca.size() - 1); i++)
	{
		for (j = i + 1; j < _pselca.size(); j++)
		{
			_dist.insert(make_pair(make_pair(i, j), _pselca[i]->distance(*_pselca[j])));
		}
	}

#if defined(DEBUG_DIST)
	for (i = 0 ; i < (_pselca.size() - 1); i++)
	{
		for (j = i + 1; j < _pselca.size(); j++)
		{
			printf("%d %d %f\n", i, j, _dist[make_pair(i, j)]);
		}
	}
	exit(-1);
#endif
}

void ProteinStruct::setPdbInfo() {
	_pdbinfo.clear();

	string resinfo;
	int i, j, count = 0;
	vector<Residue*> residues = _sys.getResidues();
	for (i = 0; i < residues.size(); i++) {
		resinfo.clear();

		Residue& r = *(residues[i]);
		for (j = 0; j < r.atomSize(); j++) {
			count++;
			resinfo += (r[j].pdbLine(r.getNum(), count) + string("\n"));
		}

		_pdbinfo.push_back(resinfo);
	}
}

void ProteinStruct::setProteinStruct()
{
	setBBCoords();
	setCACoords();	
	setDihed();
	setDist();
	setNumRes();
	setPdbInfo();
	setSeq();	

#if defined(DEBUG_PS)
	int i, j;
	cout << "CA atom coords:" << endl;
	for (i = 0; i < _cacoor.size(); i++)
	{
		printf("%f %f %f\n", _cacoor[i][0], _cacoor[i][1], _cacoor[i][2]);
	}
	cout << "dihedral angles:" << endl;
	for (i = 0; i < _dihed.size(); i++)
	{
		printf("%f %f\n", _dihed[i][PHI_IDX], _dihed[i][PSI_IDX]);
	}
	cout << "distance:" << endl;
	for (i = 0 ; i < (_numres - 1); i++)
	{
		for (j = i + 1; j < _numres; j++)
		{
			printf("%d %d %f\n", i, j, _dist[make_pair(i, j)]);
		}
	}
	cout << "number of residues: " << _numres << endl;
	cout << "PDB information:" << endl;
	for (i = 0; i < _pdbinfo.size(); i++)
	{
		cout << _pdbinfo[i];
	}
	cout << "sequence:" << endl;
	for (i = 0; i < _seq.size(); i++)
	{
		cout << _seq[i] << endl;
	}
	exit(-1);
#endif
}

void ProteinStruct::setSeq() {
	_seq.clear();

	int i;
	vector<Residue*> residues = _sys.getResidues();
	for (i = 0; i < residues.size(); i++) {
		_seq.push_back(residues[i]->getName());
	}
}

void ProteinStruct::writeBBCoords(fstream & ofs,const CreateOptions & copts)
{
	int i;
	if (!copts.getBinary())
	{
		writeString(ofs, "BBCOORDS\n", copts.getBinary());
	}
	for (i = 0; i < _bbcoor.size(); i++)
	{
		writeDatum(ofs, double(_bbcoor[i][0]), copts.getBinary());
		writeString(ofs, copts.getWordSep(), copts.getBinary());
		writeDatum(ofs, double(_bbcoor[i][1]), copts.getBinary());
		writeString(ofs, copts.getWordSep(), copts.getBinary());
		writeDatum(ofs, double(_bbcoor[i][2]), copts.getBinary());
		writeString(ofs, copts.getWordTer(), copts.getBinary());
	}
	if (!copts.getBinary())
	{
		writeString(ofs, "END\n", copts.getBinary());
	}
}

void ProteinStruct::writeCACoords(fstream & ofs,const CreateOptions & copts)
{
	int i;
	if (!copts.getBinary())
	{
		writeString(ofs, "CACOORDS\n", copts.getBinary());
	}
	for (i = 0; i < _cacoor.size(); i++)
	{
		writeDatum(ofs, double(_cacoor[i][0]), copts.getBinary());
		writeString(ofs, copts.getWordSep(), copts.getBinary());
		writeDatum(ofs, double(_cacoor[i][1]), copts.getBinary());
		writeString(ofs, copts.getWordSep(), copts.getBinary());
		writeDatum(ofs, double(_cacoor[i][2]), copts.getBinary());
		writeString(ofs, copts.getWordTer(), copts.getBinary());
	}
	if (!copts.getBinary())
	{
		writeString(ofs, "END\n", copts.getBinary());
	}
}

//	void ProteinStruct::writeDihedralAngles(fstream & ofs,const CreateOptions & copts)
//	{
//		int i;
//		if (!copts.bin)
//		{
//			writeString(ofs, "DIHEDRALANGLES\n", copts.bin);
//		}
//		for (i = 0; i < _dihed.size(); i++)
//		{
//			writeDatum(ofs, double(_dihed[i].first), copts.bin);
//			writeString(ofs, copts.sep, copts.bin);
//			writeDatum(ofs, double(_dihed[i].second), copts.bin);
//			writeString(ofs, copts.ter, copts.bin);
//		}
//		if (!copts.bin)
//		{
//			writeString(ofs, "END\n", copts.bin);
//		}
//	}
//	
//	void ProteinStruct::writeDistance(fstream & ofs,const CreateOptions & copts)
//	{
//		int i, j;
//		if (!copts.bin)
//		{
//			writeString(ofs, "DISTANCE\n", copts.bin);
//		}
//		for (i = 0; i < (_numres - 1); i++)
//		{
//			for (j = i + 1; j < _numres; j++)
//			{
//				writeDatum(ofs, int(i), copts.bin);
//				writeString(ofs, copts.sep, copts.bin);
//				writeDatum(ofs, int(j), copts.bin);
//				writeString(ofs, copts.sep, copts.bin);
//				writeDatum(ofs, double(_dist[make_pair(i, j)]), copts.bin);
//				writeString(ofs, copts.ter, copts.bin);
//			}
//		}
//		if (!copts.bin)
//		{
//			writeString(ofs, "END\n", copts.bin);
//		}
//	}

void ProteinStruct::writePdbInfo(fstream & ofs,const CreateOptions & copts)
{
	int i, countbytes;
	if (copts.getBinary())
	{
		countbytes = sizeof(int) * _pdbinfo.size();
		for (i = 0; i < _pdbinfo.size(); i++)
		{
			writeDatum(ofs, int(countbytes), copts.getBinary());
			countbytes += sizeof(int); // for number of chars
			countbytes += (sizeof(char) * _pdbinfo[i].size());
		}
		for (i = 0; i < _pdbinfo.size(); i++)
		{
			countbytes = (sizeof(char) * _pdbinfo[i].size());
			writeDatum(ofs, int(countbytes), copts.getBinary());
			writeString(ofs, _pdbinfo[i], copts.getBinary());
		}
	}
	else
	{
		writeString(ofs, "PDBINFO\n", copts.getBinary());
		for (i = 0; i < _pdbinfo.size(); i++)
		{
			writeString(ofs, _pdbinfo[i], copts.getBinary());
		}
		writeString(ofs, "END\n", copts.getBinary());
	}
}

void ProteinStruct::writeProteinStruct(fstream & ofs,const CreateOptions & copts)
{
	writeBBCoords(ofs, copts);
	writeCACoords(ofs, copts);
//		writeDihedralAngles(ofs, copts);
//		writeDistance(ofs, copts);
	writePdbInfo(ofs, copts);
	writeSeq(ofs, copts);
}

void ProteinStruct::writeProteinStructFile(fstream & ofs,const CreateOptions & copts)
{
	int countbytes = writeProteinStructFileHeader(ofs, copts, NUM_SEC_PS);
	writeProteinStruct(ofs, copts);

	if (copts.getBinary())
	{
		streampos begin, end;
		ofs.seekp(0, ios::beg);
		begin = ofs.tellp();
		ofs.seekp(0, ios::end);
		end = ofs.tellp();
		ASSERT(countbytes == (end - begin), "size of protein structure file not consistent");
	}

	cout << "Protein structure file created." << endl;
}

int ProteinStruct::writeProteinStructFileHeader(fstream & ofs, const CreateOptions & copts, const int & numsec)
{
	if (!copts.getBinary())
	{
		writeDatum(ofs, int(_numres), copts.getBinary());
		writeString(ofs, copts.getWordTer(), copts.getBinary());
		return -1;
	}

	int i, countbytes = sizeof(int) * numsec;

	// offset of BB coords section
	writeDatum(ofs, int(countbytes), copts.getBinary());
	countbytes += sizeof(double) * NUM_COORDS * _bbcoor.size();
	// offset of CA coords section
	writeDatum(ofs, int(countbytes), copts.getBinary());
	countbytes += sizeof(double) * NUM_COORDS * _cacoor.size();

//		// offset of dihed section
//		writeDatum(ofs, int(countbytes), copts.bin);
//		countbytes += sizeof(double) * 2 * _dihed.size();
//		// offset of dist section
//		writeDatum(ofs, int(countbytes), copts.bin);
//		countbytes += (sizeof(int) * 2 + sizeof(double)) * _dist.size();

	// number of residues
	writeDatum(ofs, int(_numres), copts.getBinary());	
	// offset of PDB info section
	writeDatum(ofs, int(countbytes), copts.getBinary());
	countbytes += sizeof(int) * _pdbinfo.size();
	for (i = 0; i < _pdbinfo.size(); i++)
	{
		countbytes += sizeof(int); // for number of chars
		countbytes += (sizeof(char) * _pdbinfo[i].size());
	}	
	// offset of seq section
	writeDatum(ofs, int(countbytes), copts.getBinary());
	countbytes += sizeof(char) * LEN_AA_CODE * _seq.size();

	return countbytes;
}

void ProteinStruct::writeSeq(fstream & ofs,const CreateOptions & copts)
{
	int i;
	if (!copts.getBinary())
	{
		writeString(ofs, "SEQ\n", copts.getBinary());
	}
	for (i = 0; i < _seq.size(); i++)
	{
		writeString(ofs, pad(_seq[i], LEN_AA_CODE), copts.getBinary());
		writeString(ofs, copts.getWordTer(), copts.getBinary());
	}
	if (!copts.getBinary())
	{
		writeString(ofs, "END\n", copts.getBinary());
	}
}


// function members for master only
ProteinStruct::ProteinStruct(const string & psfname)
{
	readProteinStruct(psfname);
}

double ProteinStruct::calcDistBB(const int & ri,const int & rj)
{
	if (ri == rj)
	{
		return 0.0;
	}

	int cai = ri * NUM_BBA + CA_IDX;
	int caj = rj * NUM_BBA + CA_IDX;
	return sqrt((_bbcoor[cai][0] - _bbcoor[caj][0]) * (_bbcoor[cai][0] - _bbcoor[caj][0])
		+ (_bbcoor[cai][1] - _bbcoor[caj][1]) * (_bbcoor[cai][1] - _bbcoor[caj][1])
		+ (_bbcoor[cai][2] - _bbcoor[caj][2]) * (_bbcoor[cai][2] - _bbcoor[caj][2]));
}

double ProteinStruct::calcDistCA(const int & ri, const int & rj)
{
	if (ri == rj)
	{
		return 0.0;
	}

	return sqrt((_cacoor[ri][0] - _cacoor[rj][0]) * (_cacoor[ri][0] - _cacoor[rj][0])
		+ (_cacoor[ri][1] - _cacoor[rj][1]) * (_cacoor[ri][1] - _cacoor[rj][1])
		+ (_cacoor[ri][2] - _cacoor[rj][2]) * (_cacoor[ri][2] - _cacoor[rj][2]));
} 

void ProteinStruct::closeProteinStructFile()
{
	if (_ifp != NULL)
	{
		fclose(_ifp);
	}
}

bool ProteinStruct::fullBB()
{
	for (int i = 0; i < _fullbb.size(); i++)
	{
		if (!_fullbb[i])
		{
			return false;
		}
	}
	return true;
}

void ProteinStruct::initPdbInfo()
{
	_pdbinfo.assign(_numres, string());
}

void ProteinStruct::readBBCoords()
{
	int i, j, offset, pos;
	const int coordlen = _numres * NUM_BBA * NUM_COORDS;
	double* xyz = (double*) malloc(coordlen * sizeof(double));
	_bbcoor.assign(_numres * NUM_BBA, vector<double>());
	_fullbb.assign(_numres, true);

	ASSERT(fseek(_ifp, sizeof(int) * SEC_BBCOOR, SEEK_SET) == 0, "could not seek offset of BB coords section in file %s", _fname.c_str());
	ASSERT(fread((void*) &offset, sizeof(int), 1, _ifp) == 1, "could not read offset of BB coords section in file %s", _fname.c_str());
	ASSERT(fseek(_ifp, offset, SEEK_SET) == 0, "could not seek BB coords section in file %s", _fname.c_str());
	ASSERT(fread((void*) xyz, sizeof(double), coordlen, _ifp) == coordlen, "could not read BB coords section in file %s", _fname.c_str());
	for (i = 0; i < _numres; i++)
	{
		for (j = 0; j < NUM_BBA; j++)
		{
			pos = i * NUM_BBA * NUM_COORDS + j * NUM_COORDS;
			_bbcoor[i * NUM_BBA + j].assign(xyz + pos, xyz + pos + NUM_COORDS);

			if (_fullbb[i] && (IMPOSSIBLE_COORD == _bbcoor[i * NUM_BBA + j][0]))
			{
				_fullbb[i] = false;
			}
		}
	}
	
	free((void*) xyz);

#if defined(DEBUG_BB)
	cout << "BB coords:\n";
	for (i = 0; i < _numres; i++)
	{
		printf("N\t%.3f\t%.3f\t%.3f\n", _bbcoor[4 * i + 0][0], _bbcoor[4 * i + 0][1], _bbcoor[4 * i + 0][2]);
		printf("CA\t%.3f\t%.3f\t%.3f\n", _bbcoor[4 * i + 1][0], _bbcoor[4 * i + 1][1], _bbcoor[4 * i + 1][2]);
		printf("C\t%.3f\t%.3f\t%.3f\n", _bbcoor[4 * i + 2][0], _bbcoor[4 * i + 2][1], _bbcoor[4 * i + 2][2]);
		printf("O\t%.3f\t%.3f\t%.3f\n", _bbcoor[4 * i + 3][0], _bbcoor[4 * i + 3][1], _bbcoor[4 * i + 3][2]);
		cout << _fullbb[i] << "\n";
		cout << "\n";
	}
	exit(-1);
#endif

}

void ProteinStruct::readCACoords()
{
	int i, offset, pos;
	const int coordlen = _numres * NUM_COORDS;
	double* xyz = (double*) malloc(coordlen * sizeof(double));
	_cacoor.assign(_numres, vector<double>());
	
	ASSERT(fseek(_ifp, sizeof(int) * SEC_CACOOR, SEEK_SET) == 0, "could not seek offset of CA coords section in file %s", _fname.c_str());
	ASSERT(fread((void*) &offset, sizeof(int), 1, _ifp) == 1, "could not read offset of CA coords section in file %s", _fname.c_str());
	ASSERT(fseek(_ifp, offset, SEEK_SET) == 0, "could not seek CA coords section in file %s", _fname.c_str());
	ASSERT(fread((void*) xyz, sizeof(double), coordlen, _ifp) == coordlen, "could not read CA coords section in file %s", _fname.c_str());
	for (i = 0; i < _numres; i++)
	{
		pos = i * NUM_COORDS;
		_cacoor[i].assign(xyz + pos, xyz + pos + NUM_COORDS);
	}

	free((void*) xyz);
}

void ProteinStruct::readNumRes()
{
	ASSERT(fseek(_ifp, sizeof(int) * SEC_NUMRES, SEEK_SET) == 0, "could not seek number of residues in file %s", _fname.c_str());	
	ASSERT(fread((void*) &_numres, sizeof(int), 1, _ifp) == 1, "could not read number of residues in file %s", _fname.c_str());
}

void ProteinStruct::readPdbInfo()
{
	_pdbinfo.clear();
	
	int i, offset, os_pdbinfo, pos, reslen, maxreslen = 10000; // guess as to what the max possible residue string will be
	char* resinfo = (char*) malloc((maxreslen + 1) * sizeof(char)); // 1 for '\0'
	
	ASSERT(fseek(_ifp, sizeof(int) * SEC_PDBINFO, SEEK_SET) == 0, "could not seek offset of PDB info section in file %s", _fname.c_str());
	ASSERT(fread((void*) &os_pdbinfo, sizeof(int), 1, _ifp) == 1, "could not read offset of PDB info section in file %s", _fname.c_str());
	ASSERT(fseek(_ifp, os_pdbinfo, SEEK_SET) == 0, "could not seek PDB info section in file %s", _fname.c_str());
	ASSERT(fread((void*) &offset, sizeof(int), 1, _ifp) == 1, "could not read offset of PDB info in file %s", _fname.c_str());
	pos = os_pdbinfo + offset;
	ASSERT(fseek(_ifp, pos, SEEK_SET) == 0, "could not seek PDB info in file %s", _fname.c_str());
	for (i = 0; i < _numres; i++)
	{
		ASSERT(fread((void*) &reslen, sizeof(int), 1, _ifp) == 1, "could not read length of PDB info in file %s", _fname.c_str());
		if (reslen > maxreslen)
		{
			maxreslen = reslen;
			free((void*) resinfo);
			resinfo = (char*) malloc((maxreslen + 1) * sizeof(char)); // 1 for '\0'
		}
		ASSERT(fread((void*) resinfo, sizeof(char), reslen, _ifp) == reslen, "could not read PDB info in file %s", _fname.c_str());
		resinfo[reslen] = '\0';
		_pdbinfo.push_back(string(resinfo));
	}
	free((void*) resinfo);
}

void ProteinStruct::readPdbInfo(const vector < int > & residx)
{
	int i;
	vector<int> ri;
	// check which residues need to actually be read
	for (i = 0; i < residx.size(); i++)
	{
		if (_pdbinfo[residx[i]].empty())
		{
			ri.push_back(residx[i]);
		}
	}
	if (ri.empty())
	{
		return;
	}

	vector<int> offset(ri.size());
	int os_pdbinfo, pos, reslen, maxreslen = 10000; // guess as to what the max possible residue string will be
	char* resinfo = (char*) malloc((maxreslen + 1) * sizeof(char)); // 1 for '\0'

	ASSERT(fseek(_ifp, sizeof(int) * SEC_PDBINFO, SEEK_SET) == 0, "could not seek offset of PDB info section in file %s", _fname.c_str());
	ASSERT(fread((void*) &os_pdbinfo, sizeof(int), 1, _ifp) == 1, "could not read offset of PDB info section in file %s", _fname.c_str());

	for (i = 0; i < ri.size(); i++)
	{
		pos = os_pdbinfo + sizeof(int) * ri[i];
		if (pos != ftell(_ifp))
		{
			ASSERT(fseek(_ifp, pos, SEEK_SET) == 0, "could not seek offset of PDB info in file %s", _fname.c_str());
		}
		ASSERT(fread((void*) &offset[i], sizeof(int), 1, _ifp) == 1, "could not read offset of PDB info in file %s", _fname.c_str());
	}
	for (i = 0; i < ri.size(); i++)
	{
		pos = os_pdbinfo + offset[i];
		if (pos != ftell(_ifp))
		{
			ASSERT(fseek(_ifp, pos, SEEK_SET) == 0, "could not seek PDB info in file %s", _fname.c_str());
		}
		ASSERT(fread((void*) &reslen, sizeof(int), 1, _ifp) == 1, "could not read length of PDB info in file %s", _fname.c_str());
		if (reslen > maxreslen)
		{
			maxreslen = reslen;
			free((void*) resinfo);
			resinfo = (char*) malloc((maxreslen + 1) * sizeof(char)); // 1 for '\0'
		}
		ASSERT(fread((void*) resinfo, sizeof(char), reslen, _ifp) == reslen, "could not read PDB info in file %s", _fname.c_str());
		resinfo[reslen] = '\0';
		_pdbinfo[ri[i]] = string(resinfo);
	}
	free((void*) resinfo);
}

void ProteinStruct::readProteinStruct(const string & psfname)
{
	_fname = psfname;

	openFileC(_ifp, _fname, "rb");

	readNumRes();

#if defined(DEBUG_PS)
	cout << "file name: " << _fname << "\n";
	cout << "number of residues: " << _numres << "\n";
	exit(-1);
#endif
}

void ProteinStruct::readSeq(char * & seq)
{
	delete [] seq;

	int offset;
	const int seqlen = LEN_AA_CODE * _numres;
	seq = new char [seqlen];

	ASSERT(fseek(_ifp, sizeof(int) * SEC_SEQ, SEEK_SET) == 0, "could not seek offset of sequence section in file %s", _fname.c_str());
	ASSERT(fread((void*) &offset, sizeof(int), 1, _ifp) == 1, "could not read offset of sequence section in file %s", _fname.c_str());
	ASSERT(fseek(_ifp, offset, SEEK_SET) == 0, "could not seek sequence section in file %s", _fname.c_str());
	ASSERT(fread((void*) seq, sizeof(char), seqlen, _ifp) == seqlen, "could not read sequence section in file %s", _fname.c_str());
}

