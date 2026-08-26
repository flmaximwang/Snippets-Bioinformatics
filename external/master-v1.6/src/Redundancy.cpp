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

#include "Redundancy.h"

bool alignChains(const vector<string> & acceptAA, const vector<string> & seqA, const vector<string> & seqB, const double & seqIdenThresh, AtomPointerVector & caA, AtomPointerVector & caB, const int & mode, const double & optRmsdThresh, vector<double> & t, vector<vector<double> > & u, const double & contRmsdThresh, const vector<double> & ct, const vector<vector<double> > & cu)
{
	// True: two chains are redundant
	// False: two chains are non-redundant

	// mode=1: align central chains, and calculate optimal RMSD as well as rotation and translation matrices
	// mode=0: align corresponding neighbors, and calculate optimal RMSD as well as contact RMSD

	vector<string> seqAAligned, seqBAligned;
	if (alignSeqNeedlemanWunsch(acceptAA ,seqA, seqB, seqIdenThresh, seqAAligned, seqBAligned))
	{
		AtomPointerVector caAAligned, caBAligned;
		int i, j = 0, k = 0;
		for (i = 0; i < seqAAligned.size(); i++)
		{
			if ((seqAAligned[i] != "---") && (seqBAligned[i] != "---"))
			{
				caAAligned.push_back(caA[j]);
				caBAligned.push_back(caB[k]);
			}
			if (seqAAligned[i] != "---")
			{
				j++;
			}
			if (seqBAligned[i] != "---")
			{
				k++;
			}
		}
		ASSERT((j == caA.size()) && (k == caB.size()), "number of Ca atoms not consistent");
		if (calcRmsdKabsch(caAAligned, caBAligned, mode, t, u) <= optRmsdThresh)
		{
			if (mode)
			{
				return true;
			}
			double contRmsd = 0.0, xa, ya, za, xb, yb, zb;
			for (i = 0; i < caAAligned.size(); i++)
			{
				xa = ct[0] + cu[0][0] * caAAligned[i]->getX() + cu[0][1] * caAAligned[i]->getY() + cu[0][2] * caAAligned[i]->getZ();
				ya = ct[1] + cu[1][0] * caAAligned[i]->getX() + cu[1][1] * caAAligned[i]->getY() + cu[1][2] * caAAligned[i]->getZ();
				za = ct[2] + cu[2][0] * caAAligned[i]->getX() + cu[2][1] * caAAligned[i]->getY() + cu[2][2] * caAAligned[i]->getZ();
				xb = caBAligned[i]->getX();
				yb = caBAligned[i]->getY();
				zb = caBAligned[i]->getZ();
				contRmsd += ((xa - xb) * (xa - xb) + (ya - yb) * (ya - yb) + (za - zb) * (za - zb));
			}
			contRmsd = sqrt(contRmsd / double(caAAligned.size()));
			if (contRmsd <= contRmsdThresh)
			{
				return true;
			}
		}
	}

	return false;
}

bool alignSeqNeedlemanWunsch(const vector<string> & acceptAA, const vector<string> & seqA, const vector<string> & seqB, const double & seqIdenThresh, vector<string> & seqAAligned, vector<string> & seqBAligned)
{
	// True: two sequences are redundant
	// False: two sequences are not redundant
	if ((0 == seqA.size()) && (0 == seqB.size()))
	{
		return false;
	}

	seqAAligned.clear();
	seqBAligned.clear();

	map<string, string> aa2aa;
	aa2aa["ALA"] = "A";
	aa2aa["ARG"] = "R";
	aa2aa["ASN"] = "N";
	aa2aa["ASP"] = "D";
	aa2aa["CYS"] = "C";
	aa2aa["GLN"] = "Q";
	aa2aa["GLU"] = "E";
	aa2aa["GLY"] = "G";
	aa2aa["HIS"] = "H"; aa2aa["HSC"] = "H"; aa2aa["HSD"] = "H"; aa2aa["HSE"] = "H"; aa2aa["HSP"] = "H";
	aa2aa["ILE"] = "I";
	aa2aa["LEU"] = "L";
	aa2aa["LYS"] = "K";
	aa2aa["MET"] = "M"; aa2aa["MSE"] = "M";
	aa2aa["PHE"] = "F";
	aa2aa["PRO"] = "P";
	aa2aa["SER"] = "S";
	aa2aa["THR"] = "T";
	aa2aa["TRP"] = "W";
	aa2aa["TYR"] = "Y";
	aa2aa["VAL"] = "V";

	aa2aa["CSO"] = "CSO";
	aa2aa["HIP"] = "HIP";
	aa2aa["PTR"] = "PTR";
	aa2aa["SEC"] = "SEC";
	aa2aa["SEP"] = "SEP";
	aa2aa["TPO"] = "TPO";

	for (unsigned int i = 0; i < acceptAA.size(); i++)
	{
		if (aa2aa.find(acceptAA[i]) == aa2aa.end())
		{
			aa2aa[acceptAA[i]] = acceptAA[i];
		}
	}

	map<string, int> aa2idx;
	aa2idx["ALA"] = 0;
	aa2idx["ARG"] = 1;
	aa2idx["ASN"] = 2;
	aa2idx["ASP"] = 3;
	aa2idx["CYS"] = 4;
	aa2idx["GLN"] = 5;
	aa2idx["GLU"] = 6;
	aa2idx["GLY"] = 7;
	aa2idx["HIS"] = 8; aa2idx["HSC"] = 8; aa2idx["HSD"] = 8; aa2idx["HSE"] = 8; aa2idx["HSP"] = 8;
	aa2idx["ILE"] = 9;
	aa2idx["LEU"] = 10;
	aa2idx["LYS"] = 11;
	aa2idx["MET"] = 12; aa2idx["MSE"] = 12;
	aa2idx["PHE"] = 13;
	aa2idx["PRO"] = 14;
	aa2idx["SER"] = 15;
	aa2idx["THR"] = 16;
	aa2idx["TRP"] = 17;
	aa2idx["TYR"] = 18;
	aa2idx["VAL"] = 19;

	aa2idx["CSO"] = 22; aa2idx["HIP"] = 22; aa2idx["PTR"] = 22; aa2idx["SEC"] = 22; aa2idx["SEP"] = 22; aa2idx["TPO"] = 22;

	for (unsigned int i = 0; i < acceptAA.size(); i++)
	{
		if (aa2idx.find(acceptAA[i]) == aa2idx.end())
		{
			aa2idx[acceptAA[i]] = 22;
		}
	}

	const int matSub[24][24] = // PAM250
	{
		{ 2, -2,  0,  0, -2,  0,  0,  1, -1, -1, -2, -1, -1, -3,  1,  1,  1, -6, -3,  0,  0,  0,  0, -8}, // A: ALA
		{-2,  6,  0, -1, -4,  1, -1, -3,  2, -2, -3,  3,  0, -4,  0,  0, -1,  2, -4, -2, -1,  0, -1, -8}, // R: ARG
		{ 0,  0,  2,  2, -4,  1,  1,  0,  2, -2, -3,  1, -2, -3,  0,  1,  0, -4, -2, -2,  2,  1,  0, -8}, // N: ASN
		{ 0, -1,  2,  4, -5,  2,  3,  1,  1, -2, -4,  0, -3, -6, -1,  0,  0, -7, -4, -2,  3,  3, -1, -8}, // D: ASP
		{-2, -4, -4, -5, 12, -5, -5, -3, -3, -2, -6, -5, -5, -4, -3,  0, -2, -8,  0, -2, -4, -5, -3, -8}, // C: CYS
		{ 0,  1,  1,  2, -5,  4,  2, -1,  3, -2, -2,  1, -1, -5,  0, -1, -1, -5, -4, -2,  1,  3, -1, -8}, // Q: GLN
		{ 0, -1,  1,  3, -5,  2,  4,  0,  1, -2, -3,  0, -2, -5, -1,  0,  0, -7, -4, -2,  3,  3, -1, -8}, // E: GLU
		{ 1, -3,  0,  1, -3, -1,  0,  5, -2, -3, -4, -2, -3, -5,  0,  1,  0, -7, -5, -1,  0,  0, -1, -8}, // G: GLY
		{-1,  2,  2,  1, -3,  3,  1, -2,  6, -2, -2,  0, -2, -2,  0, -1, -1, -3,  0, -2,  1,  2, -1, -8}, // H: HIS, HSC, HSD, HSE, HSP
		{-1, -2, -2, -2, -2, -2, -2, -3, -2,  5,  2, -2,  2,  1, -2, -1,  0, -5, -1,  4, -2, -2, -1, -8}, // I: ILE
		{-2, -3, -3, -4, -6, -2, -3, -4, -2,  2,  6, -3,  4,  2, -3, -3, -2, -2, -1,  2, -3, -3, -1, -8}, // L: LEU
		{-1,  3,  1,  0, -5,  1,  0, -2,  0, -2, -3,  5,  0, -5, -1,  0,  0, -3, -4, -2,  1,  0, -1, -8}, // K: LYS
		{-1,  0, -2, -3, -5, -1, -2, -3, -2,  2,  4,  0,  6,  0, -2, -2, -1, -4, -2,  2, -2, -2, -1, -8}, // M: MET, MSE
		{-3, -4, -3, -6, -4, -5, -5, -5, -2,  1,  2, -5,  0,  9, -5, -3, -3,  0,  7, -1, -4, -5, -2, -8}, // F: PHE
		{ 1,  0,  0, -1, -3,  0, -1,  0,  0, -2, -3, -1, -2, -5,  6,  1,  0, -6, -5, -1, -1,  0, -1, -8}, // P: PRO
		{ 1,  0,  1,  0,  0, -1,  0,  1, -1, -1, -3,  0, -2, -3,  1,  2,  1, -2, -3, -1,  0,  0,  0, -8}, // S: SER
		{ 1, -1,  0,  0, -2, -1,  0,  0, -1,  0, -2,  0, -1, -3,  0,  1,  3, -5, -3,  0,  0, -1,  0, -8}, // T: THR
		{-6,  2, -4, -7, -8, -5, -7, -7, -3, -5, -2, -3, -4,  0, -6, -2, -5, 17,  0, -6, -5, -6, -4, -8}, // W: TRP
		{-3, -4, -2, -4,  0, -4, -4, -5,  0, -1, -1, -4, -2,  7, -5, -3, -3,  0, 10, -2, -3, -4, -2, -8}, // Y: TYR
		{ 0, -2, -2, -2, -2, -2, -2, -1, -2,  4,  2, -2,  2, -1, -1, -1,  0, -6, -2,  4, -2, -2, -1, -8}, // V: VAL
		{ 0, -1,  2,  3, -4,  1,  3,  0,  1, -2, -3,  1, -2, -4, -1,  0,  0, -5, -3, -2,  3,  2, -1, -8}, // B: ASP & ASN
		{ 0,  0,  1,  3, -5,  3,  3,  0,  2, -2, -3,  0, -2, -5,  0,  0, -1, -6, -4, -2,  2,  3, -1, -8}, // Z: GLU & GLN
		{ 0, -1,  0, -1, -3, -1, -1, -1, -1, -1, -1, -1, -1, -2, -1,  0,  0, -4, -2, -1, -1, -1, -1, -8}, // X: unnatural amino acids
		{-8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8,  1}  // gap
	};

	vector<vector<int> > M(seqA.size() + 1, vector<int>(seqB.size() + 1)); // i aligned with j
	vector<vector<int> > I(seqA.size() + 1, vector<int>(seqB.size() + 1)); // i aligned with gap
	vector<vector<int> > J(seqA.size() + 1, vector<int>(seqB.size() + 1)); // j aligned with gap
	const unsigned int diag = 1, up = 2, left = 3;
	unsigned int d, i, j, nIden;
	const int inf = 999999, openGap = -11, extGap = -1;
	int score, scoreCur;
	const string seqGap = "---";
	double seqIden;

	// initialize dynamic programming matrices
	M[0][0] = 0;
	I[0][0] = 0;
	J[0][0] = 0;
	for (j = 1; j <= seqB.size(); j++) // first row
	{
		M[0][j] = -inf;
		I[0][j] = -inf;
		J[0][j] = openGap + (j * extGap);
	}
	for (i = 1; i <= seqA.size(); i++) // first column
	{
		M[i][0] = -inf;
		I[i][0] = openGap + (i * extGap);
		J[i][0] = -inf;
	}
	// complete dynamic programming matrix
	for (i = 1; i <= seqA.size(); i++)
	{
		for (j = 1; j <= seqB.size(); j++)
		{
			M[i][j] = max(M[i - 1][j - 1], max(I[i - 1][j - 1], J[i - 1][j - 1])) + matSub[aa2idx[seqA[i - 1]]][aa2idx[seqB[j - 1]]];
			I[i][j] = max(M[i - 1][j] + openGap + extGap, I[i - 1][j] + extGap);
			J[i][j] = max(M[i][j - 1] + openGap + extGap, J[i][j - 1] + extGap);
		}
	}

	score = max(M[seqA.size()][seqB.size()], max(I[seqA.size()][seqB.size()], J[seqA.size()][seqB.size()]));
	// trace back
	if (M[seqA.size()][seqB.size()] == score)
	{
		d = diag;
	}
	else
	{
		if (I[seqA.size()][seqB.size()] == score)
		{
			d = up;
		}
		else
		{
			if (J[seqA.size()][seqB.size()] == score)
			{
				d = left;
			}
			else
			{
				error("unidentified alignment score %d", score);
			}
		}
	}
	scoreCur = score;
	i = seqA.size();
	j = seqB.size();

	do {
		if (diag == d)
		{
			if ((M[i - 1][j - 1] + matSub[aa2idx[seqA[i - 1]]][aa2idx[seqB[j - 1]]]) == scoreCur)
			{
				scoreCur = M[i - 1][j - 1];
				d = diag;
			}
			else
			{
				if ((I[i - 1][j - 1] + matSub[aa2idx[seqA[i - 1]]][aa2idx[seqB[j - 1]]]) == scoreCur)
				{
					scoreCur = I[i - 1][j - 1];
					d = up;
				}
				else
				{
					if ((J[i - 1][j - 1] + matSub[aa2idx[seqA[i - 1]]][aa2idx[seqB[j - 1]]]) == scoreCur)
					{
						scoreCur = J[i - 1][j - 1];
						d = left;
					}
					else
					{
						error("unidentified alignment score %d", scoreCur);
					}
				}
			}
			
			seqAAligned.push_back(seqA[i - 1]);
			seqBAligned.push_back(seqB[j - 1]);
			i--;
			j--;
		}
		else
		{
			if (up == d)
			{
				if ((M[i - 1][j] + openGap + extGap) == scoreCur)
				{
					scoreCur = M[i - 1][j];
					d = diag;
				}
				else
				{
					if ((I[i - 1][j] + extGap) == scoreCur)
					{
						scoreCur = I[i - 1][j];
						d = up;
					}
					else
					{
						error("unidentified alignment score %d", scoreCur);
					}
				}
				
				seqAAligned.push_back(seqA[i - 1]);
				seqBAligned.push_back(seqGap);
				i--;
			}
			else
			{
				if (left == d)
				{
					if ((M[i][j - 1] + openGap + extGap) == scoreCur)
					{
						scoreCur = M[i][j - 1];
						d = diag;
					}
					else
					{
						if ((J[i][j - 1] + extGap) == scoreCur)
						{
							scoreCur = J[i][j - 1];
							d = left;
						}
						else
						{
							error("unidentified alignment score %d", scoreCur);
						}
					}
					
					seqAAligned.push_back(seqGap);
					seqBAligned.push_back(seqB[j - 1]);
					j--;
				}
				else
				{
					error("unidentified direction %d", d);
				}
			}
		}
	} while ((i > 0) && (j > 0));
	ASSERT((0 == i) || (0 == j), "backtrace not complete");
	while (i > 0) // first column
	{
		seqAAligned.push_back(seqA[i - 1]);
		seqBAligned.push_back(seqGap);
		i--;
	}
	while (j > 0) // first row
	{
		seqBAligned.push_back(seqB[j - 1]);
		seqAAligned.push_back(seqGap);
		j--;
	}
	reverse(seqAAligned.begin(), seqAAligned.end());
	reverse(seqBAligned.begin(), seqBAligned.end());
	ASSERT(seqAAligned.size() == seqBAligned.size(), "length of aligned sequences not consistent");
	
#if defined(DEBUG_SEQALI)
	cout << "sequence A:" << endl;
	for (i = 0; i < seqA.size(); i++)
	{
		cout << " " << seqA[i];
	}
	cout << endl;
	cout << "sequence B:" << endl;
	for (j = 0; j < seqB.size(); j++)
	{
		cout << " " << seqB[j];
	}
	cout << endl;
	cout << "aligned sequence A:" << endl;
	for (i = 0; i < seqAAligned.size(); i++)
	{
		cout << " " << seqAAligned[i];
	}
	cout << endl;
	cout << "aligned sequence B:" << endl;
	for (j = 0; j < seqBAligned.size(); j++)
	{
		cout << " " << seqBAligned[j];
	}
	cout << endl;
#endif

	// calculate sequence identity
	nIden = 0;
	for (i = 0; i < seqAAligned.size(); i++)
	{
		if ((seqGap != seqAAligned[i]) && (seqGap != seqBAligned[i]) && (aa2aa[seqAAligned[i]] == aa2aa[seqBAligned[i]]))
		{
			nIden++;
		}
	}
	seqIden = double(nIden) / double(seqAAligned.size());

#if defined(DEBUG_SEQALI)
	cout << "sequence identity: " << seqIden << endl;
	exit(-1);
#endif
	
	if (seqIden < seqIdenThresh)
	{
		return false;
	}
	return true;
}

// defines whether two crystal units are in contact
// count how many pairs of atoms from two units are within a distance cutoff
bool areUnitsInContact(AtomPointerVector& A, AtomPointerVector& B) {
	double dcut = 6.0;
  	double dcut2 = dcut * dcut;
  	int i, n = 0, x, y, z, xi, yi, zi, ai;
  	int nthresh = min(20, int(ceil(0.5 * (A.size() + B.size()) / 10)));
	string atomname;

  	// create a grid and place B into it
  	map<int, map<int, map<int, AtomPointerVector> > > grid;
  	CartesianPoint c = A.getGeometricCenter();
  	for (i = 0; i < B.size(); i++)
  	{
		atomname = B[i]->getName();
		if (atomname[0] == 'H')
		{
			continue;
		}
		
		gridPoint(B[i], c, dcut, &x, &y, &z);
    	if (grid.find(x) == grid.end())
    	{
			grid[x] = map<int, map<int, AtomPointerVector> >();
    	}
    	if (grid[x].find(y) == grid[x].end())
    	{
			grid[x][y] = map<int, AtomPointerVector> ();
    	}
    	if (grid[x][y].find(z) == grid[x][y].end())
    	{
			grid[x][y][z] = AtomPointerVector ();
    	}
    	grid[x][y][z].push_back(B[i]);
  	}

  	// visit all atoms in A
  	for (i = 0; i < A.size(); i++)
  	{
		atomname = A[i]->getName();
		if (atomname[0] == 'H')
		{
			continue;
		}
		
    	// what grid box does this atom fit into?
    	gridPoint(A[i], c, dcut, &x, &y, &z);
    	for (xi = x - 1; xi <= x + 1; xi++)
    	{
      		if (grid.find(xi) == grid.end())
      		{
	  			continue;
      		}
      		map<int, map<int, AtomPointerVector> >& gridX = grid[xi];
      		for (yi = y - 1; yi <= y + 1; yi++)
      		{
        		if (gridX.find(yi) == gridX.end())
        		{
					continue;
        		}
        		map<int, AtomPointerVector>& gridXY = gridX[yi];
        		for (zi = z - 1; zi <= z + 1; zi++)
        		{
          			if (gridXY.find(zi) == gridXY.end())
          			{
		  				continue;
          			}
          			AtomPointerVector& gridXYZ = gridXY[zi];
          			for (ai = 0; ai < gridXYZ.size(); ai++)
          			{
            			if (A[i]->distance2(*(gridXYZ[ai])) < dcut2)
            			{
							n++;
            			}
          			}
        		}
      		}
    	}
    	if (n > nthresh)
    	{
			return true;
    	}
  	}
  	return false;
}

int calcNumOvlpChains(Chain* curChain, vector<Chain*> & curChainSet, map<Chain*, vector<Chain*> > & chainContacts)
{
	int i, j, numOvlpChains = 0;
	map<Chain*, bool> chainInc;

	for (i = 0; i < curChainSet.size(); i++)
	{
		chainInc[curChainSet[i]] = true;
		for (j = 0; j < chainContacts[curChainSet[i]].size(); j++)
		{
			chainInc[chainContacts[curChainSet[i]][j]] = true;
		}
	}

	if (chainInc.find(curChain) != chainInc.end())
	{
		numOvlpChains++;
	}
	for (i = 0; i < chainContacts[curChain].size(); i++)
	{
		if (chainInc.find(chainContacts[curChain][i]) != chainInc.end())
		{
			numOvlpChains++;
		}
	}

	return numOvlpChains;
}

/**************************************************************************
  Implemetation of Kabsch algoritm for finding the best rotation matrix
---------------------------------------------------------------------------
  x    - x(i,m) are coordinates of atom m in set x            (input)
  y    - y(i,m) are coordinates of atom m in set y            (input)
  n    - n is number of atom pairs                            (input)
  mode  - 0:calculate rmsd only                               (input)
          1:calculate rmsd,u,t                                (takes longer)
  rms   - sum of w*(ux+t-y)**2 over all atom pairs            (output)
  u    - u(i,j) is   rotation  matrix for best superposition  (output)
  t    - t(i)   is translation vector for best superposition  (output)
**************************************************************************/
double calcRmsdKabsch(AtomPointerVector &_align, AtomPointerVector &_ref, int mode, vector<double> & t, vector<vector<double> > & u)
{
	int i, j, m, m1, l, k;
	double e0, rms1, d, h, g;
	double cth, sth, sqrth, p, det, sigma;  
	double xc[3], yc[3];
	double a[3][3], b[3][3], r[3][3], e[3], rr[6], ss[6];
	double sqrt3=1.73205080756888, tol=0.01;
	int ip[]={0, 1, 3, 1, 2, 4, 3, 4, 5};
	int ip2312[]={1, 2, 0, 1};
	
	int a_failed=0, b_failed=0;
	double epsilon=0.00000001;
	
	int n=_ref.size();
	ASSERT((n >= 0) && (n == _align.size()), "invalid protein length for calculating RMSD");

	if (0 == n)
	{
		return 0.0;
	}
	
	//initializtation
	double rmsd=0;
	rms1=0;
	e0=0;
	for (i=0; i<3; i++) {
		xc[i]=0.0;
		yc[i]=0.0;
		t[i]=0.0;
		for (j=0; j<3; j++) {
			u[i][j]=0.0;
			r[i][j]=0.0;
			a[i][j]=0.0;
			if (i==j) {
				u[i][j]=1.0;
				a[i][j]=1.0;
			}
		}
	} 

	//compute centers for vector sets x, y
	for(i=0; i<n; i++){
		xc[0] += _align[i]->getX();
		xc[1] += _align[i]->getY();
		xc[2] += _align[i]->getZ();
		
		yc[0] += _ref[i]->getX();
		yc[1] += _ref[i]->getY();
		yc[2] += _ref[i]->getZ();
	}
	for(i=0; i<3; i++){
		xc[i] = xc[i]/(double)n;
		yc[i] = yc[i]/(double)n;        
	}
	
	//compute e0 and matrix r
	for (m=0; m<n; m++) {
		e0 += (_align[m]->getX()-xc[0])*(_align[m]->getX()-xc[0]) \
		  +(_ref[m]->getX()-yc[0])*(_ref[m]->getX()-yc[0]);
		e0 += (_align[m]->getY()-xc[1])*(_align[m]->getY()-xc[1]) \
		  +(_ref[m]->getY()-yc[1])*(_ref[m]->getY()-yc[1]);
		e0 += (_align[m]->getZ()-xc[2])*(_align[m]->getZ()-xc[2]) \
		  +(_ref[m]->getZ()-yc[2])*(_ref[m]->getZ()-yc[2]);
		r[0][0] += (_ref[m]->getX() - yc[0])*(_align[m]->getX() - xc[0]);
		r[0][1] += (_ref[m]->getX() - yc[0])*(_align[m]->getY() - xc[1]);
		r[0][2] += (_ref[m]->getX() - yc[0])*(_align[m]->getZ() - xc[2]);
		r[1][0] += (_ref[m]->getY() - yc[1])*(_align[m]->getX() - xc[0]);
		r[1][1] += (_ref[m]->getY() - yc[1])*(_align[m]->getY() - xc[1]);
		r[1][2] += (_ref[m]->getY() - yc[1])*(_align[m]->getZ() - xc[2]);
		r[2][0] += (_ref[m]->getZ() - yc[2])*(_align[m]->getX() - xc[0]);
		r[2][1] += (_ref[m]->getZ() - yc[2])*(_align[m]->getY() - xc[1]);
		r[2][2] += (_ref[m]->getZ() - yc[2])*(_align[m]->getZ() - xc[2]);
	}
	//compute determinat of matrix r
	det = r[0][0] * ( r[1][1]*r[2][2] - r[1][2]*r[2][1] )		\
	- r[0][1] * ( r[1][0]*r[2][2] - r[1][2]*r[2][0] )		\
	+ r[0][2] * ( r[1][0]*r[2][1] - r[1][1]*r[2][0] ); 
	sigma=det;
	
	//compute tras(r)*r
	m = 0;
	for (j=0; j<3; j++) {
		for (i=0; i<=j; i++) {            
			rr[m]=r[0][i]*r[0][j]+r[1][i]*r[1][j]+r[2][i]*r[2][j];
			m++;
		}
	}
	
	double spur=(rr[0]+rr[2]+rr[5]) / 3.0;
	double cof = (((((rr[2]*rr[5] - rr[4]*rr[4]) + rr[0]*rr[5])	\
		  - rr[3]*rr[3]) + rr[0]*rr[2]) - rr[1]*rr[1]) / 3.0;
	det = det*det; 
	
	for (i=0; i<3; i++){
		e[i]=spur;
	}
	
	if (spur>0) {
		d = spur*spur;
		h = d - cof;
		g = (spur*cof - det)/2.0 - spur*h;

		if (h>0) {
			sqrth = sqrt(h);
			d = h*h*h - g*g;
			if(d<0.0) d=0.0;
			d = atan2( sqrt(d), -g ) / 3.0;			
			cth = sqrth * cos(d);
			sth = sqrth*sqrt3*sin(d);
			e[0]= (spur + cth) + cth;
			e[1]= (spur - cth) + sth;            
			e[2]= (spur - cth) - sth;

			if (mode!=0) {//compute a                
				for (l=0; l<3; l=l+2) {
					d = e[l];  
					ss[0] = (d-rr[2]) * (d-rr[5])  - rr[4]*rr[4];
					ss[1] = (d-rr[5]) * rr[1]      + rr[3]*rr[4];
					ss[2] = (d-rr[0]) * (d-rr[5])  - rr[3]*rr[3];
					ss[3] = (d-rr[2]) * rr[3]      + rr[1]*rr[4];
					ss[4] = (d-rr[0]) * rr[4]      + rr[1]*rr[3];                
					ss[5] = (d-rr[0]) * (d-rr[2])  - rr[1]*rr[1]; 

					if (fabs(ss[0])<=epsilon) ss[0]=0.0;
					if (fabs(ss[1])<=epsilon) ss[1]=0.0;
					if (fabs(ss[2])<=epsilon) ss[2]=0.0;
					if (fabs(ss[3])<=epsilon) ss[3]=0.0;
					if (fabs(ss[4])<=epsilon) ss[4]=0.0;
					if (fabs(ss[5])<=epsilon) ss[5]=0.0;

					if (fabs(ss[0]) >= fabs(ss[2])) {
						j=0;                    
						if( fabs(ss[0]) < fabs(ss[5])){
							j = 2;
						}
					} else if ( fabs(ss[2]) >= fabs(ss[5]) ){
						j = 1;
					} else {
						j = 2;
					}

					d = 0.0;
					j = 3 * j;
					for (i=0; i<3; i++) {
						k=ip[i+j];
						a[i][l] = ss[k];
						d = d + ss[k]*ss[k];						
					} 


					//if( d > 0.0 ) d = 1.0 / sqrt(d);
					if (d > epsilon) d = 1.0 / sqrt(d);
					else d=0.0;
					for (i=0; i<3; i++) {
						a[i][l] = a[i][l] * d;
					}               
				}//for l

				d = a[0][0]*a[0][2] + a[1][0]*a[1][2] + a[2][0]*a[2][2];
				if ((e[0] - e[1]) > (e[1] - e[2])) {
					m1=2;
					m=0;
				} else {
					m1=0;
					m=2;                
				}
				p=0;
				for(i=0; i<3; i++){
					a[i][m1] = a[i][m1] - d*a[i][m];
					p = p + a[i][m1]*a[i][m1];
				}
				if (p <= tol) {
					p = 1.0;
					for (i=0; i<3; i++) {
						if (p < fabs(a[i][m])){
							continue;
						}
						p = fabs( a[i][m] );
						j = i;                    
					}
					k = ip2312[j];
					l = ip2312[j+1];
					p = sqrt( a[k][m]*a[k][m] + a[l][m]*a[l][m] ); 
					if (p > tol) {
						a[j][m1] = 0.0;
						a[k][m1] = -a[l][m]/p;
						a[l][m1] =  a[k][m]/p;                                                       
					} else {//goto 40
						a_failed=1;
					}     
				} else {//if p<=tol
					p = 1.0 / sqrt(p);
					for(i=0; i<3; i++){
						a[i][m1] = a[i][m1]*p;
					}                                  
				}//else p<=tol  
				if (a_failed!=1) {
					a[0][1] = a[1][2]*a[2][0] - a[1][0]*a[2][2];
					a[1][1] = a[2][2]*a[0][0] - a[2][0]*a[0][2];
					a[2][1] = a[0][2]*a[1][0] - a[0][0]*a[1][2];       
				}                                   
			}//if(mode!=0)       
		}//h>0

		//compute b anyway
		if (mode!=0 && a_failed!=1) {//a is computed correctly
			//compute b
			for (l=0; l<2; l++) {
				d=0.0;
				for(i=0; i<3; i++){
					b[i][l] = r[i][0]*a[0][l] + r[i][1]*a[1][l] + r[i][2]*a[2][l];
					d = d + b[i][l]*b[i][l];
				}
				//if( d > 0 ) d = 1.0 / sqrt(d);
				if (d > epsilon) d = 1.0 / sqrt(d);
				else d=0.0;
				for (i=0; i<3; i++) {
					b[i][l] = b[i][l]*d;
				}                
			}            
			d = b[0][0]*b[0][1] + b[1][0]*b[1][1] + b[2][0]*b[2][1];
			p=0.0;

			for (i=0; i<3; i++) {
				b[i][1] = b[i][1] - d*b[i][0];
				p += b[i][1]*b[i][1];
			}

			if (p <= tol) {
				p = 1.0;
				for (i=0; i<3; i++) {
					if (p<fabs(b[i][0])) {
						continue;
					}
					p = fabs( b[i][0] );
					j=i;
				}
				k = ip2312[j];
				l = ip2312[j+1];
				p = sqrt( b[k][0]*b[k][0] + b[l][0]*b[l][0] ); 
				if (p > tol) {
					b[j][1] = 0.0;
					b[k][1] = -b[l][0]/p;
					b[l][1] =  b[k][0]/p;        
				} else {
					//goto 40
					b_failed=1;
				}                
			} else {//if( p <= tol )
				p = 1.0 / sqrt(p);
				for(i=0; i<3; i++){
					b[i][1]=b[i][1]*p;
				}
			}            
			if (b_failed!=1){
				b[0][2] = b[1][0]*b[2][1] - b[1][1]*b[2][0];
				b[1][2] = b[2][0]*b[0][1] - b[2][1]*b[0][0];
				b[2][2] = b[0][0]*b[1][1] - b[0][1]*b[1][0]; 
				//compute u
				for (i=0; i<3; i++){
					for(j=0; j<3; j++){
						u[i][j] = b[i][0]*a[j][0] + b[i][1]*a[j][1]	\
								+ b[i][2]*a[j][2];
					}
				}
			}

			//compute t
			for(i=0; i<3; i++){
				t[i] = ((yc[i] - u[i][0]*xc[0]) - u[i][1]*xc[1])	\
						- u[i][2]*xc[2];
			}            
		}//if(mode!=0 && a_failed!=1)
	} else {//spur>0, just compute t and errors
		//compute t
		for (i=0; i<3; i++) {
			t[i] = ((yc[i] - u[i][0]*xc[0]) - u[i][1]*xc[1]) - u[i][2]*xc[2];
		}
	} //else spur>0 

	//compute rmsd
	for(i=0; i<3; i++){
		if( e[i] < 0 ) e[i] = 0;
		e[i] = sqrt( e[i] );           
	}            
	d = e[2];
	if( sigma < 0.0 ){
		d = - d;
	}
	d = (d + e[1]) + e[0];
	rms1 = (e0 - d) - d; 
	if( rms1 < 0.0 ) rms1 = 0.0;  

	rmsd=sqrt(rms1/(double)n);

	return rmsd;
}

int findMostOvlpChainSet(const int & minGrpIdx, const int & curChainIdx, vector<vector<Chain*> > & chainGrpsRefined, map<Chain*, vector<Chain*> > & chainContacts, vector<Chain*> & curChainSet)
{
	curChainSet.clear();

	int i, j, maxi, maxj, maxNumOvlpChains, numGrpsLeft = chainGrpsRefined.size(), numOvlpChains, ttlNumChains;
	vector<bool> chosen(chainGrpsRefined.size(), false);
	curChainSet.push_back(chainGrpsRefined[minGrpIdx][curChainIdx]);
	chosen[minGrpIdx] = true;
	numGrpsLeft--;
	ttlNumChains = chainContacts[chainGrpsRefined[minGrpIdx][curChainIdx]].size() + 1; // 1 for central chain

	while (numGrpsLeft > 0)
	{
		maxNumOvlpChains = -1;
		maxi = -1;
		maxj = -1;
		for (i = 0; i < chainGrpsRefined.size(); i++)
		{
			if (chosen[i])
			{
				continue;
			}

			for (j = 0; j < chainGrpsRefined[i].size(); j++)
			{
				numOvlpChains = calcNumOvlpChains(chainGrpsRefined[i][j], curChainSet, chainContacts);
				if (numOvlpChains > maxNumOvlpChains)
				{
					maxNumOvlpChains = numOvlpChains;
					maxi = i;
					maxj = j;
				}
			}
		}

		curChainSet.push_back(chainGrpsRefined[maxi][maxj]);
		chosen[maxi] = true;
		numGrpsLeft--;
		ttlNumChains += chainContacts[chainGrpsRefined[maxi][maxj]].size() + 1; // 1 for central chain
		ttlNumChains -= maxNumOvlpChains; // overlapping chains have been counted twice
	}

	return ttlNumChains;
}

map<Chain*, AtomPointerVector> getChainCa(Structure & S)
{
	map<Chain*, AtomPointerVector> chainCa;
	int i, j;
	for (i = 0; i < S.chainSize(); i++)
	{
		Chain& chn = S.getChain(i);
		for (j = 0; j < chn.residueSize(); j++)
		{
			Residue& r = chn.getResidue(j);
			Atom* ca = r.findAtom("CA", true);
		  	chainCa[&chn].push_back(ca);
		}
	}
	return chainCa;
}

map<Chain*, vector<double> > getChainCaCentroid(map<Chain*, AtomPointerVector> & chainCa)
{
	map<Chain*, vector<double> > chainCaCentroid;
	map<Chain*, AtomPointerVector>::iterator it;
	double x, y, z;
	int i;	
	for (it = chainCa.begin(); it != chainCa.end(); ++it)
	{
		x = 0.0;
		y = 0.0;
		z = 0.0;
		for (i = 0; i < it->second.size(); i++)
		{
			x += it->second[i]->getX();
			y += it->second[i]->getY();
			z += it->second[i]->getZ();
		}
		x = x / (it->second.size());
		chainCaCentroid[it->first].push_back(x);
		y = y / (it->second.size());
		chainCaCentroid[it->first].push_back(y);
		z = z / (it->second.size());
		chainCaCentroid[it->first].push_back(z);
	}

	return chainCaCentroid;
}

map<Chain*, vector<Chain*> > getChainContacts(Structure& S) {
	map<Chain*, vector<Chain*> > chainContacts;
	int i, j;
	for (i = 0; i < S.chainSize(); i++) {
		AtomPointerVector atomsA = S.getChain(i).getAtoms();
		if (chainContacts.find(&(S.getChain(i))) == chainContacts.end()) {
			// deal with the situation that there is only one chain
			chainContacts[&(S.getChain(i))] = vector<Chain*>();
		}
		for (j = i + 1; j < S.chainSize(); j++) {
			AtomPointerVector atomsB = S.getChain(j).getAtoms();
			if (areUnitsInContact(atomsA, atomsB)) {
				chainContacts[&(S.getChain(i))].push_back(&(S.getChain(j)));
				chainContacts[&(S.getChain(j))].push_back(&(S.getChain(i)));
			}
		}
	}
	return chainContacts;
}

map<Chain*, vector<string> > getChainSeq(Structure & S)
{
	map<Chain*, vector<string> > chainSeq;
	int i, j;
	for (i = 0; i < S.chainSize(); i++)
	{
		Chain& chn = S.getChain(i);
		for (j = 0; j < chn.positionSize(); j++)
		{
			chainSeq[&chn].push_back((chn.getResidue(j)).getName());
		}
	}
	return chainSeq;
}

vector<CorrChains> getCorrChains(vector<Chain*> & contactA, vector<Chain*> & contactB, map<Chain*, vector<double> > & chnCaCentroid, const vector<double> & t, const vector<vector<double> > & u)
{
	vector<CorrChains> coChains;
	set<CorrChains, centDistComp> ccEns;
	set<CorrChains, centDistComp>::iterator it;
	double xa, ya, za, xb, yb, zb;
	int i, j;
	CorrChains cc;
	for (i = 0; i < contactA.size(); i++)
	{
		xa = t[0] + u[0][0] * chnCaCentroid[contactA[i]][0] + u[0][1] * chnCaCentroid[contactA[i]][1] + u[0][2] * chnCaCentroid[contactA[i]][2];
		ya = t[1] + u[1][0] * chnCaCentroid[contactA[i]][0] + u[1][1] * chnCaCentroid[contactA[i]][1] + u[1][2] * chnCaCentroid[contactA[i]][2];
		za = t[2] + u[2][0] * chnCaCentroid[contactA[i]][0] + u[2][1] * chnCaCentroid[contactA[i]][1] + u[2][2] * chnCaCentroid[contactA[i]][2];
		for (j = 0; j < contactB.size(); j++)
		{
			xb = chnCaCentroid[contactB[j]][0];
			yb = chnCaCentroid[contactB[j]][1];
			zb = chnCaCentroid[contactB[j]][2];

			cc.chnA = contactA[i];
			cc.chnB = contactB[j];
			cc.centDist = sqrt((xa - xb) * (xa - xb) + (ya - yb) * (ya - yb) + (za - zb) * (za - zb));
			ccEns.insert(cc);
		}
	}
	while (!ccEns.empty())
	{
		it = ccEns.begin();
		coChains.push_back(*it);
		for (it = ccEns.begin(); it != ccEns.end(); ++it)
		{
			if ((it->chnA == coChains[coChains.size() - 1].chnA) || (it->chnB == coChains[coChains.size() - 1].chnB))
			{
				ccEns.erase(it);
			}
		}
	}	

	return coChains;
}

void gridPoint(Atom* a, CartesianPoint& c, double gs, int* ip, int* jp, int* kp) {
  CartesianPoint d = (CartesianPoint) a->getCoor() - c;
  int i = int(d.getX() < 0 ? ceil(d.getX()/gs): floor(d.getX()/gs));
  int j = int(d.getY() < 0 ? ceil(d.getY()/gs): floor(d.getY()/gs));
  int k = int(d.getZ() < 0 ? ceil(d.getZ()/gs): floor(d.getZ()/gs));
  if (ip != NULL) *ip = i;
  if (jp != NULL) *jp = j;
  if (kp != NULL) *kp = k;
}

void mergeChainGroups(Chain* chnA, Chain* chnB, map<Chain*, int> & chainGrpIdx, vector<vector<Chain*> > & chainGrps)
{
	if (chainGrpIdx[chnA] == chainGrpIdx[chnB])
	{
		return;
	}
	ASSERT((chainGrps[chainGrpIdx[chnA]].size() > 0) && (chainGrps[chainGrpIdx[chnB]].size() > 0), "number of chains in a group should be > 0");
	int i, idxA = chainGrpIdx[chnA], idxB = chainGrpIdx[chnB];
	vector<Chain*>::iterator it, beg, end;
	if (chainGrps[idxA].size() < chainGrps[idxB].size())
	{		
		it = chainGrps[idxB].end();
		beg = chainGrps[idxA].begin();
		end = chainGrps[idxA].end();
		chainGrps[idxB].insert(it, beg, end);
		for (i = 0; i < chainGrps[idxA].size(); i++)
		{
			chainGrpIdx[chainGrps[idxA][i]] = idxB;
		}
		chainGrps[idxA].clear();
	}
	else
	{		
		it = chainGrps[idxA].end();
		beg = chainGrps[idxB].begin();
		end = chainGrps[idxB].end();
		chainGrps[idxA].insert(it, beg, end);
		for (i = 0; i < chainGrps[idxB].size(); i++)
		{
			chainGrpIdx[chainGrps[idxB][i]] = idxA;
		}
		chainGrps[idxB].clear();
	}
}
