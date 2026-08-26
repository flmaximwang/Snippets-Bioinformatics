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

#if !defined(REDUNDANCY_H)
#define REDUNDANCY_H

#include "Common.h"

class CorrChains;

bool alignChains(const vector<string> &, const vector<string> &, const vector<string> &, const double &, AtomPointerVector &, AtomPointerVector &, const int &, const double &, vector<double> &, vector<vector<double> > &, const double &, const vector<double> &, const vector<vector<double> > &);
bool alignSeqNeedlemanWunsch(const vector<string> &, const vector<string> &, const vector<string> &, const double &, vector<string> &, vector<string> &);
bool areUnitsInContact(AtomPointerVector*, AtomPointerVector*);
int calcNumOvlpChains(Chain*, vector<Chain*> &, map<Chain*, vector<Chain*> > &);
double calcRmsdKabsch(AtomPointerVector &, AtomPointerVector &, int, vector<double> &, vector<vector<double> > &);
int findMostOvlpChainSet(const int &, const int &, vector<vector<Chain*> > &, map<Chain*, vector<Chain*> > &, vector<Chain*> &);
map<Chain*, AtomPointerVector> getChainCa(Structure &);
map<Chain*, vector<double> > getChainCaCentroid(map<Chain*, AtomPointerVector> &);
map<Chain*, vector<Chain*> > getChainContacts(Structure &);
map<Chain*, vector<string> > getChainSeq(Structure &);
vector<CorrChains> getCorrChains(vector<Chain*> &, vector<Chain*> &, map<Chain*, vector<double> > &, const vector<double> &, const vector<vector<double> > &);
void gridPoint(Atom* a, CartesianPoint& c, double gs, int* ip=NULL, int* jp=NULL, int* kp=NULL);
void mergeChainGroups(Chain*, Chain*, map<Chain*, int> &, vector<vector<Chain*> > &);

class CorrChains
{
	public:
		Chain* chnA;
		Chain* chnB;
		double centDist;

		CorrChains() {}
		~CorrChains() {}
};

struct centDistComp
{
	bool operator() (const CorrChains & ccA, const CorrChains & ccB)
	{
		if (ccA.centDist != ccB.centDist)
		{
			return (ccA.centDist < ccB.centDist);
		}
		else
		{
			return ((&ccA) < (&ccB));
		}
	}
};

#endif
