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

#if !defined(SEARCHRESULTS_H)
#define SEARCHRESULTS_H

#include "Common.h"
#include "QueryStruct.h"

class Match;

class SearchResults
{
	public:
		SearchResults(const int & most = 0, const int & least = 0);
		SearchResults(QueryStruct& qs, const int & most = 0, const int & least = 0);
		SearchResults(const string & mifname, QueryStruct & qs, const vector<vector<int> > & gaplen = vector<vector<int> >(), const bool & ko = false, const int & most = 0, const int & least = 0, const double & rmsdcut = MAX_DOUBLE);
		~SearchResults();

		int addTsFile(const string &); // NOTE: (low priority) may want to make this function protected and make Match a friend of SearchResults
		void clear();
		bool empty() const { return _L.empty(); }
		bool full();
		int getMatchLeast() const { return _least; }
		int getMatchMost() const { return _most; }
		unsigned long int getNumMatInst() const { return _nMatInst; }
		int getNumTs() const { return _tsfiles.size(); }
		QueryStruct& getQs() { return (*_qs); }
		int getQsNumSeg() const { return _qs->getNumSeg(); }
		int getQsSegLen(const int & i) const { return (_qs->getResBefBrk(i + 1) - _qs->getResBefBrk(i)); }
		vector<char*>& getSeq() { return _seq; }
		char* getSeq(const int & i) { return _seq[i]; }
		string getTsFile(const int & i) const { return *(_tsfiles[i]); }
		int getTsNumRes(const int & i) const { return _tsnr[i]; }
		void initTsNumRes() { _tsnr.assign(_tsfiles.size(), -1); }
		void insertMatch(Match &);
		set<Match>& matchSet() { return _L; } // a getter for a member object
		vector<Match*> matchVector(); // utility function, in case the user wants a vector
		friend ostream& operator<<(ostream &, SearchResults &);
		void removeWorstMatch();
		void setMatchLeast(const int & least) { _least = least; }
		void setMatchList(const string & mifname, const vector<vector<int> > & gaplen = vector<vector<int> >(), const bool & ko = false, const double & rmsdcut = MAX_DOUBLE);
		void setMatchMost(const int & most) { _most = most; }
		void setQs(const string & fn);
		void setQs(QueryStruct & qs) { _qs = &qs; }
		void setTsFiles(const vector<string> &);
		void setTsNumRes(const int &, const int &);
		unsigned int size() const { return _L.size(); }	
		double worstRmsd() const;

	private:
		set<Match> _L;
		int _least; // at least this many matches
		int _most; // at most this many matches
		unsigned long int _nMatInst; // number of matches inserted so far, including those removed due to the 'most' limit
		QueryStruct* _qs;
		bool _qsfree;
		vector<char*> _seq;
		vector<string*> _tsfiles;
		map<string, int> _tsfmap;
		vector<int> _tsnr;
};

#endif
