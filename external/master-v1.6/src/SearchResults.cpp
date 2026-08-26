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

#include "Match.h"
#include "SearchResults.h"

SearchResults::SearchResults(const int & most, const int & least)
{
	_least = least;
	_most = most;
	_nMatInst = 0;
	_qs = NULL;
	_qsfree = false;
}

SearchResults::SearchResults(QueryStruct& qs, const int & most, const int & least)
{
	_least = least;
	_most = most;
	_nMatInst = 0;
	_qs = &qs;
	_qsfree = false;
}

SearchResults::SearchResults(const string & mifname, QueryStruct & qs, const vector<vector<int> > & gaplen, const bool & ko, const int & most, const int & least, const double & rmsdcut)
{
	_least = least;
	_most = most;
	_nMatInst = 0;
	_qs = &qs;
	_qsfree = false;
	setMatchList(mifname, gaplen, ko, rmsdcut);
}

SearchResults::~SearchResults()
{
	clear();
}

int SearchResults::addTsFile(const string & fn)
{
	if (_tsfmap.find(fn) == _tsfmap.end())
	{
		int i = _tsfmap.size();
		_tsfmap[fn] = i;
		_tsfiles.push_back((string *) &(_tsfmap.find(fn)->first));
	}
	return _tsfmap[fn];
}

void SearchResults::clear()
{
	_L.clear();
	_least = 0;
	_most = 0;
	_nMatInst = 0;
	_qs->closeProteinStructFile();
	if (_qsfree)
	{
		delete _qs;
	}
	for (int i = 0; i < _seq.size(); i++)
	{
		delete [] _seq[i];
	}
	_seq.clear();
	_tsfiles.clear();
	_tsfmap.clear();
	_tsnr.clear();
}

bool SearchResults::full()
{
	if ((_most > 0) && (_L.size() >= _most))
	{
		return true;
	}
	return false;
}

void SearchResults::insertMatch(Match & m)
{
	m.setInsertTime(_nMatInst);
	_L.insert(m);
	_nMatInst++;
	if ((_most > 0) && (_L.size() > _most))
	{
		_L.erase(--_L.end());
	}
}

vector<Match*> SearchResults::matchVector()
{
	vector<Match*> matchVec(_L.size(), NULL);
	set<Match>::iterator it;
	int i = 0;
	for (it = _L.begin(); it != _L.end(); ++it)
	{
		matchVec[i++] = (Match *) &(*it);
	}
	return matchVec;
}

ostream& operator<<(ostream & os, SearchResults & L)
{
	set<Match>::iterator it;
	for (it = L.matchSet().begin(); it != L.matchSet().end(); ++it)
	{
		os << (*it) << "\n";
	}

	return os;
}

void SearchResults::removeWorstMatch()
{
	_L.erase(--_L.end());
}

void SearchResults::setMatchList(const string & mifname, const vector<vector<int> > & gaplen, const bool & ko, const double & rmsdcut)
{
	ASSERT(_L.empty(), "match list already exists");
	
	FILE* ifp;
	int maxline = 1000;
	char *line = (char*) malloc(maxline * sizeof(char));
	Match m(this, -1, ko);

	openFileC(ifp, mifname, "r");

	while (fgets(line, maxline, ifp) != NULL)
	{
		ASSERT(line[strlen(line) - 1] == '\n', "lines in file %s are over %d long. Please increase max line limit and recompile", mifname.c_str(), maxline);
		if (strlen(line) > 0)
		{
			if (0 == m.parseMatch(line, gaplen))
			{
				if ((rmsdcut >= MAX_DOUBLE) || ((rmsdcut < MAX_DOUBLE) && (m.getRmsd() <= rmsdcut)))
				{
					insertMatch(m);
				}
			}
		}
		if (feof(ifp))
		{
			break;
		}
	}

	ASSERT(0 == ferror(ifp), "could not read file %s", mifname.c_str());	
	fclose(ifp);
	free(line);
}

void SearchResults::setQs(const string & fn)
{
	_qs = new QueryStruct;
	_qsfree = true;
	_qs->readQueryStruct(fn);
}

void SearchResults::setTsFiles(const vector < string > & list)
{
	int i;
	_tsfiles.assign(list.size(), NULL);
	for (i = 0; i < list.size(); i++)
	{
		_tsfiles[i] = (string *) &(list[i]);
	}
}

void SearchResults::setTsNumRes(const int & i, const int & nr)
{
	ASSERT(_tsnr[i] <= 0, "each target should set its number of residues only once");
	_tsnr[i] = nr;
}

double SearchResults::worstRmsd() const
{
	if (_L.empty())
	{
		return MAX_DOUBLE;
	}
	else
	{
		return _L.rbegin()->getRmsd();
	}
}

