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

#if !defined(SEARCH_H)
#define SEARCH_H

#include "CenResGrp.h"
#include "Common.h"
#include "Match.h"
#include "QueryStruct.h"
#include "SearchResults.h"
#include "TargetStruct.h"

void auxSearchByDistDistrBB(SearchResults &, CenResGrp &, QueryStruct &, TargetStruct &, const int &, const double &, const double &, const int &, const double &, const int &, const double &, const double &);
void auxSearchByDistDistrCA(SearchResults &, CenResGrp &, QueryStruct &, TargetStruct &, const int &, const double &, const double &, const int &, const double &, const int &, const double &, const double &);
void beforeOut(SearchResults & mlist, vector<pair<Match*, int> > & toSort, const string & mifname = "", const bool & bbrmsd = false, const string & sofname = "", const string & sodname = "", const string & otype = "", const bool & ddzscore = false, const vector<vector<int> > & gaplen = vector<vector<int> >(), const bool & keeporder = false, const bool & skiprmsd = false);
double calcMaxDistDev(const vector<vector<double> > &, const vector<int> &, const vector<vector<double> > &, const vector<int> &);
double calcRmsdKabsch(const vector<vector<double> > &, const vector<int> &, const vector<vector<double> > &, const vector<int> &, const int &, double*, double**);
double distDevCut(const int &, const double &, const int &, const double &, const int &, const double &);
void extractSeq(string &, const Match &, char*, const string &, const vector<vector<int> > &);
void masterSearch(QueryStruct & qs, TargetStruct & ts, SearchResults & mlist, const bool & bbrmsd, const double & phieps = 180.0, const double & psieps = 180.0, const int & dmode = 0, const double & deps = 2.0);
void outDistDevZscore(SearchResults &, const bool &, const double &);
void outMatch(SearchResults &, const string &, const vector<pair<Match*, int> > &);
void outSeq(SearchResults &, const string &, const string &, const vector<vector<int> > &, const vector<pair<Match*, int> > &);
void renameStruct(SearchResults &, const string &, const string &, vector<pair<Match*, int> > &);
double rmsdCut(const int &, const double &, const int &, const double &, const int &, const int &, const double &);
void searchByDistDistr(SearchResults &, const bool &, const double &, const double &, const int &, const double &, const int &, const double &, const double &, const vector<vector<int> > &);
void transformStructByRMSD(vector<string> &, const vector<string> &, const vector<int> &, double*, double**);

#endif
