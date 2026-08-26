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

#include "SearchOptions.h"

SearchOptions::SearchOptions()
{
	_bbrmsd = false;
	_ddzscore = false;
	_dmode = 0;
	_keeporder = false;
	_minn = 0;
	_skiprmsd = false;
	_phieps = 180.0;
	_psieps = 180.0;
	_rmode = 0;
	_rthresh = MAX_DOUBLE;
	_otype = "match";
	_topn = 0;
	_tune = 0.5;
}

void SearchOptions::setDistDevCut(const char * s)
{
	if (!((sscanf(s, "%lf", &(_deps)) == 1) && (_deps >= 0.0)))
	{
		usage();
		error("bad greedy distance deviation cutoff value");
	}
	_dmode = 1;
}

void SearchOptions::setGapLen(const string & str)
{
	ASSERT(_gaplen.size() == 0, "gap length constraints have been specified");
	
	string subs, subss;
	size_t cur, next, subc, subn;

	// the first gap length constraint
	_gaplen.push_back(vector<int>());
	// looking for the first delimiter (i.e., semicolon)
	cur = 0;
	next = str.find_first_of(";", cur);
	while (next != string::npos)
	{
		// found a delimiter, indicating the end of the current gap length constraint		
		subs = str.substr(cur, next - cur);
		if (subs.length() > 0)
		{
			// looking for the connector (i.e., minus)
			subc = 0;
			subn = subs.find_first_of("-", subc);
			while (subn != string::npos)
			{
				// found a connector, indicating the end of the current component of the current gap length constraint
				subss = subs.substr(subc, subn - subc);
				if (!isDigit(subss))
				{
					usage();
					error("bad loop length value");
				}
				// found a valid component of the current gap length constraint
				_gaplen[_gaplen.size() - 1].push_back(atoi(subss.c_str()));
				// looking for the next connector, though only one connector is allowed here
				subc = subn + 1;
				subn = subs.find_first_of("-", subc);				
			}
			// found no connector
			if (subc < subs.length())
			{
				// found the last component of the current gap length constraint
				subss = subs.substr(subc);
				if (!isDigit(subss))
				{
					usage();
					error("bad loop length value");
				}
				// found a valid component of the current gap length constraint
				_gaplen[_gaplen.size() - 1].push_back(atoi(subss.c_str()));
			}
			// a connector should be followed by a vaild component
			if ('-' == subs[subs.length() - 1])
			{
				usage();
				error("bad loop length value");
			}
			// at most two valid components are allowed here
			if (_gaplen[_gaplen.size() - 1].size() > 2)
			{
				usage();
				error("bad loop length value");
			}
		}
		// a delimiter should be followed by the next gap length constraint		
		_gaplen.push_back(vector<int>());
		// looking for the next delimiter		
		cur = next + 1;
		next = str.find_first_of(";", cur);		
	}
	// found no delimiter
	if (cur < str.length())
	{
		// found the last gap length constraint
		subs = str.substr(cur);
		// looking for the connector (i.e., minus)
		subc = 0;
		subn = subs.find_first_of("-", subc);
		while (subn != string::npos)
		{
			// found a connector, indicating the end of the current component of the current gap length constraint
			subss = subs.substr(subc, subn - subc);
			if (!isDigit(subss))
			{
				usage();
				error("bad loop length value");
			}
			// found a valid component of the current gap length constraint
			_gaplen[_gaplen.size() - 1].push_back(atoi(subss.c_str()));
			// looking for the next connector, though only one connector is allowed here
			subc = subn + 1;
			subn = subs.find_first_of("-", subc);				
		}
		// found no connector
		if (subc < subs.length())
		{
			// found the last component of the current gap length constraint
			subss = subs.substr(subc);
			if (!isDigit(subss))
			{
				usage();
				error("bad loop length value");
			}
			// found a valid component of the current gap length constraint
			_gaplen[_gaplen.size() - 1].push_back(atoi(subss.c_str()));
		}
		// a connector should be followed by a vaild component
		if ('-' == subs[subs.length() - 1])
		{
			usage();
			error("bad loop length value");
		}
		// at most two valid components are allowed here
		if (_gaplen[_gaplen.size() - 1].size() > 2)
		{
			usage();
			error("bad loop length value");
		}
	}
#if defined(DEBUG_GAPLEN)
	int i, j;
	for (i = 0; i < _gaplen.size(); i++)
	{
		for (j = 0; j < _gaplen[i].size(); j++)
		{
			cout << " " << _gaplen[i][j];
		}
		cout << "\n";
	}
	exit(-1);
#endif
}

void SearchOptions::setMinN(const char * s)
{
	if (!((sscanf(s, "%d", &(_minn)) == 1) && (_minn >= 0)))
	{
		usage();
		error("bad minN value");
	}
}

void SearchOptions::setPhiDevCut(const char * s)
{
	if (!((sscanf(s, "%lf", &(_phieps)) == 1) && (_phieps >= 0.0)))
	{
		usage();
		error("bad Phi angle deviation cutoff value");
	}
}

void SearchOptions::setPsiDevCut(const char * s)
{
	if (!((sscanf(s, "%lf", &(_psieps)) == 1) && (_psieps >= 0.0)))
	{
		usage();
		error("bad Psi angle deviation cutoff value");
	}
}

void SearchOptions::setRmsdBoundMode(const char * s)
{
	if (!((sscanf(s, "%d", &(_rmode)) == 1) && ((_rmode == 0) || (_rmode == 1) || (_rmode == 2))))
	{
		usage();
		error("bad CA RMSD bounding mode value");
	}
}

void SearchOptions::setRmsdCut(const char * s)
{
	if (!((sscanf(s, "%lf", &(_rthresh)) == 1) && (_rthresh > 0.0)))
	{
		usage();
		error("bad CA RMSD cutoff value");
	}
}

void SearchOptions::setStructOutDir(const string & dir)
{
	_sodname = dir;
	if ((!_sodname.empty()) && (!exist(_sodname.c_str())))
  	{
  		ASSERT(mkdir(_sodname.c_str(), 0755) == 0, "could not make directory %s", _sodname.c_str());
  	}
}

void SearchOptions::setOutType(const string & type)
{
	_otype = type;
	if ((_otype.compare("full") != 0) && (_otype.compare("match") != 0) && (_otype.compare("wgap") != 0))
	{
		usage();
		error("bad output type value");
	}
}

void SearchOptions::setTopN(const char * s)
{
	if (!((sscanf(s, "%d", &(_topn)) == 1) && (_topn >= 0)))
	{
		usage();
		error("bad topN value");
	}
}

void SearchOptions::setTsFile(const string & fn)
{
	ASSERT(_tsfnames.size() == 0, "targets already exist");
	_tsfnames.push_back(fn);
}

void SearchOptions::setTsFiles(const string & list)
{
	ASSERT(_tsfnames.size() == 0, "targets already exist");
	file2array(list, _tsfnames);
}

void SearchOptions::setTuningParam(const char * s)
{
	if (!((sscanf(s, "%lf", &(_tune)) == 1) && (_tune >= 0.0) && (_tune <= 1.0)))
	{
		usage();
		error("bad tuning parameter value");
	}
}

void SearchOptions::usage()
{
	int w = 100, p1 = 1, p2 = p1 + 18;
	string description = 
	    "This program performs a structure-based search for "
	    "an arbitrary input structure (the query), composed of one or more disjoint segments, against a given structural database. Both the query and the "
	    "database need to be specified, the former via the --query option, and the latter via --target or --targetList options (see below). The search is "
	    "performed by backbone root-mean-square-deviation (RMSD), with the program returning provably all structural matches that are within the user-specified RMSD "
	    "cutoff (i.e., --rmsdCut option) of the query. Either CA-only or full-backbone (i.e., N, C, CA, and O atoms) RMSD can be used as the metric (see option "
	    "--bbRMSD). If you use this tool for your research, please cite:"
	    "\n\nZhou, G. Grigoryan, Protein Science, 24(4): 508-524, 2015";
	cout << endl << optionUsage("", suiteName(), w, 0, 0) << endl << endl;
    cout << optionUsage("", description, w, 0, 0) << "\n\n";
    cout << optionUsage("Query and database specification options:", "", w, 0, 0) << "\n";
	cout << optionUsage("--query", "required: query PDS file. To generate this from a PDB file, use the accompanying program createPDS.", w, p1, p2) << "\n";
	cout << optionUsage("--target", "target PDS file. If given, the query will be searched against just this one file. To search against a database, use the --targetList option below.", w, p1, p2) << "\n";
	cout << optionUsage("--targetList", "a file with a list of paths to PDS files (i.e., the database). Each target will be searched. Either --target or --targetList must be given.", w, p1, p2) << "\n";

    cout << optionUsage("\nSearch criteria:", "", w, 0, 0) << "\n";
	cout << optionUsage("--rmsdCut", "required: RMSD cutoff for defining a match (in Angstroms).", w, p1, p2) << "\n";
    cout << optionUsage("--bbRMSD", "a flag that changes the metric from the default CA-based RMSD to full-backbone RMSD.", w, p1, p2) << "\n";
    cout << optionUsage("--topN", "keep the best this many matches in terms of the search metric (must be integer); default is 0 (no limit). Can be combined with --rmsdCut to given useful behaviors. E.g., --rmsdCut 1.5 --topN 5000 will find all matches within 1.5 A of the query, unless there are over 5,000, in which case only the top 5,000 (by RMSD) will be found and reported. Note, --topN can improve performance, because it effective tells the search to lower the RMSD threshold once the desired number of matches are found.", w, p1, p2) << "\n";
	cout << optionUsage("--minN", "return at least this many of the best matches (in terms of the search metric), regardless of the specified RMSD cutoff. Default (i.e., 0) is to have no minimum. Must be a non-negative integer. If both --minN and --topN are specified, the former must be no greater than the latter.", w, p1, p2) << "\n";
	cout << optionUsage("--gapLen", "constrains matches by the length of sequence that maps between adjacent segments in the query (i.e., gap lengths). For example, "
	                                "suppose a query consists of two disjoint segments, A and B, which align onto regions A' and B' in a particular "
	                                "match M. Suppose further that A' and B' belong to the same chain in M. Then, the region in M after A' and before B' maps "
	                                "to the gap between A and B in the query. This option allows matches to be filtered by the number of residues in this gap region in M. "
	                                "Each gap restraint has the form 'min-max', where min and max are the minimal and maximal gap length, respectively. The number of "
	                                "these restraints must match the number of gaps in the query (i.e., number of disjoint segments in the query minus one), and gap "
	                                "restraints should be combined with a semicolon. E.g., for a three-segment query, one might say --gapLen '1-10;0-3', which would mean that the "
	                                "first gap (i.e., between segments 1 and 2 in the query) should be between 1 and 10 residues, and the second gap (i.e., "
	                                "between segments 2 and 3 in the query) should be between 0 and 3 residues. Another example: --gapLen '3;4-5' will require a length of 3 for "
	                                "the first gap and between 4 and 5 for the second. Important: with this option, the order of segments in the query is assumed "
	                                "to be as listed in the original PDB file of the query. Thus, specifying --gapLen stipulates the particular order of segments "
	                                "(in the N-to-C sense) in the matches. This is not the case in general, as MASTER will consider all orderings of disjoint segments "
	                                "when --gapLen is not specified.", w, p1, p2) << "\n";

    cout << optionUsage("", "\nOutput options (note, all output forms list matches in the same order):", w, 0, 0) << "\n";
    cout << optionUsage("--outType", "controls the output region of matches, affecting --seqOut and --structOut options below. Three "
                                     "output regions are possible: 'match' (default) -- just the matching region (i.e., just those residues that align onto the query), "
                                     "'full' -- the entire database entry, and 'wgap' -- the matching region plus any gap regions that map between disjoint segments in "
                                     "the query (only works when --gapLen is specified).", w, p1, p2) << "\n";
	cout << optionUsage("--matchOut", "file name for outputting match addresses. Each line contains all information about a single match (database entry, specific "
	                                  "alignment location, and RMSD). Although just a plain text file, it is mostly intended as a record of the result that can be "
	                                  "be used to later re-generate the matches in various ways (see reprocessing options below).", w, p1, p2) << "\n";
	cout << optionUsage("--seqOut", "file name for outputting match sequences (one per line); sequences of output regions are written.", w, p1, p2) << "\n";
	cout << optionUsage("--structOut", "name of directory for writing match structures in PDB format (one PDB file per match). The directory is created if it does "
	                                   "not exist. The region of structure output is controlled by --outType, but structures are always written optimally "
	                                   "superimposed onto the query over the matching region (by the search metric).", w, p1, p2) << "\n";

    cout << optionUsage("\nReprocessing of a previous search:", "", w, 0, 0) << "\n";
	cout << optionUsage("--matchIn","a file produced by --matchOut in a previous run. If specified with --query, will skip the search and instead will work with the "
	                                "given matches. Note that the same search database as in the original run should be specified (via either --target or --targetList), "
	                                "as the match file does not store the database but only match addresses.", w, p1, p2) << "\n";
	cout << optionUsage("--skipRMSD","optional: if specified with --matchIn, will produce outputs using the same RMSD metric as the input, skipping RMSD re-calculation, which would normally be the default behavior.", w, p1, p2) << "\n";
	cout << optionUsage("--keepOrder","optional: if specified with --matchIn, will produce outputs in the same order as the input and therefore skip match resorting, which would normally be the default behavior.", w, p1, p2) << "\n";

    cout << optionUsage("", "\nHeuristic search options. Using these is NOT recommended and will sacrifice the provable completeness of search results for (potential) search speed-up:", w, 0, 0) << "\n";
    cout << optionUsage("--rmsdMode", "RMSD bounding mode. 0 -- provable RMSD bounds will be calculated, guaranteeing that all matches within --rmsdCut will be found "
                                      "(default). 1 -- greedy bound that enforces some uniformity of RMSD residuals (see paper). 2 -- uses --rmsdCut for both the "
                                      "overall RMSD cutoff as well as the cutoff for partial matches (see paper).", w, p1, p2) << "\n";
    cout << optionUsage("--tune", "tuning parameter for greedy RMSD cutoff (i.e., when --rmsdMode is 1); 0.5 by default (see paper).", w, p1, p2) << "\n";
    cout << optionUsage("--dEps", "user-defined greedy distance deviation cutoff (in Angstroms). If given, rather than applying a provable bound on inter-segment "
                                  "distances, this cutoff will be greedily applied as a first filter.", w, p1, p2) << "\n";
    cout << optionUsage("--phiEps", "phi angle deviation cutoff (in degrees); default is 180.0, meaning no cutoff is applied.", w, p1, p2) << "\n";
    cout << optionUsage("--psiEps", "psi angle deviation cutoff (in degrees); default is 180.0, meaning no cutoff is applied.", w, p1, p2) << "\n";
    cout << optionUsage("--ddZscore", "output a Z-score that describes the distribution of inter-segment distance deviations (between query and matches) relative to "
                                      "the greedy cutoff --dEps. High Z-scores (> 3.5, in our experience) indicate a good choice of --dEps and suggests that all or "
                                      "nearly all matches were found despite the greedy constraint.", w, p1, p2) << "\n";

//	cout << optionUsage("--ddOutFile", "optional: file name to output distance deviations of matches to (only works when Z score is output). By default, no distance deviations are output. ", w, p1, p2) << "\n";
    cout << endl;
}
