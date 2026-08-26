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
// POST V1.1
#include "Common.h"
#include "Search.h"
#include "SearchOptions.h"
#include "SearchResults.h"

void parseCmdLine(int, char**, SearchOptions &);

int main(int argc, char* argv[])
{
	double elapse;
	const long int microPerSec = 1000000L;
	vector<pair<Match*, int> > toSort;
	SearchResults mlist;
	SearchOptions sopts;
	struct timeval start, end;
	clock_t t;

	parseCmdLine(argc, argv, sopts);

	// initialize match list
	mlist.setMatchLeast(sopts.getMinN());
	mlist.setMatchMost(sopts.getTopN());
	mlist.setQs(sopts.getQsFile());	
	if (sopts.getGapLen().size() > 0)
	{
		ASSERT(mlist.getQs().getNumSeg() == (sopts.getGapLen().size() + 1), "numbers of query segments and gap length constraints not compatible");
	}	

	if (!sopts.getMatchInFile().empty())
	{
		ASSERT(gettimeofday(&start, NULL) == 0, "could not obtain the parse start time");
		t = clock();

		if (sopts.getBBRmsd())
		{
			mlist.getQs().readBBCoords();
			ASSERT(mlist.getQs().fullBB(), "backbone atoms missing in file %s", mlist.getQs().getFileName().c_str());
		}
		mlist.setMatchList(sopts.getMatchInFile(), sopts.getGapLen(), sopts.getKeepOrder(), sopts.getRmsdCut());
		
		t = clock() - t;
		ASSERT(gettimeofday(&end, NULL) == 0, "could not obtain the parse finish time");
		elapse = double((end.tv_sec - start.tv_sec) * microPerSec + end.tv_usec - start.tv_usec) / double(microPerSec);
		printf("Parse completed in %.2f seconds, and used %.2f seconds of CPU time.\n", elapse, float(t) / CLOCKS_PER_SEC);
	}
	else
	{
		ASSERT(gettimeofday(&start, NULL) == 0, "could not obtain the search start time");
		t = clock();

		mlist.setTsFiles(sopts.getTsFiles());
		// NOTE: no need to pass all entries of sopts separately
		searchByDistDistr(mlist, sopts.getBBRmsd(), sopts.getPhiDevCut(), sopts.getPsiDevCut(), sopts.getDistDevBoundMode(), sopts.getDistDevCut(), sopts.getRmsdBoundMode(), sopts.getRmsdCut(), sopts.getTuningParam(), sopts.getGapLen());

		t = clock() - t;
		ASSERT(gettimeofday(&end, NULL) == 0, "could not obtain the search finish time");
		elapse = double((end.tv_sec - start.tv_sec) * microPerSec + end.tv_usec - start.tv_usec) / double(microPerSec);
		printf("Search completed in %.2f seconds, and used %.2f seconds of CPU time.\n", elapse, float(t) / CLOCKS_PER_SEC);
	}
	ASSERT(gettimeofday(&start, NULL) == 0, "could not obtain the output start time");
	t = clock();

	beforeOut(mlist, toSort, sopts.getMatchInFile(), sopts.getBBRmsd(), sopts.getSeqOutFile(), sopts.getStructOutDir(),  sopts.getOutType(), (sopts.getDistDevZscore() && (mlist.getQs().getNumSeg() > 1)), sopts.getGapLen(), sopts.getKeepOrder(), sopts.getSkipRMSD());
	outMatch(mlist, sopts.getMatchOutFile(), toSort);
	outSeq(mlist, sopts.getSeqOutFile(), sopts.getOutType(), sopts.getGapLen(), toSort);
	renameStruct(mlist, sopts.getStructOutDir(), sopts.getOutType(), toSort);
	outDistDevZscore(mlist, (sopts.getDistDevZscore() && (mlist.getQs().getNumSeg() > 1)), sopts.getDistDevCut());
	
	t = clock() - t;
	ASSERT(gettimeofday(&end, NULL) == 0, "could not obtain the output finish time");
	elapse = double((end.tv_sec - start.tv_sec) * microPerSec + end.tv_usec - start.tv_usec) / double(microPerSec);
	printf("Output completed in %.2f seconds, and used %.2f seconds of CPU time.\n", elapse, float(t) / CLOCKS_PER_SEC);
	
	return 0;
}

void parseCmdLine(int argc, char** argv, SearchOptions & sopts)
{
	map<string, bool> spec;

	while (1)
	{
		int oind = 0;
		static struct option opts[] =
		{
			{"query", 1, 0, 1},
			{"target", 1, 0, 2},
			{"targetList", 1, 0, 3},
      		{"topN", 1, 0, 4},
      		{"phiEps", 1, 0, 5},
      		{"psiEps", 1, 0, 6},
      		{"dEps", 1, 0, 7},
      		{"rmsdCut", 1, 0, 8},
      		{"structOut", 1, 0, 9},
      		{"matchOut", 1, 0, 10},
      		{"rmsdMode", 1, 0, 11},
      		{"ddZscore", 0, 0, 12},
      		{"tune", 1, 0, 13},
      		{"seqOut", 1, 0, 14},
			{"bbRMSD", 0, 0, 15},
			{"outType", 1, 0, 16},
			{"matchIn", 1, 0, 17},
			{"gapLen", 1, 0, 18},
			{"skipRMSD", 0, 0, 19},
			{"keepOrder", 0, 0, 20},
			{"minN", 1, 0, 21},
      		{0, 0, 0, 0}
		};

		int c = getopt_long (argc, argv, "", opts, &oind);
		if (c == -1)
		{
			break;
		}

		switch (c) {
			case 1:
				sopts.setQsFile(string(optarg));
        		spec[string(opts[oind].name)] = true;
        		break;

      		case 2:
				sopts.setTsFile(string(optarg));
        		spec[string(opts[oind].name)] = true;
        		break;

      		case 3:
				sopts.setTsFiles(string(optarg));
        		spec[string(opts[oind].name)] = true;
        		break;

      		case 4:
				sopts.setTopN(optarg);
        		spec[string(opts[oind].name)] = true;
        		break;

	  		case 5:
				sopts.setPhiDevCut(optarg);
	  			spec[string(opts[oind].name)] = true;
				break;

	  		case 6:
				sopts.setPsiDevCut(optarg);
	  			spec[string(opts[oind].name)] = true;
				break;

	  		case 7:
				sopts.setDistDevCut(optarg);
	  			spec[string(opts[oind].name)] = true;
				break;

      		case 8:
				sopts.setRmsdCut(optarg);
        		spec[string(opts[oind].name)] = true;
        		break;

	  		case 9:
				sopts.setStructOutDir(string(optarg));
				spec[string(opts[oind].name)] = true;
	  			break;

	  		case 10:
				sopts.setMatchOutFile(string(optarg));
				spec[string(opts[oind].name)] = true;
	  			break;
			
	  		case 11:
				sopts.setRmsdBoundMode(optarg);
        		spec[string(opts[oind].name)] = true;
        		break;

	  		case 12:
				sopts.setDistDevZscore(true);
        		spec[string(opts[oind].name)] = true;
        		break;
			
	  		case 13:
				sopts.setTuningParam(optarg);
				spec[string(opts[oind].name)] = true;
		        break;
			
	  		case 14:
				sopts.setSeqOutFile(string(optarg));
				spec[string(opts[oind].name)] = true;
			  	break;

			case 15:
				sopts.setBBRmsd(true);
				spec[string(opts[oind].name)] = true;
			  	break;

			case 16:
				sopts.setOutType(string(optarg));
				spec[string(opts[oind].name)] = true;
			  	break;

			case 17:
				sopts.setMatchInFile(string(optarg));
				spec[string(opts[oind].name)] = true;
	  			break;

			case 18:
				sopts.setGapLen(string(optarg));
				spec[string(opts[oind].name)] = true;
	  			break;

			case 19:
				sopts.setSkipRMSD(true);
				spec[string(opts[oind].name)] = true;
	  			break;

			case 20:
				sopts.setKeepOrder(true);
				spec[string(opts[oind].name)] = true;
	  			break;

			case 21:
				sopts.setMinN(optarg);
				spec[string(opts[oind].name)] = true;
				break;

      		case '?':
				sopts.usage();
				exit(-1);

      		default:
        		printf ("?? getopt returned character code %d ??\n", c);
				sopts.usage();
				exit(-1);
		}
	}  	

	if (optind < argc)
	{
		sopts.usage();
		printf ("non-option ARGV-elements: ");
		while (optind < argc)
    	{
			printf ("%s ", argv[optind++]);
    	}
		printf ("\n");
		exit(-1);
	}

	// make sure all spec options have been specified
  	if (!((spec.find(string("query")) != spec.end())
		&& ((((spec.find(string("target")) != spec.end()) || (spec.find(string("targetList")) != spec.end())) && (spec.find(string("rmsdCut")) != spec.end()))
		|| (spec.find(string("matchIn")) != spec.end()))))
  	{
		sopts.usage();
		error("not all required options specified");
	}
	if ((spec.find(string("target")) != spec.end()) && (spec.find(string("targetList")) != spec.end()))
	{
		sopts.usage();
		error("'--target' can not be specified with '--targetList'");
	}
	if (((spec.find(string("target")) != spec.end()) || (spec.find(string("targetList")) != spec.end())) 
		&& (spec.find(string("matchIn")) != spec.end()))
	{
		sopts.usage();
		error("'--matchIn' can not be specified with either '--target' or '--targetList'");
	}

	if ((sopts.getTopN() > 0) && (sopts.getTopN() < sopts.getMinN()))
	{
		sopts.usage();
		error("'--topN' should have been >= '--minN'");
	}

	if ((sopts.getOutType().compare("wgap") == 0) && (sopts.getGapLen().size() <= 0))
	{
		sopts.usage();
		error("'wgap' output type only works when '--gapLen' option has been specified");
	}

	if ((spec.find(string("matchIn")) != spec.end()) && (spec.find(string("minN")) != spec.end()) && (sopts.getMinN() > 0))
	{
		warning("'--minN' does not work with '--matchIn'");
	}
	
#if !defined(DEBUG_ZSCORE)
	if (0 == sopts.getDistDevBoundMode())
	{
		sopts.setDistDevZscore(false);
	}
#endif

	if ((!sopts.getMatchInFile().empty()) && (!sopts.getStructOutDir().empty()))
	{
		sopts.setSkipRMSD(false);
	}

#if defined(DEBUG_CMDLINE)
	cout << "output backbone RMSDs? " << sopts.getBBRmsd() << "\n";
	cout << "distance deviation cutoff: " << sopts.getDistDevCut() << "\n";
	cout << "distance deviation mode: " << sopts.getDistDevBoundMode() << "\n";
	cout << "output distance deviation Z-score? " << sopts.getDistDevZscore() << "\n";
	cout << "minN: " << sopts.getMinN() << "\n";
	cout << "Phi angle deviation cutoff: " << sopts.getPhiDevCut() << "\n";
	cout << "Psi angle deviation cutoff: " << sopts.getPsiDevCut() << "\n";
	cout << "CA RMSD mode: " << sopts.getRmsdBoundMode() << "\n";
	cout << "CA RMSD cutoff: " << sopts.getRmsdCut() << "\n";
	cout << "output type: " << sopts.getOutType() << "\n";
	cout << "topN: " << sopts.getTopN() << "\n";
	cout << "tuning parameter for CA RMSD bounding: " << sopts.getTuningParam() << "\n";
	exit(-1);
#endif
}
